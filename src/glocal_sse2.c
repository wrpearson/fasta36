/******************************************************************
  Copyright 2010 by Michael Farrar.  All rights reserved.
  This program may not be sold or incorporated into a commercial product,
  in whole or in part, without written consent of Michael Farrar.  For 
  further information regarding permission for use or reproduction, please 
  contact: Michael Farrar at farrar.michael@gmail.com.
*******************************************************************/

/*
  Written by Michael Farrar, 2010.
  Please send bug reports and/or suggestions to farrar.michael@gmail.com.
*/

#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
#include "param.h"
#include "dropgsw2.h"
#include "global_sse2.h"

#ifdef __SUNPRO_C
#include <sunmedia_intrin.h>
#else
#include <emmintrin.h>
#endif

#ifdef SW_SSE2

static inline __m128i
max_epu16(__m128i a, __m128i b)
{
  a = _mm_subs_epu16 (a, b);
  b = _mm_adds_epu16 (b, a);
  return b;
}

int
glocal_sse2_word(int                  queryLength,
                 unsigned short      *profile,
                 const unsigned char *dbSeq,
                 int                  dbLength,
                 unsigned short       gapOpen,
                 unsigned short       gapExtend,
                 unsigned short       ceiling,
                 struct f_struct     *f_str)
{
  int     i, j;

  int     max;
  int     score;
  int     scale;
  int     temp;
  int     distance;
  int     initScale;
  int     hinit;
  int     zero;

  int     offset;
  int     position;

  int     cmp;
  int     iter;
    
  __m128i *pvH;
  __m128i *pvE;

  __m128i vE, vF, vH;
  __m128i vHNext;
  __m128i vFPrev;

  __m128i vGapOpen;
  __m128i vGapExtend;
  __m128i vCeiling;

  __m128i vScale;
  __m128i vScaleAmt;
  __m128i vScaleTmp;

  __m128i vTemp;
  __m128i vNull;

  __m128i *pvScore;

  scale = 0;
  initScale = 0;

  max = 0x80000000;
  iter = (queryLength + 7) / 8;
  offset = (queryLength - 1) % iter;
  position = 7 - (queryLength - 1) / iter;

  pvH = (__m128i *)f_str->workspace;
  pvE = pvH + iter;

  /* Load gap opening penalty to all elements of a constant */
  vGapOpen = _mm_setzero_si128();	/* transfered from Apple Devel smith_waterman_sse2.c fix */
  vGapOpen = _mm_insert_epi16 (vGapOpen, gapOpen, 0);
  vGapOpen = _mm_shufflelo_epi16 (vGapOpen, 0);
  vGapOpen = _mm_shuffle_epi32 (vGapOpen, 0);

  /* Load gap extension penalty to all elements of a constant */
  vGapExtend = _mm_setzero_si128();	/* transfered from Apple Devel smith_waterman_sse2.c fix */
  vGapExtend = _mm_insert_epi16 (vGapExtend, gapExtend, 0);
  vGapExtend = _mm_shufflelo_epi16 (vGapExtend, 0);
  vGapExtend = _mm_shuffle_epi32 (vGapExtend, 0);

  /* Generate the ceiling before scaling */
  vTemp = _mm_setzero_si128();	/* transfered from Apple Devel smith_waterman_sse2.c fix */
  vTemp = _mm_insert_epi16 (vTemp, ceiling, 0);
  vTemp = _mm_shufflelo_epi16 (vTemp, 0);
  vTemp = _mm_shuffle_epi32 (vTemp, 0);
  vCeiling = _mm_cmpeq_epi16 (vTemp, vTemp);
  vCeiling = _mm_srli_epi16 (vCeiling, 1);
  vCeiling = _mm_subs_epi16 (vCeiling, vTemp);
  vCeiling = _mm_subs_epi16 (vCeiling, vGapOpen);

  vGapExtend = _mm_srli_si128 (vGapExtend, 14);
  vNull = _mm_cmpeq_epi16 (vTemp, vTemp);
  vNull = _mm_slli_epi16 (vNull, 15);
  vScaleAmt = _mm_xor_si128 (vNull, vNull);

  hinit = gapOpen * 2 - SHORT_BIAS;
  zero = hinit;

  /* Zero out the storage vector */
  vTemp = _mm_adds_epi16 (vNull, vGapOpen);
  for (i = 0; i < iter; i++) {
    _mm_store_si128 (pvH + i, vTemp);
    _mm_store_si128 (pvE + i, vNull);
  }

  /* initialize F */
  vF = vNull;
  vFPrev = vNull;

  /* load and scale H for the next round */
  vH = _mm_load_si128 (pvH + iter - 1);
  vH = _mm_slli_si128 (vH, 2);
  vH = _mm_insert_epi16 (vH, zero, 0);

  for (i = 0; i < dbLength; ++i) {
    /* fetch first data asap. */
    pvScore = (__m128i *) profile + dbSeq[i] * iter;

    vF = _mm_insert_epi16 (vNull, hinit, 0);
    vF = _mm_adds_epi16 (vF, vGapExtend);
    vF = _mm_subs_epi16 (vF, vGapOpen);

    vH = _mm_max_epi16 (vH, vFPrev);
    for (j = 0; j < iter; j++) {
      /* correct H from the previous columns F */
      vHNext = _mm_load_si128 (pvH + j);
      vHNext = _mm_max_epi16 (vHNext, vFPrev);

      /* load and correct E value */
      vE = _mm_load_si128 (pvE + j);
      vTemp = _mm_subs_epi16 (vHNext, vGapOpen);
      vE = _mm_max_epi16 (vE, vTemp);
      _mm_store_si128 (pvE + j, vE);

      /* add score to vH */
      vH = _mm_adds_epi16 (vH, *pvScore++);

      /* get max from vH, vE and vF */
      vH = _mm_max_epi16 (vH, vE);
      vH = _mm_max_epi16 (vH, vF);
      _mm_store_si128 (pvH + j, vH);

      /* update vF value */
      vH = _mm_subs_epi16 (vH, vGapOpen);
      vF = _mm_max_epi16 (vF, vH);

      /* load the next h values */
      vH = vHNext;
    }

    /* check if we need to scale before the next round */
    vTemp = _mm_cmpgt_epi16 (vF, vCeiling);
    cmp  = _mm_movemask_epi8 (vTemp);

    /* broadcast F values */
    vF = _mm_xor_si128 (vF, vNull);

    vTemp  = _mm_slli_si128 (vF, 2);
    vTemp = _mm_subs_epu16 (vTemp, vScaleAmt);
    vF = max_epu16 (vF, vTemp);

    vTemp  = _mm_slli_si128 (vF, 4);
    vScaleTmp = _mm_slli_si128 (vScaleAmt, 2);
    vScaleTmp = _mm_adds_epu16 (vScaleTmp, vScaleAmt);
    vTemp = _mm_subs_epu16 (vTemp, vScaleTmp);
    vF = max_epu16 (vF, vTemp);

    vTemp = _mm_slli_si128 (vScaleTmp, 4);
    vScaleTmp = _mm_adds_epu16 (vScaleTmp, vTemp);
    vTemp  = _mm_slli_si128 (vF, 8);
    vTemp = _mm_subs_epu16 (vTemp, vScaleTmp);
    vF = max_epu16 (vF, vTemp);

    /* scale if necessary */
    if (cmp != 0x0000)  {
      __m128i vScale1;
      __m128i vScale2;

      scale = hinit - gapOpen * 2 + SHORT_BIAS;
      initScale = initScale + scale;

      vScale = _mm_slli_si128 (vF, 2);
      vScale = _mm_subs_epu16 (vScale, vGapOpen);
      vScale = _mm_subs_epu16 (vScale, vScaleAmt);
      vScale = _mm_insert_epi16 (vScale, scale, 0);

      vTemp = _mm_slli_si128 (vScale, 2);
      vTemp = _mm_subs_epu16 (vScale, vTemp);
      vScaleAmt = _mm_adds_epu16 (vScaleAmt, vTemp);
      vTemp = _mm_slli_si128 (vScale, 2);
      vTemp = _mm_subs_epu16 (vTemp, vScale);
      vScaleAmt = _mm_subs_epu16 (vScaleAmt, vTemp);
      vTemp = _mm_subs_epu8 (vTemp, vTemp);
      vTemp = _mm_insert_epi16 (vTemp, scale, 0);
      vScaleAmt = _mm_subs_epu16 (vScaleAmt, vTemp);

      /* rescale the previous F */
      vF = _mm_subs_epu16 (vF, vScale);

      /* rescale the initial H value */
      hinit = zero;

      /* check if we can continue in signed 16-bits */
      vTemp = _mm_xor_si128 (vF, vNull);
      vTemp = _mm_cmpgt_epi16 (vTemp, vCeiling);
      cmp  = _mm_movemask_epi8 (vTemp);
      if (cmp != 0x0000) {
        return OVERFLOW_SCORE;
      }

      vTemp   = _mm_adds_epi16 (vCeiling, vCeiling);
      vScale1 = _mm_subs_epu16 (vScale, vTemp);
      vScale2 = _mm_subs_epu16 (vScale, vScale1);

      /* scale all the vectors */
      for (j = 0; j < iter; j++) {
        /* load H and E */
        vH = _mm_load_si128 (pvH + j);
        vE = _mm_load_si128 (pvE + j);

        /* get max from vH, vE and vF */
        vH = _mm_subs_epi16 (vH, vScale1);
        vH = _mm_subs_epi16 (vH, vScale2);
        vE = _mm_subs_epi16 (vE, vScale1);
        vE = _mm_subs_epi16 (vE, vScale2);

        /* save the H and E */
        _mm_store_si128 (pvH + j, vH);
        _mm_store_si128 (pvE + j, vE);
      }

      vScale = vScaleAmt;
      for (j = 0; j < position; ++j) {
        vScale = _mm_slli_si128 (vScale, 2);
      }

      /* calculate the final scaling amount */
      /* vTemp   = _mm_xor_si128 (vTemp, vTemp); */
      vTemp = _mm_setzero_si128();	/* transfered from Apple Devel fix for smith_waterman_sse2.c */
      vScale1 = _mm_unpacklo_epi16 (vScale, vTemp);
      vScale2 = _mm_unpackhi_epi16 (vScale, vTemp);
      vScale  = _mm_add_epi32 (vScale1, vScale2);
      vTemp  = _mm_srli_si128 (vScale, 8);
      vScale = _mm_add_epi32 (vScale, vTemp);
      vTemp  = _mm_srli_si128 (vScale, 4);
      vScale = _mm_add_epi32 (vScale, vTemp);
      scale  = (int) (unsigned short) _mm_extract_epi16 (vScale, 0);
      temp   = (int) (unsigned short) _mm_extract_epi16 (vScale, 1);
      scale  = scale + (temp << 16) + initScale;
    }

    /* scale the F value for the next round */
    vFPrev = _mm_slli_si128 (vF, 2);
    vFPrev = _mm_subs_epu16 (vFPrev, vScaleAmt);
    vFPrev = _mm_xor_si128 (vFPrev, vNull);

    vF = _mm_xor_si128 (vF, vNull);

    vH = _mm_load_si128 (pvH + offset);
    vH = _mm_max_epi16 (vH, vFPrev);
    for (j = 0; j < position; ++j) {
      vH = _mm_slli_si128 (vH, 2);
    }
    score = (int) (signed short) _mm_extract_epi16 (vH, 7);
    score = score + SHORT_BIAS;

    /* return largest score */
    distance = (queryLength + i + 1) * gapExtend;
    score = score - (gapOpen * 2) - distance + scale;
    max = (max > score) ? max : score;

    /* load and scale H for the next round */
    hinit += gapExtend;
    vH = _mm_load_si128 (pvH + iter - 1);
    vH = _mm_slli_si128 (vH, 2);
    vH = _mm_xor_si128 (vH, vNull);
    vH = _mm_subs_epu16 (vH, vScaleAmt);
    vH = _mm_xor_si128 (vH, vNull);
    vH = _mm_insert_epi16 (vH, hinit, 0);
  }

  return max;
}

int
glocal_sse2_byte(int                  queryLength,
                 unsigned char       *profile,
                 const unsigned char *dbSeq,
                 int                  dbLength,
                 unsigned short       gapOpen,
                 unsigned short       gapExtend,
                 unsigned short       ceiling,
                 unsigned short       bias,
                 struct f_struct     *f_str)
{
  int     i, j;

  int     max;
  int     score;
  int     scale;
  int     distance;
  int     initScale;

  int     offset;
  int     position;

  int     dup;
  int     cmp;
  int     iter;
    
  __m128i *pvH;
  __m128i *pvE;

  __m128i vE, vF, vH;
  __m128i vHInit;
  __m128i vHNext;
  __m128i vFPrev;

  __m128i vBias;
  __m128i vGapOpen;
  __m128i vGapExtend;
  __m128i vCeiling;

  __m128i vScale;
  __m128i vScaleAmt;
  __m128i vScaleTmp;

  __m128i vTemp;
  __m128i vZero;
  __m128i vNull;

  __m128i *pvScore;

  scale = 0;
  initScale = 0;

  max = 0x80000000;
  iter = (queryLength + 15) / 16;
  offset = (queryLength - 1) % iter;
  position = 15 - (queryLength - 1) / iter;

  pvH = (__m128i *)f_str->workspace;
  pvE = pvH + iter;

  /* Load the bias to all elements of a constant */
  dup    = (bias << 8) | (bias & 0x00ff);
  vBias = _mm_setzero_si128();	/* transfered from Apple Devel smith_waterman_sse2.c fix */
  vBias = _mm_insert_epi16 (vBias, dup, 0);
  vBias = _mm_shufflelo_epi16 (vBias, 0);
  vBias = _mm_shuffle_epi32 (vBias, 0);

  /* Load gap opening penalty to all elements of a constant */
  dup      = (gapOpen << 8) | (gapOpen & 0x00ff);
  vGapOpen = _mm_setzero_si128();	/* transfered from Apple Devel smith_waterman_sse2.c fix */
  vGapOpen = _mm_insert_epi16 (vGapOpen, dup, 0);
  vGapOpen = _mm_shufflelo_epi16 (vGapOpen, 0);
  vGapOpen = _mm_shuffle_epi32 (vGapOpen, 0);

  /* Load gap extension penalty to all elements of a constant */
  dup    = (gapExtend << 8) | (gapExtend & 0x00ff);
  vGapExtend = _mm_setzero_si128();	/* transfered from Apple Devel smith_waterman_sse2.c fix */
  vGapExtend = _mm_insert_epi16 (vGapExtend, dup, 0);
  vGapExtend = _mm_shufflelo_epi16 (vGapExtend, 0);
  vGapExtend = _mm_shuffle_epi32 (vGapExtend, 0);

  /* Generate the ceiling before scaling */
  dup    = (ceiling << 8) | (ceiling & 0x00ff);
  vTemp = _mm_setzero_si128();	/* transfered from Apple Devel smith_waterman_sse2.c fix */
  vTemp = _mm_insert_epi16 (vTemp, dup, 0);
  vTemp = _mm_shufflelo_epi16 (vTemp, 0);
  vTemp = _mm_shuffle_epi32 (vTemp, 0);
  vCeiling = _mm_cmpeq_epi8 (vTemp, vTemp);
  vCeiling = _mm_subs_epu8 (vCeiling, vTemp);
  vCeiling = _mm_subs_epu8 (vCeiling, vGapOpen);

  /* since we want to use the full range, zero is redefined as */
  /* 2 * gapOpen.  the lowest scaled score will an insert followed */
  /* by a delete. */
  vHInit = _mm_adds_epu8 (vGapOpen, vGapOpen);
  vHInit = _mm_srli_si128 (vHInit, 15);
  vZero  = vHInit;

  vGapExtend = _mm_srli_si128 (vGapExtend, 15);
  /*   vNull = _mm_xor_si128 (vNull, vNull); */
  vNull = _mm_setzero_si128();	/* transfered from Apple Devel smith_waterman_sse2.c fix */
  vScaleAmt = vNull;

  /* Zero out the storage vector */
  for (i = 0; i < iter; i++) {
    _mm_store_si128 (pvH + i, vGapOpen);
    _mm_store_si128 (pvE + i, vNull);
  }

  /* initialize F */
  vF = vNull;
  vFPrev = vNull;

  /* load and scale H for the next round */
  vH = _mm_load_si128 (pvH + iter - 1);
  vH = _mm_slli_si128 (vH, 1);
  vH = _mm_or_si128 (vH, vZero);

  for (i = 0; i < dbLength; ++i) {
    /* fetch first data asap. */
    pvScore = (__m128i *) profile + dbSeq[i] * iter;

    vF = _mm_adds_epu8 (vHInit, vGapExtend);
    vF = _mm_subs_epu8 (vF, vGapOpen);

    vH = _mm_max_epu8 (vH, vFPrev);
    for (j = 0; j < iter; j++) {
      /* correct H from the previous columns F */
      vHNext = _mm_load_si128 (pvH + j);
      vHNext = _mm_max_epu8 (vHNext, vFPrev);

      /* load and correct E value */
      vE = _mm_load_si128 (pvE + j);
      vTemp = _mm_subs_epu8 (vHNext, vGapOpen);
      vE = _mm_max_epu8 (vE, vTemp);
      _mm_store_si128 (pvE + j, vE);

      /* add score to vH */
      vH = _mm_adds_epu8 (vH, *pvScore++);
      vH = _mm_subs_epu8 (vH, vBias);

      /* get max from vH, vE and vF */
      vH = _mm_max_epu8 (vH, vE);
      vH = _mm_max_epu8 (vH, vF);
      _mm_store_si128 (pvH + j, vH);

      /* update vF value */
      vH = _mm_subs_epu8 (vH, vGapOpen);
      vF = _mm_max_epu8 (vF, vH);

      /* load the next h values */
      vH = vHNext;
    }

    /* check if we need to scale before the next round */
    vTemp = _mm_subs_epu8 (vCeiling, vF);
    vTemp = _mm_cmpeq_epi8 (vTemp, vNull);
    cmp  = _mm_movemask_epi8 (vTemp);

    /* broadcast F values */
    vTemp  = _mm_slli_si128 (vF, 1);
    vTemp = _mm_subs_epu8 (vTemp, vScaleAmt);
    vF = _mm_max_epu8 (vF, vTemp);

    vScaleTmp = _mm_slli_si128 (vScaleAmt, 1);
    vScaleTmp = _mm_adds_epu8 (vScaleTmp, vScaleAmt);
    vTemp  = _mm_slli_si128 (vF, 2);
    vTemp = _mm_subs_epu8 (vTemp, vScaleTmp);
    vF = _mm_max_epu8 (vF, vTemp);

    vTemp = _mm_slli_si128 (vScaleTmp, 2);
    vScaleTmp = _mm_adds_epu8 (vScaleTmp, vTemp);
    vTemp  = _mm_slli_si128 (vF, 4);
    vTemp = _mm_subs_epu8 (vTemp, vScaleTmp);
    vF = _mm_max_epu8 (vF, vTemp);

    vTemp = _mm_slli_si128 (vScaleTmp, 4);
    vScaleTmp = _mm_adds_epu8 (vScaleTmp, vTemp);
    vTemp  = _mm_slli_si128 (vF, 8);
    vTemp = _mm_subs_epu8 (vTemp, vScaleTmp);
    vF = _mm_max_epu8 (vF, vTemp);

    /* scale if necessary */
    if (cmp != 0x0000) {
      vHInit = _mm_subs_epu8 (vHInit, vGapOpen);
      vHInit = _mm_subs_epu8 (vHInit, vGapOpen);
      scale = _mm_extract_epi16 (vHInit, 0);
      initScale = initScale + scale;

      vScale = _mm_slli_si128 (vF, 1);
      vScale = _mm_subs_epu8 (vScale, vGapOpen);
      vScale = _mm_subs_epu8 (vScale, vScaleAmt);
      vScale = _mm_or_si128 (vScale, vHInit);

      vTemp = _mm_slli_si128 (vScale, 1);
      vTemp = _mm_subs_epu8 (vScale, vTemp);
      vScaleAmt = _mm_adds_epu8 (vScaleAmt, vTemp);
      vTemp = _mm_slli_si128 (vScale, 1);
      vTemp = _mm_subs_epu8 (vTemp, vScale);
      vScaleAmt = _mm_subs_epu8 (vScaleAmt, vTemp);
      vScaleAmt = _mm_subs_epu8 (vScaleAmt, vHInit);

      /* rescale the previous F */
      vF = _mm_subs_epu8 (vF, vScale);

      /* rescale the initial H value */
      vHInit = vZero;

      /* check if we can continue in 8-bits */
      vTemp = _mm_subs_epu8 (vCeiling, vF);
      vTemp = _mm_cmpeq_epi8 (vTemp, vNull);
      cmp  = _mm_movemask_epi8 (vTemp);
      if (cmp != 0x0000) {
        return OVERFLOW_SCORE;
      }

      /* scale all the vectors */
      for (j = 0; j < iter; j++) {
        /* load H and E */
        vH = _mm_load_si128 (pvH + j);
        vE = _mm_load_si128 (pvE + j);

        /* get max from vH, vE and vF */
        vH = _mm_subs_epu8 (vH, vScale);
        vE = _mm_subs_epu8 (vE, vScale);

        /* save the H and E */
        _mm_store_si128 (pvH + j, vH);
        _mm_store_si128 (pvE + j, vE);
      }

      /* calculate the final scaling amount */
      vScale = vScaleAmt;
      for (j = 0; j < position; ++j) {
        vScale = _mm_slli_si128 (vScale, 1);
      }
      vTemp = _mm_unpacklo_epi8 (vScale, vNull);
      vScale = _mm_unpackhi_epi8 (vScale, vNull);
      vScale = _mm_adds_epi16 (vScale, vTemp);
      vTemp = _mm_srli_si128 (vScale, 8);
      vScale = _mm_adds_epi16 (vScale, vTemp);
      vTemp = _mm_srli_si128 (vScale, 4);
      vScale = _mm_adds_epi16 (vScale, vTemp);
      vTemp = _mm_srli_si128 (vScale, 2);
      vScale = _mm_adds_epi16 (vScale, vTemp);
      scale = (int) _mm_extract_epi16 (vScale, 0);
      scale = scale + initScale;
    }

    /* scale the F value for the next round */
    vFPrev = _mm_slli_si128 (vF, 1);
    vFPrev = _mm_subs_epu8 (vFPrev, vScaleAmt);

    /* calculate the max glocal score for this column */
    vH = _mm_load_si128 (pvH + offset);
    vH = _mm_max_epu8 (vH, vF);
    for (j = 0; j < position; ++j) {
      vH = _mm_slli_si128 (vH, 1);
    }
    score = (int) (unsigned short) _mm_extract_epi16 (vH, 7);
    score >>= 8;

    /* return largest score */
    distance = (queryLength + i + 1) * gapExtend;
    score = score - (gapOpen * 2) - distance + scale;
    max = (max > score) ? max : score;

    /* load and scale H for the next round */
    vHInit = _mm_adds_epu8 (vHInit, vGapExtend);
    vH = _mm_load_si128 (pvH + iter - 1);
    vH = _mm_slli_si128 (vH, 1);
    vH = _mm_subs_epu8 (vH, vScaleAmt);
    vH = _mm_or_si128 (vH, vHInit);
  }

  return max;
}
#else

/* No SSE2 support. Avoid compiler complaints about empty object */

int nw_dummy;

#endif

