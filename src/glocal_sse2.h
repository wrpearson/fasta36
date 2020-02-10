/******************************************************************
  Copyright 2010 Michael Farrar.  
  Freely distributed under the BSD open source license.
  See the ../COPYRIGHT.sse2 file for details.
*******************************************************************/

/*
  Written by Michael Farrar, 2010.
*/

#ifndef INCLUDE_GLOCAL_SSE2_H
#define INCLUDE_GLOCAL_SSE2_H

#define SHORT_BIAS     32768
#define OVERFLOW_SCORE 0x7f000000

int
glocal_sse2_word(int                  queryLength,
                 unsigned short      *profile,
                 const unsigned char *dbSeq,
                 int                  dbLength,
                 unsigned short       gapOpen,
                 unsigned short       gapExtend,
                 unsigned short       ceiling,
                 struct f_struct     *f_str);

int
glocal_sse2_byte(int                  queryLength,
                 unsigned char       *profile,
                 const unsigned char *dbSeq,
                 int                  dbLength,
                 unsigned short       gapOpen,
                 unsigned short       gapExtend,
                 unsigned short       ceiling,
                 unsigned short       bias,
                 struct f_struct     *f_str);

#endif /* INCLUDE_GLOCAL_SSE2_H */
