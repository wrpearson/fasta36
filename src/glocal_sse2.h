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
