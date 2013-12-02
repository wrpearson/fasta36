
/* $Id: smith_waterman_sse2.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

/******************************************************************
  Copyright 2006 by Michael Farrar.  All rights reserved.
  This program may not be sold or incorporated into a commercial product,
  in whole or in part, without written consent of Michael Farrar.  For 
  further information regarding permission for use or reproduction, please 
  contact: Michael Farrar at farrar.michael@gmail.com.
*******************************************************************/

/*
  Written by Michael Farrar, 2006.
  Please send bug reports and/or suggestions to farrar.michael@gmail.com.
*/

#ifndef SMITH_WATERMAN_SSE2_H
#define SMITH_WATERMAN_SSE2_H

int
smith_waterman_sse2_word(const unsigned char *     query_sequence,
                         unsigned short *    query_profile_word,
                         const int                 query_length,
                         const unsigned char *     db_sequence,
                         const int                 db_length,
                         unsigned short      gap_open,
                         unsigned short      gap_extend,
                         struct f_struct *   f_str);


int
smith_waterman_sse2_byte(const unsigned char *     query_sequence,
                         unsigned char *     query_profile_byte,
                         const int                 query_length,
                         const unsigned char *     db_sequence,
                         const int                 db_length,
                         unsigned char       bias,
                         unsigned char       gap_open,
                         unsigned char       gap_extend,
                         struct f_struct *   f_str);

#endif /* SMITH_WATERMAN_SSE2_H */
