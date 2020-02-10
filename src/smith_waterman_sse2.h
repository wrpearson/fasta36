
/* $Id: smith_waterman_sse2.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

/******************************************************************
  Copyright 2006 Michael Farrar.  
  Freely distributed under the BSD open source license.
  See the ../COPYRIGHT.sse2 file for details.
*******************************************************************/

/*
  Written by Michael Farrar, 2006.
  For information, contact Sean Eddy (sean@eddylab.org).
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
