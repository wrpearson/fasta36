/* $Id: smith_waterman_altivec.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

int
smith_waterman_altivec_word(const unsigned char *     query_sequence,
                            unsigned short *    query_profile_word,
                            const int                 query_length,
                            const unsigned char *     db_sequence,
                            const int                 db_length,
                            unsigned short      bias,
                            unsigned short      gap_open,
                            unsigned short      gap_extend,
                            struct f_struct *   f_str);


int
smith_waterman_altivec_byte(const unsigned char *     query_sequence,
                            unsigned char *     query_profile_byte,
                            const int                 query_length,
                            const unsigned char *     db_sequence,
                            const int                 db_length,
                            unsigned char       bias,
                            unsigned char       gap_open,
                            unsigned char       gap_extend,
                            struct f_struct *   f_str);

