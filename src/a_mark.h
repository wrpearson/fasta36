/* a_mark.h - symbols used to indicate match/mismatch alignment code */

/* $Id: a_mark.h 1024 2012-08-07 18:08:45Z wrp $ */

/* copyright (c) 2003 by William R. Pearson and The Rector & Vistors
   of the University of Virginia */

/* Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing,
   software distributed under this License is distributed on an "AS
   IS" BASIS, WITHOUT WRRANTIES OR CONDITIONS OF ANY KIND, either
   express or implied.  See the License for the specific language
   governing permissions and limitations under the License. 
*/

/* character types in aln_code strings */
#define M_BLANK 0
#define M_NEG 1
#define M_ZERO 2
#define M_POS 3
#define M_IDENT 4
#define M_DEL 5

#define MX_A0 0
#define MX_A1 1
#define MX_A2 2
#ifdef M10_CONS_L
#define MX_A10 3
#else
#define MX_A10 4
#endif
#define MX_ACC 5
#define MX_ABLAST 6

static char *
aln_map_sym[] = {"  ..: ",	/* 0 */
		 " Xxx  ",	/* 1 */
		 "    . ",	/* 2 */
		 " mzp=-",	/* 3: 10a */
		 "  ..:-",	/* 4: 10b */
		 " <z>=-",	/* 5: calc_code */
		 "   += "	/* 6: MX_MBLAST blast */
		  };



