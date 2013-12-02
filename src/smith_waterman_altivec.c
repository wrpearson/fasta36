
/* Implementation of the Wozniak "anti-diagonal" vectorization
   strategy for Smith-Waterman comparison, Wozniak (1997) Comp.
   Appl. Biosci. 13:145-150

   November, 2004
*/

/*
  Written by Erik Lindahl, Stockholm Bioinformatics Center, 2004.
  Please send bug reports and/or suggestions to lindahl@sbc.su.se.
*/

#include <stdio.h>

#include "defs.h"
#include "param.h"
#include "dropgsw2.h"

#ifdef SW_ALTIVEC

int
smith_waterman_altivec_word(unsigned char *     query_sequence,
                            unsigned short *    query_profile_word,
                            int                 query_length,
                            unsigned char *     db_sequence,
                            int                 db_length,
                            unsigned short      bias,
                            unsigned short      gap_open,
                            unsigned short      gap_extend,
                            struct f_struct *   f_str)
{
    int                     i,j,k;
    unsigned short *        p;
    unsigned short          score;   
    unsigned char *         p_dbseq;
    int                     alphabet_size = f_str->alphabet_size;
    unsigned short *        workspace     = (unsigned short *)f_str->workspace;

    vector unsigned short   Fup,Hup1,Hup2,E,F,H,tmp;
    vector unsigned char    perm;
    vector unsigned short   v_maxscore;
    vector unsigned short   v_bias,v_gapopen,v_gapextend;
    vector unsigned short   v_score;
    vector unsigned short   v_score_q1;
    vector unsigned short   v_score_q2;
    vector unsigned short   v_score_q3;
    vector unsigned short   v_score_load; 
    vector unsigned char    queue1_to_score  = (vector unsigned char)(16,17,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    vector unsigned char    queue2_to_queue1 = (vector unsigned char)(0,1,18,19,4,5,6,7,8,9,10,11,12,13,14,15);
    vector unsigned char    queue3_to_queue2 = (vector unsigned char)(16,16,16,16,16,21,16,0,16,1,16,2,16,3,16,4);
    vector unsigned char    queue3_with_load = (vector unsigned char)(23,5,6,7,8,25,9,10,11,27,12,13,29,14,31,16);
        
    /* Load the bias to all elements of a constant */
    v_bias           = vec_lde(0,&bias);
    perm             = vec_lvsl(0,&bias);
    v_bias           = vec_perm(v_bias,v_bias,perm);
    v_bias           = vec_splat(v_bias,0);
    
    /* Load gap opening penalty to all elements of a constant */
    v_gapopen        = vec_lde(0,&gap_open);
    perm             = vec_lvsl(0,&gap_open);
    v_gapopen        = vec_perm(v_gapopen,v_gapopen,perm);
    v_gapopen        = vec_splat(v_gapopen,0);

    /* Load gap extension penalty to all elements of a constant */
    v_gapextend      = vec_lde(0,&gap_extend);  
    perm             = vec_lvsl(0,&gap_extend);
    v_gapextend      = vec_perm(v_gapextend,v_gapextend,perm);
    v_gapextend      = vec_splat(v_gapextend,0);
    
    v_maxscore = vec_xor(v_maxscore,v_maxscore);
   
    // Zero out the storage vector 
    k = 2*(db_length+7);
        
    for(i=0,j=0;i<k;i++,j+=16)
    {
        // borrow the zero value in v_maxscore to have something to store
        vec_st(v_maxscore,j,workspace);
    }
    
    for(i=0;i<query_length;i+=8)
    {
        // fetch first data asap.
        p_dbseq    = db_sequence;
        k          = *p_dbseq++;
        v_score_load = vec_ld(16*k,query_profile_word);

        // zero lots of stuff. 
        // We use both the VPERM and VSIU unit to knock off some cycles.
        
        E          = vec_splat_u16(0);
        F          = vec_xor(F,F);
        H          = vec_splat_u16(0);
        Hup2       = vec_xor(Hup2,Hup2);
        v_score_q1 = vec_splat_u16(0);
        v_score_q2 = vec_xor(v_score_q2,v_score_q2);
        v_score_q3 = vec_splat_u16(0);

        // reset pointers to the start of the saved data from the last row
        p = workspace;
                
        // PROLOGUE 1
        // prefetch next residue
        k          = *p_dbseq++;
        
        // Create the actual diagonal score vector
        // and update the queue of incomplete score vectors
        
        v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
        v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
        v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
        v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
        
        // prefetch score for next step 
        v_score_load = vec_ld(16*k,query_profile_word);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 16; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,14);
        Hup1   = vec_sld(Hup1,H,14);            
        
        // do the dynamic programming 

        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Save value to use for next diagonal H 
        Hup2 = Hup1;
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        // PROLOGUE 2
        // prefetch next residue
        k          = *p_dbseq++;
        
        // Create the actual diagonal score vector
        // and update the queue of incomplete score vectors
        
        v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
        v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
        v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
        v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
        
        // prefetch score for next step 
        v_score_load = vec_ld(16*k,query_profile_word);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 16; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,14);
        Hup1   = vec_sld(Hup1,H,14);            
        
        // do the dynamic programming 

        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Save value to use for next diagonal H 
        Hup2 = Hup1;
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        

        // PROLOGUE 3
        // prefetch next residue
        k          = *p_dbseq++;
        
        // Create the actual diagonal score vector
        // and update the queue of incomplete score vectors
        
        v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
        v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
        v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
        v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);

        // prefetch score for next step 
        v_score_load = vec_ld(16*k,query_profile_word);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 16; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,14);
        Hup1   = vec_sld(Hup1,H,14);            
        
        // do the dynamic programming 

        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Save value to use for next diagonal H 
        Hup2 = Hup1;
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        

        // PROLOGUE 4
        // prefetch next residue
        k          = *p_dbseq++;
        
        // Create the actual diagonal score vector
        // and update the queue of incomplete score vectors
        
        v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
        v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
        v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
        v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
        
        // prefetch score for next step 
        v_score_load = vec_ld(16*k,query_profile_word);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 16; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,14);
        Hup1   = vec_sld(Hup1,H,14);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Save value to use for next diagonal H 
        Hup2 = Hup1;
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        

        // PROLOGUE 5
        // prefetch next residue
        k          = *p_dbseq++;
        
        // Create the actual diagonal score vector
        // and update the queue of incomplete score vectors
        
        v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
        v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
        v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
        v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
        
        // prefetch score for next step 
        v_score_load = vec_ld(16*k,query_profile_word);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 16; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,14);
        Hup1   = vec_sld(Hup1,H,14);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Save value to use for next diagonal H 
        Hup2 = Hup1;
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        

        // PROLOGUE 6
        // prefetch next residue
        k          = *p_dbseq++;
        
        // Create the actual diagonal score vector
        // and update the queue of incomplete score vectors
        
        v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
        v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
        v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
        v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
        
        // prefetch score for next step 
        v_score_load = vec_ld(16*k,query_profile_word);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 16; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,14);
        Hup1   = vec_sld(Hup1,H,14);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Save value to use for next diagonal H 
        Hup2 = Hup1;
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);

        
        // PROLOGUE 7
        // prefetch next residue
        k          = *p_dbseq++;
        
        // Create the actual diagonal score vector
        // and update the queue of incomplete score vectors
        
        v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
        v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
        v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
        v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
        
        // prefetch score for next step 
        v_score_load = vec_ld(16*k,query_profile_word);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 16; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,14);
        Hup1   = vec_sld(Hup1,H,14);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Save value to use for next diagonal H 
        Hup2 = Hup1;
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        

        // PROLOGUE 8
        // prefetch next residue
        k          = *p_dbseq++;
        
        // Create the actual diagonal score vector
        // and update the queue of incomplete score vectors
        
        v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
        v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
        v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
        v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
        
        // prefetch score for next step 
        v_score_load = vec_ld(16*k,query_profile_word);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 16; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,14);
        Hup1   = vec_sld(Hup1,H,14);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Save value to use for next diagonal H 
        Hup2 = Hup1;
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
    

        // reset pointers to the start of the saved data from the last row
        p = workspace;

        for(j=8;j<db_length;j+=8)
        {           
            // STEP 1
            
            // prefetch next residue
            k          = *p_dbseq++;
            
            // Create the actual diagonal score vector
            // and update the queue of incomplete score vectors

            v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
            v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
            v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
            v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
            
            // prefetch score for next step
            v_score_load = vec_ld(16*k,query_profile_word);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(256, p);
            Hup1   = vec_ld(272, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 16; // move ahead 32 bytes
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,14);
            Hup1   = vec_sld(Hup1,H,14);            

            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);

            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F); 
            
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
 
 
            // STEP 2
            
            // prefetch next residue
            k          = *p_dbseq++;
            
            // Create the actual diagonal score vector
            // and update the queue of incomplete score vectors
            
            v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
            v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
            v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
            v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
            
            // prefetch score for next step
            v_score_load = vec_ld(16*k,query_profile_word);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(256, p);
            Hup2   = vec_ld(272, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 16; // move ahead 32 bytes
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,14);
            Hup2   = vec_sld(Hup2,H,14);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F); 
            
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            


            // STEP 3
            
            // prefetch next residue
            k          = *p_dbseq++;
            
            // Create the actual diagonal score vector
            // and update the queue of incomplete score vectors
            
            v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
            v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
            v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
            v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
            
            // prefetch score for next step
            v_score_load = vec_ld(16*k,query_profile_word);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(256, p);
            Hup1   = vec_ld(272, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 16; // move ahead 32 bytes
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,14);
            Hup1   = vec_sld(Hup1,H,14);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F); 
            

            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            

            
            // STEP 4
            
            // prefetch next residue
            k          = *p_dbseq++;
            
            // Create the actual diagonal score vector
            // and update the queue of incomplete score vectors
            
            v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
            v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
            v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
            v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
            
            // prefetch score for next step
            v_score_load = vec_ld(16*k,query_profile_word);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(256, p);
            Hup2   = vec_ld(272, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 16; // move ahead 32 bytes
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,14);
            Hup2   = vec_sld(Hup2,H,14);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F); 

            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            


            // STEP 5
            
            // prefetch next residue
            k          = *p_dbseq++;
            
            // Create the actual diagonal score vector
            // and update the queue of incomplete score vectors
            
            v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
            v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
            v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
            v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
            
            // prefetch score for next step
            v_score_load = vec_ld(16*k,query_profile_word);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(256, p);
            Hup1   = vec_ld(272, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 16; // move ahead 32 bytes
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,14);
            Hup1   = vec_sld(Hup1,H,14);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F); 
            
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            


            // STEP 6
            
            // prefetch next residue
            k          = *p_dbseq++;
            
            // Create the actual diagonal score vector
            // and update the queue of incomplete score vectors
            
            v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
            v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
            v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
            v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
            
            // prefetch score for next step
            v_score_load = vec_ld(16*k,query_profile_word);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(256, p);
            Hup2   = vec_ld(272, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 16; // move ahead 32 bytes
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,14);
            Hup2   = vec_sld(Hup2,H,14);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F); 
            

            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            

            
            // STEP 7
            
            // prefetch next residue
            k          = *p_dbseq++;
            
            // Create the actual diagonal score vector
            // and update the queue of incomplete score vectors
            
            v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
            v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
            v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
            v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
            
            // prefetch score for next step
            v_score_load = vec_ld(16*k,query_profile_word);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(256, p);
            Hup1   = vec_ld(272, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 16; // move ahead 32 bytes
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,14);
            Hup1   = vec_sld(Hup1,H,14);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F); 
            

            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            


            // STEP 8
            
            // prefetch next residue
            k          = *p_dbseq++;
            
            // Create the actual diagonal score vector
            // and update the queue of incomplete score vectors
            
            v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
            v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
            v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
            v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
            
            // prefetch score for next step
            v_score_load = vec_ld(16*k,query_profile_word);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(256, p);
            Hup2   = vec_ld(272, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 16; // move ahead 32 bytes
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,14);
            Hup2   = vec_sld(Hup2,H,14);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F); 
            
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
        }
        
        v_score_load = vec_splat_u16(0);
        
        for(;j<db_length+7;j++)
        {
            // Create the actual diagonal score vector
            // and update the queue of incomplete score vectors
            //
            // This could of course be done with only vec_perm or vec_sel,
            // but since they use different execution units we have found
            // it to be slightly faster to mix them.
            v_score    = vec_perm(v_score_q1, v_score_load, queue1_to_score);
            v_score_q1 = vec_perm(v_score_q2, v_score_load, queue2_to_queue1);
            v_score_q2 = vec_perm(v_score_q3, v_score_load, queue3_to_queue2);
            v_score_q3 = vec_perm(v_score_q3, v_score_load, queue3_with_load);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 16; // move ahead 32 bytes
            
            // v_score_load contains all zeros
            Fup    = vec_sld(v_score_load,F,14);
            Hup1   = vec_sld(v_score_load,H,14);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Save value to use for next diagonal H 
            Hup2 = Hup1;
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
        }
        vec_st(F, 0,  p);
        vec_st(H, 16, p);

        query_profile_word += 8*alphabet_size;
    }

    // find largest score in the v_maxscore vector
    tmp = vec_sld(v_maxscore,v_maxscore,8);
    v_maxscore = vec_max(v_maxscore,tmp);
    tmp = vec_sld(v_maxscore,v_maxscore,4);
    v_maxscore = vec_max(v_maxscore,tmp);
    tmp = vec_sld(v_maxscore,v_maxscore,2);
    v_maxscore = vec_max(v_maxscore,tmp);

    // store in temporary variable
    vec_ste(v_maxscore,0,&score);
    
    // return largest score
    return score;
}

int
smith_waterman_altivec_byte(unsigned char *     query_sequence,
                            unsigned char *     query_profile_byte,
                            int                 query_length,
                            unsigned char *     db_sequence,
                            int                 db_length,
                            unsigned char       bias,
                            unsigned char       gap_open,
                            unsigned char       gap_extend,
                            struct f_struct *   f_str)
{
    int                     i,j,k,k8;
    int                     overflow;
    unsigned char *         p;
    unsigned char           score;   
    int                     alphabet_size = f_str->alphabet_size;
    unsigned char *         workspace     = (unsigned char *)f_str->workspace;
    
    vector unsigned char    Fup,Hup1,Hup2,E,F,H,tmp;
    vector unsigned char    perm;
    vector unsigned char    v_maxscore;
    vector unsigned char    v_bias,v_gapopen,v_gapextend;
    vector unsigned char    v_score;
    vector unsigned char    v_score_q1;
    vector unsigned char    v_score_q2;
    vector unsigned char    v_score_q3;
    vector unsigned char    v_score_q4;
    vector unsigned char    v_score_q5;
    vector unsigned char    v_score_load1;
    vector unsigned char    v_score_load2;  
    vector unsigned char    v_zero;  

    vector unsigned char    queue1_to_score  = (vector unsigned char)(16,1,2,3,4,5,6,7,24,9,10,11,12,13,14,15);
    vector unsigned char    queue2_to_queue1 = (vector unsigned char)(16,17,2,3,4,5,6,7,24,25,10,11,12,13,14,15);
    vector unsigned char    queue3_to_queue2 = (vector unsigned char)(16,17,18,3,4,5,6,7,24,25,26,11,12,13,14,15);
    vector unsigned char    queue4_to_queue3 = (vector unsigned char)(16,17,18,19,4,5,6,7,24,25,26,27,12,13,14,15);
    vector unsigned char    queue5_to_queue4 = (vector unsigned char)(16,17,18,19,20,2,3,4,24,25,26,27,28,10,11,12);
    vector unsigned char    queue5_with_load = (vector unsigned char)(19,20,21,5,6,22,7,23,27,28,29,13,14,30,15,31);
    vector unsigned char    merge_score_load = (vector unsigned char)(0,1,2,3,4,5,6,7,24,25,26,27,28,29,30,31);

    v_zero           = vec_splat_u8(0);
        
    /* Load the bias to all elements of a constant */
    v_bias           = vec_lde(0,&bias);
    perm             = vec_lvsl(0,&bias);
    v_bias           = vec_perm(v_bias,v_bias,perm);
    v_bias           = vec_splat(v_bias,0);
    
    /* Load gap opening penalty to all elements of a constant */
    v_gapopen        = vec_lde(0,&gap_open);
    perm             = vec_lvsl(0,&gap_open);
    v_gapopen        = vec_perm(v_gapopen,v_gapopen,perm);
    v_gapopen        = vec_splat(v_gapopen,0);

    /* Load gap extension penalty to all elements of a constant */
    v_gapextend      = vec_lde(0,&gap_extend);  
    perm             = vec_lvsl(0,&gap_extend);
    v_gapextend      = vec_perm(v_gapextend,v_gapextend,perm);
    v_gapextend      = vec_splat(v_gapextend,0);
    
    v_maxscore = vec_xor(v_maxscore,v_maxscore);
   
    // Zero out the storage vector 
    k = (db_length+15);
    for(i=0,j=0;i<k;i++,j+=32)
    {
        // borrow the zero value in v_maxscore to have something to store
        vec_st(v_maxscore,j,workspace);
        vec_st(v_maxscore,j+16,workspace);
    }
    
    for(i=0;i<query_length;i+=16)
    {
        // zero lots of stuff. 
        // We use both the VPERM and VSIU unit to knock off some cycles.
        
        E          = vec_splat_u8(0);
        F          = vec_xor(F,F);
        H          = vec_splat_u8(0);
        Hup2      = vec_xor(Hup2,Hup2);
        v_score_q1 = vec_splat_u8(0);
        v_score_q2 = vec_xor(v_score_q2,v_score_q2);
        v_score_q3 = vec_splat_u8(0);
        v_score_q4 = vec_xor(v_score_q4,v_score_q4);
        v_score_q5 = vec_splat_u8(0);

        // reset pointers to the start of the saved data from the last row
        p = workspace;
        
        // start directly and prefetch score column
        k             = db_sequence[0];
        k8            = k;
        v_score_load1 = vec_ld(16*k,query_profile_byte);
        v_score_load2 = v_score_load1;
        v_score_load1 = vec_perm(v_score_load1,v_zero,merge_score_load);

        // PROLOGUE 1
        // prefetch next residue
        k                = db_sequence[1];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        
        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup1    = vec_sld(Hup1,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_zero,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        
        // PROLOGUE 2
        // prefetch next residue
        k                = db_sequence[2];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        
  
        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup2   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup2   = vec_sld(Hup2,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup2,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_zero,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup1,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
     
        
        // PROLOGUE 3
        // prefetch next residue
        k                = db_sequence[3];
  
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup1    = vec_sld(Hup1,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_zero,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        // PROLOGUE 4
        // prefetch next residue
        k                = db_sequence[4];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup2   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup2   = vec_sld(Hup2,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup2,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_zero,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup1,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        // PROLOGUE 5
        // prefetch next residue
        k                = db_sequence[5];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
     

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup1    = vec_sld(Hup1,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_zero,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        // PROLOGUE 6
        // prefetch next residue
        k                = db_sequence[6];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup2   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup2   = vec_sld(Hup2,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup2,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_zero,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup1,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        // PROLOGUE 7
        // prefetch next residue
        k                = db_sequence[7];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup1    = vec_sld(Hup1,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_zero,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        // PROLOGUE 8
        // prefetch next residue
        k                = db_sequence[8];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup2   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup2   = vec_sld(Hup2,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup2,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup1,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        
        // PROLOGUE 9
        // prefetch next residue
        k                = db_sequence[9];
        k8               = db_sequence[1];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        v_score_load2 = vec_ld(16*k8,query_profile_byte);
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1    = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup1    = vec_sld(Hup1,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        // PROLOGUE 10
        // prefetch next residue
        k                = db_sequence[10];
        k8               = db_sequence[2];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        v_score_load2 = vec_ld(16*k8,query_profile_byte);
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup2   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup2   = vec_sld(Hup2,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup2,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup1,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        
        // PROLOGUE 11
        // prefetch next residue
        k                = db_sequence[11];
        k8               = db_sequence[3];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        v_score_load2 = vec_ld(16*k8,query_profile_byte);
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1    = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup1    = vec_sld(Hup1,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        // PROLOGUE 12
        // prefetch next residue
        k                = db_sequence[12];
        k8               = db_sequence[4];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        v_score_load2 = vec_ld(16*k8,query_profile_byte);
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup2   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup2   = vec_sld(Hup2,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup2,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup1,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        
        // PROLOGUE 13
        // prefetch next residue
        k                = db_sequence[13];
        k8               = db_sequence[5];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        v_score_load2 = vec_ld(16*k8,query_profile_byte);
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1    = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup1    = vec_sld(Hup1,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        // PROLOGUE 14
        // prefetch next residue
        k                = db_sequence[14];
        k8               = db_sequence[6];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        v_score_load2 = vec_ld(16*k8,query_profile_byte);
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup2   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup2   = vec_sld(Hup2,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup2,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup1,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        // PROLOGUE 15
        // prefetch next residue
        k                = db_sequence[15];
        k8               = db_sequence[7];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        v_score_load2 = vec_ld(16*k8,query_profile_byte);
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup1    = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup1    = vec_sld(Hup1,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup1,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup2,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        
        
        // PROLOGUE 16
        // prefetch next residue
        k                = db_sequence[16];
        k8               = db_sequence[8];
        
        v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
        v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
        v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
        v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
        v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
        v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
        

        // prefetch score for next step 
        v_score_load1 = vec_ld(16*k,query_profile_byte);            
        v_score_load2 = vec_ld(16*k8,query_profile_byte);
        
        // load values of F and H from previous row (one unit up)
        Fup    = vec_ld(0,  p);
        Hup2   = vec_ld(16, p);
        p += 32; // move ahead 32 bytes
        
        // shift into place so we have complete F and H vectors
        // that refer to the values one unit up from each cell
        // that we are currently working on.
        Fup    = vec_sld(Fup,F,15);
        Hup2   = vec_sld(Hup2,H,15);            
        
        // do the dynamic programming 
        
        // update E value
        E   = vec_subs(E,v_gapextend);
        tmp = vec_subs(H,v_gapopen);
        E   = vec_max(E,tmp);
        
        // update F value
        F   = vec_subs(Fup,v_gapextend);
        tmp = vec_subs(Hup2,v_gapopen);
        F   = vec_max(F,tmp);
        
        v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
        
        // add score to H
        H   = vec_adds(Hup1,v_score);
        H   = vec_subs(H,v_bias);
        
        // set H to max of H,E,F
        H   = vec_max(H,E);
        H   = vec_max(H,F);
        
        // Update highest score encountered this far
        v_maxscore = vec_max(v_maxscore,H);
        
        p = workspace;
        
        for(j=16;j<db_length;j+=16)
        { 
            // STEP 1
            
            // prefetch next residue 
            k                = db_sequence[j+1];
            k8               = db_sequence[j-7];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
       
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup1   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup1    = vec_sld(Hup1,H,15);            

            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);

            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            

            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
          

            
            
            
            // STEP 2
            
            // prefetch next residue
            k                = db_sequence[j+2];
            k8               = db_sequence[j-6];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup2   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup2   = vec_sld(Hup2,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            

            
            
            
            // STEP 3
            
            // prefetch next residue
            k                = db_sequence[j+3];
            k8               = db_sequence[j-5];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup1   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup1    = vec_sld(Hup1,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
      
            

            
            
            // STEP 4
            
            // prefetch next residue
            k                = db_sequence[j+4];
            k8               = db_sequence[j-4];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup2   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup2   = vec_sld(Hup2,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            
            

            
            
            // STEP 5
            
            // prefetch next residue
            k                = db_sequence[j+5];
            k8               = db_sequence[j-3];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup1    = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup1   = vec_sld(Hup1,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            

            
            
            
            // STEP 6
            
            // prefetch next residue
            k                = db_sequence[j+6];
            k8               = db_sequence[j-2];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup2   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup2   = vec_sld(Hup2,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            

            
            
            
            // STEP 7
            
            // prefetch next residue
            k                = db_sequence[j+7];
            k8               = db_sequence[j-1];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup1    = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup1    = vec_sld(Hup1,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            
            

            
            
            // STEP 8
            
            // prefetch next residue
            k                = db_sequence[j+8];
            k8               = db_sequence[j];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup2   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup2   = vec_sld(Hup2,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            
            

            
            
            // STEP 9
            
            // prefetch next residue
            k                = db_sequence[j+9];
            k8               = db_sequence[j+1];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup1   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup1   = vec_sld(Hup1,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            // STEP 10
            
            // prefetch next residue
            k                = db_sequence[j+10];
            k8               = db_sequence[j+2];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup2   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup2   = vec_sld(Hup2,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
        
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            // STEP 11
            
            // prefetch next residue
            k                = db_sequence[j+11];
            k8               = db_sequence[j+3];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup1   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup1   = vec_sld(Hup1,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            // STEP 12
            
            // prefetch next residue
            k                = db_sequence[j+12];
            k8               = db_sequence[j+4];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup2   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup2   = vec_sld(Hup2,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            // STEP 13
            
            // prefetch next residue
            k                = db_sequence[j+13];
            k8               = db_sequence[j+5];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup1   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup1   = vec_sld(Hup1,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            // STEP 14
            
            // prefetch next residue
            k                = db_sequence[j+14];
            k8               = db_sequence[j+6];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup2   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup2   = vec_sld(Hup2,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            // STEP 15
            
            // prefetch next residue
            k                = db_sequence[j+15];
            k8               = db_sequence[j+7];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
                        
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup1   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup1   = vec_sld(Hup1,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
            // STEP 16
            
            // prefetch next residue
            k                = db_sequence[j+16];
            k8               = db_sequence[j+8];
            
            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load1 = vec_ld(16*k,query_profile_byte);
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            
            // load values of F and H from previous row (one unit up)
            Fup    = vec_ld(512, p);
            Hup2   = vec_ld(528, p);
            
            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32;
            
            // shift into place so we have complete F and H vectors
            // that refer to the values one unit up from each cell
            // that we are currently working on.
            Fup    = vec_sld(Fup,F,15);
            Hup2   = vec_sld(Hup2,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup2,v_gapopen);
            F   = vec_max(F,tmp);
            
            v_score_load1 = vec_perm(v_score_load1,v_score_load2,merge_score_load);
            
            // add score to H
            H   = vec_adds(Hup1,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
            
        }
        
        for(;j<db_length+15;j++)
        {
            k8               = db_sequence[j-7];

            v_score     = vec_perm(v_score_q1,  v_score_load1,  queue1_to_score);
            v_score_q1  = vec_perm(v_score_q2,  v_score_load1,  queue2_to_queue1);
            v_score_q2  = vec_perm(v_score_q3,  v_score_load1,  queue3_to_queue2);
            v_score_q3  = vec_perm(v_score_q4,  v_score_load1,  queue4_to_queue3);
            v_score_q4  = vec_perm(v_score_q5,  v_score_load1,  queue5_to_queue4);
            v_score_q5  = vec_perm(v_score_q5,  v_score_load1,  queue5_with_load);
            
            
            // prefetch scores for next step
            v_score_load2 = vec_ld(16*k8,query_profile_byte);
            v_score_load1 = vec_perm(v_zero,v_score_load2,merge_score_load);

            // save old values of F and H to use on next row
            vec_st(F, 0,  p);
            vec_st(H, 16, p);
            p += 32; // move ahead 32 bytes
            
            Fup    = vec_sld(v_zero,F,15);
            Hup1   = vec_sld(v_zero,H,15);            
            
            // do the dynamic programming 
            
            // update E value
            E   = vec_subs(E,v_gapextend);
            tmp = vec_subs(H,v_gapopen);
            E   = vec_max(E,tmp);
            
            // update F value
            F   = vec_subs(Fup,v_gapextend);
            tmp = vec_subs(Hup1,v_gapopen);
            F   = vec_max(F,tmp);
            
            // add score to H
            H   = vec_adds(Hup2,v_score);
            H   = vec_subs(H,v_bias);
            
            // set H to max of H,E,F
            H   = vec_max(H,E);
            H   = vec_max(H,F);
            
            // Save value to use for next diagonal H 
            Hup2 = Hup1;

            // Update highest score encountered this far
            v_maxscore = vec_max(v_maxscore,H);
        }
        vec_st(F, 512, p);
        vec_st(H, 528, p);

        query_profile_byte += 16*alphabet_size;

        // End of this row (actually 16 rows due to SIMD).
        // Before we continue, check for overflow.
        tmp      = vec_subs(vec_splat_u8(-1),v_bias);
        overflow = vec_any_ge(v_maxscore,tmp);
        

    }

    if(overflow)
    {
        return 255;
    }
    else
    {
        // find largest score in the v_maxscore vector
        tmp = vec_sld(v_maxscore,v_maxscore,8);
        v_maxscore = vec_max(v_maxscore,tmp);
        tmp = vec_sld(v_maxscore,v_maxscore,4);
        v_maxscore = vec_max(v_maxscore,tmp);
        tmp = vec_sld(v_maxscore,v_maxscore,2);
        v_maxscore = vec_max(v_maxscore,tmp);
        tmp = vec_sld(v_maxscore,v_maxscore,1);
        v_maxscore = vec_max(v_maxscore,tmp);
        
        // store in temporary variable
        vec_ste(v_maxscore,0,&score);
        
        // return largest score
        return score;
    }}


#else

/* No Altivec support. Avoid compiler complaints about empty object */

int sw_dummy;

#endif
