#!/usr/bin/env Rscript --vanilla

################################################################
# copyright (c) 2022 by William R. Pearson and The Rector &
# Visitors of the University of Virginia */
################################################################
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under this License is distributed on an "AS
# IS" BASIS, WITHOUT WRRANTIES OR CONDITIONS OF ANY KIND, either
# express or implied.  See the License for the specific language
# governing permissions and limitations under the License.
################################################################


## m8CBl_to_plot2.R alignment.m8CBl_file
##
## take an m8CBl output of the form:
##
## # lalign36 -m 8CBl ../seq/mchu.aa ../seq/mchu.aa
## # LALIGN 36.3.8h Aug, 2019
## # Query: sp|P62158|CALM_HUMAN Calmodulin; CaM - 149 aa
## # Database: ../seq/mchu.aa
## # Fields: query id, query length, subject id, subject length, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, BTOP
## # 1 hits found
## sp|P62158|CALM_HUMAN	149	sp|P62158|CALM_HUMAN	149	47.95	73	38	3	1	76	77	149	3.1e-11	49.3	1AK1QTLDTS2QE1A-E-F-KR3SRLV5DN1TY1TSTAKA2GRTH2RTSN2QENKPL1EDAE1LVQDDE2NR1VA1AI2ND1TQIVDNFYPE2LVTQ2ATRA1
## sp|P62158|CALM_HUMAN	149	sp|P62158|CALM_HUMAN	149	37.11	97	61	8	11	111	47	147	6.7e-09	41.6	2FLKQEDAMFISNLEFV1KA2DN3TDTFKP1LFGL1VM1-A1SKLMGKQDNTPDTS1AE1LIQRDEMAIFNREVVF1AK4TY1DSFAPA1FL-R-H-VLM1MNMLAGRE1MLKT1TEDESVED1EM4F-R-V-F-1KI2ND1YQIVSNAYAE1LFRVHQVM2
## sp|P62158|CALM_HUMAN	149	sp|P62158|CALM_HUMAN	149	39.39	33	20	5	1	37	113	146	0.065	18.3	MLAGDEQK2ED1QEIVAD1FM-IKR2F-S-L-F-1KI4TQIVTNTYKE1LFGVTQVM1
## # LALIGN processed 1 queries

## and convert to an x-y plot using ggplot

################
## This is a second attempt to provide an alternative to
## lav2_plt.pl, which takes an lalign.lav file and converts it to
## postscript or SVG.  Unlike m8CBl_to_plot.R, this version uses more vector math and fewer for() loops
##
## This program does not use the .lav file -- it used the BTOP alignment encoding provided by lalign36 -m8CBl
## ( the 'l' in -m8CBl is required to get thte sequence lengths)
##
################
##
## future versions could include expectation values on the alignments
## and domain diagrams
##

library(ggplot2)
library(RColorBrewer)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)

p.name<-'m8CBl_to_plot2.R'
plabel=paste(c(p.name,args,date()),collapse=' ', sep=' ')

if (is.na(args[1])) {
    print("usage m8CBl_to_plot.R alignment.m8CBl_file")
    exit(1)
}

btop_split <- function(btop_str) {

    x = c()
    y = c()

    ## print(btop_str)

    parts <- strsplit(gsub("([:digit:]+)","~\\1~",btop_str),"~")[[1]]
    ## print(parts)

    parts2 <- strsplit(gsub("([A-Z\\-][A-Z\\-])","~\\1~",parts),"~")
    parts2d <- lapply(parts2, function(z){ z[!is.na(z) & z != ""]})
    partsv <- Reduce(c,parts2d)
    ## print(partsv)
    partsv
}

btop2inc <- function(btop_vec) {

    ## convert vector to numbers if it has numbers

    btop_len = length(btop_vec)
    ## print(btop_vec)

    xy_inc = ifelse(grepl('[[:digit:]]+',btop_vec), strtoi(btop_vec), 0)
    xy_mat = ifelse(grepl('[A-Z][A-Z]',btop_vec), 1, 0)

    x_inc  = ifelse(substr(btop_vec,2,2)=='-', 1, 0)
    y_inc  = ifelse(substr(btop_vec,1,1)=='-', 1, 0)

    x_inc = x_inc + xy_inc + xy_mat
    y_inc = y_inc + xy_inc + xy_mat

    xy_list = list(x=x_inc, y=y_inc)
    ## print(xy_list)
    xy_list
}

btop2xy <- function(btop_list) {

    this_x = 0
    this_y = 0

    x_list = c()
    y_list = c()
    ## print(length(btop_list))
    ## print(btop_list)

    xy_inc = btop2inc(btop_list)

    x_list = cumsum(xy_inc$x)
    y_list = cumsum(xy_inc$y)

    xy_list = list(x=x_list, y=y_list)
    ## print(xy_list)
    xy_list
}

eval_type <- function(eval, thresh, pref) {

    thresh_len = length(thresh)
    for (ix in 1:thresh_len) {
        if (eval < thresh[ix]) {
            etype = ix
            break
        }
    }

    if (eval >= thresh[thresh_len]) {
        etype = thresh_len+1
    }

    return(etype)
}

## end of functions // main program

elinval = c(1e-4, 1e-2, 1.0, 10.0)
blinval = c(40.0,30.0,20.0, 10.0)
line_pref = 'l'

m8_col_names = c("qseqid","qlen","sseqid","slen","percid","aln_len","mismat","gaps","q_start","q_end","s_start","s_end","evalue","bits","BTOP")
m8_col_names_ann = c("qseqid","qlen","sseqid","slen","percid","aln_len","mismat","gaps","q_start","q_end","s_start","s_end","evalue","bits","BTOP","ann_str")
m8_col_names_ann2 = c("qseqid","qlen","sseqid","slen","percid","aln_len","mismat","gaps","q_start","q_end","s_start","s_end","evalue","bits","BTOP","ann_str","seq_ann_str")

m8_df <- read.table(args[1],header=FALSE)

m8_df_ncols <- dim(m8_df)[2]

if (m8_df_ncols == length(m8_col_names)) {
    colnames(m8_df) = m8_col_names
} else if (m8_df_ncols == length(m8_col_names_ann)) {
    colnames(m8_df) = m8_col_names_ann
} else if (m8_df_ncols == length(m8_col_names_ann2)) {
    colnames(m8_df) = m8_col_names_ann2
} else {
    stop("ERROR --  do not recognize number of input columns")
}

m8_df$row_num = seq.int(nrow(m8_df))

## identical sequences, plot mirror alignment
do_mirror = (m8_df[1,]$qseqid == m8_df[1,]$sseqid)

## take the BTOP encoding and convert it to a vector of tokens

btop_list <- lapply(m8_df$BTOP, btop_split)

## take the list of token vectors (one for each alignment) and convert to xy-offsets
xy_list <- lapply(btop_list, btop2xy)

## create a vector of names for each alignment -- required to keep x,y pairs separate for each alignment
diag_name = paste0(line_pref,"0")

## keep a list of the alignment names used
lname_list = c(diag_name)

## line_colors <- brewer.pal(length(elinval)+1,'Dark2')
line_colors <- c('darkblue','brown','darkgreen','green','lightgreen')

## line_color for identity  alignemnt
line_colors_n <- c()

## create x,y,aln_name dataframe for plotting, including diagonal if do_mirror
if (do_mirror) {
    xy_df = data.frame(name='aln0',x=c(1,m8_df[1,]$qlen), y=c(1,m8_df[1,]$slen),ev_ix=0, u.eline=diag_name)
    xy_lab = data.frame(u.eline=diag_name,lab='',x=m8_df[1,]$qlen/2.0,y=m8_df[1,]$slen/2.0)
    line_colors_n <- c('black')
} else {
    xy_df = NULL
    xy_lab = NULL
}

## go through each of the alignments produced by btop_split/btop2xy -- add to dataframe for plotting
for (ix in 1:length(xy_list)) {

    q_off = m8_df[ix,]$q_start-1
    s_off = m8_df[ix,]$s_start-1

    ## data frame for this alignment
    tmp_df = data.frame(name=paste0('aln',ix), x=xy_list[[ix]]$x, y=xy_list[[ix]]$y)
    tmp_df$x = tmp_df$x + q_off
    tmp_df$y = tmp_df$y + s_off

    ## add linetype/alignment color
    ev_ix = eval_type(m8_df[ix,]$evalue, elinval)
    ## print(sprintf("ix: %d evalue: %g ev_ix: %d",ix, m8_df[ix,]$evalue, ev_ix))

    tmp_df$ev_ix = ev_ix
    aln_name = paste0(line_pref,ev_ix,"_",ix)
    tmp_df$u.eline= aln_name
    line_colors_n <- append(line_colors_n,line_colors[ev_ix])

    ## add this alignment to the alignemnt dataframe
    xy_df = rbind(xy_df, tmp_df)

    q_mid = (m8_df[ix,]$q_end - m8_df[ix,]$q_start+1)/2.0 + m8_df[ix,]$q_start-1
    s_mid = (m8_df[ix,]$s_end - m8_df[ix,]$s_start+1)/2.0 + m8_df[ix,]$s_start-1

    tmp_lab = data.frame(u.eline=aln_name,lab=sprintf("%0.2g",m8_df[ix,]$evalue),x=q_mid,y=s_mid)

    xy_lab = rbind(xy_lab, tmp_lab)

    if (do_mirror) {
        tmp_df = data.frame(name=paste0('aln',ix,'c'), y=xy_list[[ix]]$x, x=xy_list[[ix]]$y)
        tmp_df$x = tmp_df$x + s_off
        tmp_df$y = tmp_df$y + q_off
        tmp_df$ev_ix = ev_ix
        tmp_df$u.eline= paste0(line_pref,ev_ix,"_",ix,'c')
        xy_df = rbind(xy_df, tmp_df)

        line_colors_n <- append(line_colors_n,line_colors[ev_ix])
    }
}

x_label = sprintf("%s %d aa",m8_df[1,]$qseqid, m8_df[1,]$qlen)
y_label = sprintf("%s %d aa",m8_df[1,]$sseqid, m8_df[1,]$slen)

pdf(file=paste0(args[1],'.pdf'),width=5, height=6)

## set colors based on evalue's
e_colors = scale_color_manual(values=line_colors_n)

## print(xy_lab)

elab_size = 4
elab_size = 4 - (nrow(xy_lab)-4)/4 * 0.75
if (elab_size < 2) {
    elab_size = 2
}

## pA -- draw alignment plot

pA <- ggplot(xy_df) + geom_line(aes(x=x,y=y,color=u.eline)) + geom_text(data=xy_lab, aes(x=x,y=y,label=lab),angle=45,nudge_x= -2.0, nudge_y=2.0, hjust=0.5, vjust=0,size=elab_size) +
    coord_fixed(ratio=1) + xlab(x_label) + ylab(y_label) + e_colors  +
    theme_bw() + guides(color='none') + theme(plot.margin=unit(c(2.0, 0.5, 0.5, 0.5), "cm"))

## pB -- draw evalue color legend
xlab_pos = c(5, 45, 5, 45, 80)
ylab_pos = c(15, 15, 10, 10, 12.5)
l_types = c('l1','l2','l3','l4','l5')
lab_df <- data.frame(x1 = xlab_pos, x2=xlab_pos + 15, y1 = ylab_pos, y2=ylab_pos, u.eline=l_types, etext=c("< 0.0001", "< 0.01","< 1.0", "< 10","> 10"))

e_colors2 = scale_color_manual(values=line_colors)

pB <- ggplot(lab_df) + geom_segment(aes(x=x1, xend=x2, y=y1,yend=y2, color=u.eline)) + geom_text(aes(x=x2+2,y=y1+0.03, label=etext),hjust=0) + xlim(c(-20,120)) + ylim(c(5,18)) + coord_fixed() + e_colors2 + guides(color='none') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),plot.margin=unit(c(-0.5, 0.5, 0.5, 0.5), "cm")) + labs(caption=plabel)

## plot out pA (the data) above pB (the legend)
plot_grid(pA,pB,ncol=1,rel_heights=c(0.8,0.2),rel_widths=c(1.0,1.0))

