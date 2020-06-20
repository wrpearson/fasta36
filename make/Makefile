#
# makefile for fasta3, fasta3_t.  Use makefile.pvm for pvcompxx.
#
#
# Dec 8, 2005 - with gcc4.0.2 (or .1) under Redhat Linux Fedora FC4 -03 breaks the alignment code
#

CC ?= gcc -g -O2

#CC=gcc -Wall -pedantic -ansi -g -O
#CC = gcc -g -DDEBUG
#CC= /usr/local/parasoft/bin.linux2/insure -g -DDEBUG

# EBI uses the following with pgcc, -O3 does not work:
# CC= pgcc -O2 -pipe -mcpu=pentiumpro -march=pentiumpro -fomit-frame-pointer

# this file works for x86 LINUX

# standard options
CFLAGS += -DSHOWSIM -DUNIX -DTIMES -DHZ=100 -DSFCHAR="':'" -DMAX_WORKERS=2 -DTHR_EXIT=pthread_exit -DPROGRESS  -DUSE_MMAP -D_REENTRANT -DHAS_INTTYPES -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -DUSE_FSEEKO -DSAMP_STATS -DPGM_DOC
# -I/usr/local/include/mysql -DMYSQL_DB 
#
#(for mySQL databases)  (also requires change to Makefile34.common)

LIB_M = -lm
#LIB_M = -L/usr/local/lib/mysql -lmysqlclient -lm
# for mySQL databases

HFLAGS= -o
NFLAGS= -o

# for Linux
THR_SUBS = pthr_subs2
THR_LIBS = -lpthread
THR_CC =

BIN = ../bin
XDIR = /seqprg/bin

DROPGSW_NA_O = dropgsw2.o wm_align.o calcons_sw.o
DROPGSW_SSE_O = dropgsw2_sse.o smith_waterman_sse2.o wm_align.o calcons_sw.o
DROPGSW_ALT_O = dropgsw2_alt.o smith_waterman_altivec.o wm_align.o calcons_sw.o
DROPGSW_O = $(DROPGSW_SSE_O)

DROPLAL_NA_O = droplal2.o lsim4.o calcons_la.o
DROPLAL_SSE_O = droplal2_sse.o smith_waterman_sse2.o lsim4.o calcons_la.o
DROPLAL_ALT_O = droplal2_sse.o smith_waterman_altivec.o lsim4.o calcons_la.o
DROPLAL_O = $(DROPLAL_SSE_O)

# renamed (fasta36)  programs
include ../make/Makefile36m.common
# conventional (fasta3) names
# include ../make/Makefile.common

