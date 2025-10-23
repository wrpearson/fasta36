
CC=cc

all : pseg pmerge

pseg : pseg.c lnfac.h pgenwin.h pgenwin.o
	$(CC) -g -o pseg pseg.c pgenwin.o -lm

pmerge : pmerge.c pgenwin.h pgenwin.o
	$(CC) -g -o pmerge pmerge.c pgenwin.o -lm

pgenwin.o : pgenwin.c pgenwin.h
	$(CC) -g -c pgenwin.c

clean:
	rm -f pseg pseg.o pgenwin.o
