#!/bin/bash

cmd="";
for i in "$@"
do
    case $i in
	-o=*|--outname=*)
	    OUTNAME="${i#*=}"
	    shift # past argument=value
	    ;;
	*)
	    cmd="$cmd $i"
	    ;;
    esac
done

bl_asn=${OUTNAME}.asn
bl0_out="$OUTNAME.html"
blm_out="$OUTNAME.msa"
blt_out="$OUTNAME.bl_tab"

# echo "OUTFILE = ${OUTNAME}"

#export BLAST_PATH="/ebi/extserv/bin/ncbi-blast+/bin"
export BLAST_PATH="/seqprg/bin"

$BLAST_PATH/blastp -outfmt 11 $cmd > $bl_asn
$BLAST_PATH/blast_formatter -archive $bl_asn -outfmt 0 -html > $bl0_out
$BLAST_PATH/blast_formatter -archive $bl_asn -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score btop'  > $blt_out
$BLAST_PATH/blast_formatter -archive $bl_asn -outfmt 2  > $blm_out

