#!/bin/bash

cmd="";
for i in "$@"
do
    case $i in
	-o=*|--outname=*)
	    OUTNAME="${i#*=}"
	    shift # past argument=value
	    ;;
	-q=*|--query=*)
	    QUERY="${i#*=}"
	    cmd="$cmd -query $QUERY"
	    shift # past argument=value
	    ;;
	--ann_script=*)
	    ANN_SCRIPT="${i#*=}"
	    shift
	    ;;
	--q_ann_script=*)
	    Q_ANN_SCRIPT="${i#*=}"
	    shift
	    ;;
	*)
	    cmd="$cmd $i"
	    ;;
    esac
done

# echo "OUTNAME: " $OUTNAME
# echo "CMD: " $cmd

if [[ $OUTNAME == '' ]]; then
    OUTNAME=${QUERY}_out
fi

#if [[ $ANN_SCRIPT == '' ]]; then
#    ANN_SCRIPT="/seqprg/bin/ann_pfam30.pl --db=pfam31_qfo --host=localhost --neg --vdoms --acc_comment"
#fi


# echo "OUTNAME2: " $OUTNAME

bl_asn="$OUTNAME.asn"
bl0_out="$OUTNAME.html"
bla_out="${OUTNAME}_an.html"
blm_out="$OUTNAME.msa"
blt_out="$OUTNAME.bl_tab"
blt_ann="$OUTNAME.bl_tab_ann"
blr_out="$OUTNAME.bl_tab_rn"

# echo "tmp_files:"
# echo $bl_asn $bl0_out $bla_out $blt_out

# echo "OUTFILE = ${OUTNAME}"

#export BLAST_PATH="/ebi/extserv/bin/ncbi-blast+/bin"
export BLAST_PATH="/seqprg/bin"

$BLAST_PATH/blastp -outfmt 11 $cmd > $bl_asn
$BLAST_PATH/blast_formatter -archive $bl_asn -outfmt 0 -html > $bl0_out
$BLAST_PATH/blast_formatter -archive $bl_asn -outfmt '7 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score btop' > $blt_out
annot_blast_btop2.pl --query $QUERY --have_qslen --dom_info --ann_script "$ANN_SCRIPT" --q_ann_script "$Q_ANN_SCRIPT" $blt_out > $blt_ann

rename_exons.py --have_qslen --dom_info $blt_ann > $blr_out
merge_blast_btab.pl --plot_url="plot_domain6t.cgi" --have_qslen --dom_info --btab $blr_out $bl0_out

# $BLAST_PATH/blast_formatter -archive $bl_asn -outfmt 2  > $blm_out
