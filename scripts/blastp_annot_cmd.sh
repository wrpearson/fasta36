#!/bin/bash

################
## blastp_annot_cmd.sh --query=seq/gstm1_human.aa --ann_script=scripts/ann_pfam_sql.pl --q_ann_script=scripts/ann_pfam_sql.pl -db /slib2/bl_dbs/pfam34_qfo78
##
## blastp reseults and intermediate files are saved in the same location as the query file
## annotation scripts must fully qualified
## --html=0 produces domain annotations but no domain graphics
## --html=1 provides exon mapping (if available) and domain graphics
##

cmd="";
DO_HTML=1

PGM='blastp'

for i in "$@"
do
    case $i in
	-p=*|--pgm=*)
	    PGM="${i#*=}"
	    shift # past argument=value
	    ;;
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
	--html=*)
	    DO_HTML="${i#*=}"
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
bla_out="${OUTNAME}_ann.html"
blm_out="$OUTNAME.msa"
blt_out="$OUTNAME.bl_tab"
blt_ann="$OUTNAME.bl_tab_ann"
blr_out="$OUTNAME.bl_tab_rn"

# echo "tmp_files:"
# echo $bl_asn $bl0_out $blt_out

# echo "OUTFILE = ${OUTNAME}"

#export BLAST_PATH="/ebi/extserv/bin/ncbi-blast+/bin"
export BLAST_PATH="/seqprg/bin"

$BLAST_PATH/$PGM -num_threads=8 -outfmt 11 $cmd > $bl_asn
$BLAST_PATH/blast_formatter -archive $bl_asn -outfmt '7 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore score btop' > $blt_out

# annot_cmd="annot_blast_btop2.pl --query $QUERY --raw --have_qslen --dom_info --ann_script "$ANN_SCRIPT" --q_ann_script "$Q_ANN_SCRIPT" $blt_out > $blt_ann"
# echo "# $annot_cmd"
annot_blast_btop2.pl --query $QUERY --raw --have_qslen --dom_info --ann_script "$ANN_SCRIPT" --q_ann_script "$Q_ANN_SCRIPT" $blt_out > $blt_ann

if [[ $DO_HTML == 1 ]]; then
    ## rename_exons.py --have_qslen --dom_info $blt_ann > $blr_out
    $BLAST_PATH/blast_formatter -archive $bl_asn -outfmt 0 -html > $bl0_out
    ## echo "# merge_blast_btab.pl --plot_url=plot_domain6t.cgi --have_qslen --dom_info --btab $blt_out $bl0_out"
    merge_blast_btab.pl --plot_url="plot_domain6t.cgi" --have_qslen --dom_info --btab $blt_ann $bl0_out

else
    $BLAST_PATH/blast_formatter -archive $bl_asn -outfmt 0 > $bl0_out
    merge_cmd="merge_blast_btab.pl --have_qslen --dom_info --btab $blt_ann $bl0_out"
    # echo "# $merge_cmd"
    merge_blast_btab.pl --have_qslen --dom_info --btab $blt_ann $bl0_out
fi

# $BLAST_PATH/blast_formatter -archive $bl_asn -outfmt 2  > $blm_out
