#!/bin/bash

cmd="";
for i in "$@"
do
    case $i in
	--outname=*)
	    OUTNAME="${i#*=}"
	    shift # past argument=value
	    ;;
	--query=*)
	    QUERY="${i#*=}"
	    shift # past argument=value
	    ;;
	--db=*)
	    DATABASE="${i#*=}"
	    shift # past argument=value
	    ;;
	--cmd=*)
	    SRCH_CMD="${i#*=}"
	    shift
	    ;;
	--ktup=*)
	    KTUP="${i#*=}"
	    shift
	    ;;
	*)
	    cmd="$cmd $i"
	    ;;
    esac
done


# echo "OUTNAME: " $OUTNAME
echo "# CMD: " $cmd

if [[ $OUTNAME == '' ]]; then
    OUTNAME=${QUERY}_out
fi

if [[ $SRCH_CMD == '' ]]; then
    SRCH_CMD=fasta36
fi

#if [[ $ANN_SCRIPT == '' ]]; then
#    ANN_SCRIPT="/seqprg/bin/ann_pfam30.pl --db=pfam31_qfo --host=localhost --neg --vdoms --acc_comment"
#fi


# echo "OUTNAME: " $OUTNAME

bl0_out="$OUTNAME.html"
bla_out="${OUTNAME}_an.html"
blt_out="$OUTNAME.fa_tab"
blr_out="$OUTNAME.fa_tab_rn"

export BLAST_PATH="/seqprg/bin"
# BLAST_PATH="../bin"

cmd="$cmd -mF8CBL=$blt_out $QUERY $DATABASE"

# echo "tmp_files:"
# echo $bl_asn $bl0_out $bla_out $blt_out
# echo "OUTFILE = ${OUTNAME}"

#echo "cmd: $cmd"
#echo "==="
#echo "bl0_out: $bl0_out"
#echo "==="

# echo "$BLAST_PATH/$SRCH_CMD $cmd > $bl0_out"

# run the program
$BLAST_PATH/$SRCH_CMD $cmd > $bl0_out

$BLAST_PATH/rename_exons.py --have_qslen --dom_info $blt_out > $blr_out

if [ ! -s $blr_out ]; then
    # echo "# " `ls -l $blt_out $blr_out`
    blr_out=$blt_out
    # echo "# " `ls -l $blt_out $blr_out`
fi

$BLAST_PATH/merge_fasta_btab.pl --plot_url="plot_domain6t.cgi" --have_qslen --dom_info --btab $blr_out $bl0_out
