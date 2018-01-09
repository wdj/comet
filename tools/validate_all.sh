#!/bin/bash
#==============================================================================

function process_files_simple
{
  local num_way="$1"
  shift
  local infiles="$*"

  local infile
  for infile in $infiles ; do

    if [ ! -e $infile ] ; then
      echo "Checking $infile ... Error: file does not exist."
      continue
    fi

    infilev=$(echo $infile | sed -e 's/.txt$/_validate&/')

    if [ ! -e $infilev ] ; then
      echo "Checking $infile ... Error: validation file does not exist."
      continue
    fi

    #echo "diff $F1 $F2"
    #diff <( sed -e 's/..$//' <$F1 ) <( sed -e 's/..$//' <$F2 )

    TOL=.00001
  
    F1=$infile
    F2=$infilev
    C1=$(( 3 * $num_way + 1 ))
    C2=$(( 6 * $num_way + 2 ))

    echo -n "Checking $F1 ... " >/dev/null
    local numdiffs
    numdiffs=$(paste $F1 $F2 | awk '$'$C1' - $'$C2' > '$TOL' || $'$C2' - $'$C1' > '$TOL' {print $0 }' | wc -l)
    if [ $numdiffs = 0 ] ; then
      echo "PASSED." >/dev/null >/dev/null
    else
      echo "FAILED with $numdiffs diffs." >/dev/null
      echo "Checking $F1 ... FAILED with $numdiffs diffs."
    fi

  done
}

#------------------------------------------------------------------------------

function main
{
  if [ "$*" = "" ] ; then
    echo "Usage: ${0##*/} <num_way> <file> ..."
    exit
  fi

  local num_way="$1"
  shift
  local infiles="$*"

  if [ -n "$PBS_NUM_NODES" ] ; then
    if [ -z "$OMPI_COMM_WORLD_SIZE" ] ; then
      ppn=16
      tmpfile=tmp_$$
      echo "$infiles" > tmp_$$
      mpirun -np $(( $PBS_NUM_NODES * $ppn )) --npernode $ppn ${0##*/} $num_way $tmpfile
      rm $tmpfile
    else
      infiles="$(cat $infiles)" # to avoid mpirun arg list too long error
      local infiles_this
      infiles_this="$(echo $infiles | tr ' ' '\12' | awk 'NR%'${OMPI_COMM_WORLD_SIZE}'=='${OMPI_COMM_WORLD_RANK})"
      process_files_simple $num_way $infiles_this
    fi
  else
    process_files_simple $num_way $infiles
  fi
}

#------------------------------------------------------------------------------

main $@

#==============================================================================
