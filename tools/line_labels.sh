#!/bin/bash
#==============================================================================
# Extract line labels from SNP text file, pad with blanks to uniform length.
#==============================================================================

function main
{
  if [ "$2" = "" ] ; then
    echo "${0##*/}: extract line labels from SNP text file"
    echo "Usage: ${0##*/} <snp_text_file> <label_file>"
    exit
  fi

  local infile="$1"
  local outfile="$2"

  local MAX_LABEL_LEN=$(cat $infile | \
    tr ' ' '\t' | \
    cut -f2 | \
    awk '{print length($0)}' | \
    sort -n | \
    tail -n1)

  #echo $MAX_LABEL_LEN

  tr ' ' '\t' < "$infile" | \
    cut -f2 | \
    awk '{printf("%-'$MAX_LABEL_LEN's\n", $0)}' > \
    "$outfile"

  # Check.

  local dots
  dots=$(for i in $(seq 1 $MAX_LABEL_LEN) ; do echo -n "."; done)

  local num_long_lines
  num_long_lines=$(grep ".$dots" "$outfile" | wc -l)

  if [ $num_long_lines -gt 0 ] ; then
    echo "Error: labels are too long; please adjust code." 1>&2
    exit 1
  fi

}

#------------------------------------------------------------------------------

main $@

#==============================================================================
