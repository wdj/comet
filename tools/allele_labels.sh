#!/bin/bash
#==============================================================================
# Extract the allele labels used in each line of the SNP text file.
# To conform to assumptions of other codes used here, for each line
# the labels in a given line are output in (ascending) alphabetical order.
#==============================================================================

function main
{
  if [ "$2" = "" ] ; then
    echo "${0##*/}: extract per-line allele labels from SNP text file"
    echo "Usage: ${0##*/} <snp_text_file> <allele_label_file>"
    exit
  fi

  local infile="$1"
  local outfile="$2"

#  cut -f5- < "$infile" \
#    | sed -e 's/^[T0\t]*$/tt/' -e 's/^[A0\t]*$/aa/' -e 's/^[G0\t]*$/gg/' -e 's/^[C0\t]*$/cc/' \
#          -e 's/^[CT0\t]*$/CT/' -e 's/^[AG0\t]*$/AG/' -e 's/^[AC0\t]*$/AC/' \
#          -e 's/^[CG0\t]*$/CG/' -e 's/^[AT0\t]*$/AT/' -e 's/^[GT0\t]*$/GT/' \
#    | tr 'a-z' 'A-Z' \
#    > "$outfile"

#  cut -f5- "$infile" \
#    | awk 'BEGIN{FS=OFS=ORS=""}{split($0,a); asort(a); for(i in a) 
#        if (a[i] != a[i-1] && a[i] != "\t") printf a[i] ; print "\n"}' \
#    | sed -e 's/^.$/&&/' \
#    > "$outfile"

  cut -f5- < "$infile" \
  | awk '
  function less_than(left, right) {
    return "" left <= "" right
  }
  function quicksort(data, left, right,   i, last)
  {
    if (left >= right)
      return

    quicksort_swap(data, left, int((left + right) / 2))
    last = left
    for (i = left + 1; i <= right; i++)
      if (less_than(data[i], data[left]))
        quicksort_swap(data, ++last, i)
    quicksort_swap(data, left, last)
    quicksort(data, left, last - 1)
    quicksort(data, last + 1, right)
  }
  function quicksort_swap(data, i, j,   temp)
  {
    temp = data[i]
    data[i] = data[j]
    data[j] = temp
  }
  BEGIN{FS=OFS=ORS=""} {
    n = split($0, a);
    quicksort(a, 0, n-1);
    for (i = 0; i < n; i++)  if (a[i] != a[i-1] && a[i] != "\t" && a[i] != "0") printf a[i] ; print "\n"
  }' \
  | sed -e 's/^.$/&&/' \
  > "$outfile"



}

#------------------------------------------------------------------------------

main $@

#==============================================================================
