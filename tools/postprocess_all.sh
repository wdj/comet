#!/bin/bash
#==============================================================================
# Postprocess the binary metrics files to make text files, possibly
# invoking parallel nodes with mpirun to make faster.
#==============================================================================

function process_files_simple
{
  local num_way="$1"
  shift
  local allele_label_file="$1"
  shift
  local line_label_file="$1"
  shift
  local metrics_bin_files="$*"

  local metrics_bin_file
  for metrics_bin_file in $metrics_bin_files ; do

    if [ ! -e $metrics_bin_file ] ; then
      echo "Error: file $metrics_bin_file does not exist." 1>&2
      exit 1
    fi

    local metrics_txt_file
    metrics_txt_file=$(echo $metrics_bin_file | sed -e 's/.bin$/.txt/')

    if [ "$metrics_bin_file" = "$metrics_txt_file" ] ; then
      echo "Error: invalid filename. $metrics_txt_file" 1>&2
      exit 1
    fi

    postprocess_file $num_way $allele_label_file $line_label_file \
      $metrics_bin_file $metrics_txt_file

  done
}

#------------------------------------------------------------------------------

function main
{
  if [ "$*" = "" ] ; then
    echo -n "Usage: ${0##*/} <num_way> <allele_label_file> <line_label_file> "
    echo "<metrics_bin_file> ..."
    exit
  fi

  local num_way="$1"
  if [ "$num_way" -ne 2 -a "$num_way" -ne 3 ] ; then
    echo "Error: invalid value for num_way. $num_way" 1>&2
    exit 1
  fi

  shift
  local files_spec="$*"

  if [ -n "$PBS_NUM_NODES" ] ; then # If on rhea and in a batch job ...
    if [ -z "$OMPI_COMM_WORLD_SIZE" ] ; then # if not invoked by mpirun ...
      ppn=16
      tmpfile=tmp_$$
      # Store list of files in file because mpirun breaks command line breaks
      # if too many
      echo "$files_spec" > $tmpfile
      # Invoke mpirun, on this bash script
      mpirun -np $(( $PBS_NUM_NODES * $ppn )) --npernode $ppn $0 $num_way $tmpfile
      rm $tmpfile
    else # if invoked by mpirun ...
      files="$(cat $files_spec)" # retrieve the list of files
      # Let each mpi rank own a subset of the files
      local files_thisrank
      files_thisrank="$(echo $files \
        | tr ' ' '\12' \
        | tail -n +3 \
        | awk 'NR%'${OMPI_COMM_WORLD_SIZE}'=='${OMPI_COMM_WORLD_RANK})"
      local allele_label_file
      allele_label_file=$(echo $files \
        | tr ' ' '\12' \
        | sed -n -e '1p')
      local line_label_file
      line_label_file=$(echo $files \
        | tr ' ' '\12' \
        | sed -n -e '2p')
      # Process the files serially
      process_files_simple $num_way $allele_label_file $line_label_file \
        $files_thisrank
    fi
  else # Not a rhea batch job ...
    # Process the files serially
    process_files_simple $num_way $files_spec
  fi
}

#------------------------------------------------------------------------------

main $@

#==============================================================================
