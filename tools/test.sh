#!/bin/bash
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

#==============================================================================

function main ()
{
  local num_passed=0
  local num_failed=0

  #local num_vector_num_field="91,79 17,73"
  local num_vector_num_field
  for num_vector_num_field in 5,7 17,73 ; do

  local num_vector=$(echo $num_vector_num_field | cut -f1 -d,)
  local num_field=$(echo $num_vector_num_field | cut -f2 -d,)

  echo "num_vector $num_vector num_field $num_field"

  local metric_type_prec
  for metric_type_prec in ccc duo czekanowski_single czekanowski_double ; do
  #for metric_type_prec in ccc duo ; do
  #for metric_type_prec in duo ; do
  #for metric_type_prec in czekanowski_single ; do

    local metric_type=$(echo $metric_type_prec | cut -f1 -d_)
    local prec=$(echo $metric_type_prec | cut -f2 -d_)

    local cmd=""

    # Make tped file.

    cmd="./make_tped_example.sh $metric_type ${metric_type}_test.tped \
                                $num_vector $num_field"
    echo $cmd
    $cmd

    # Convert tped file to CoMet binary input file format.

    cmd="./preprocess $metric_type_prec \
                      ${metric_type}_test.tped \
                      ${metric_type_prec}_test.bin"
    echo $cmd
    $cmd

    # Extract line labels from tped.

    cmd="./line_labels.sh ${metric_type}_test.tped \
                          ${metric_type}_test_line_labels.txt"
    echo $cmd
    $cmd

    # Extract (binary) line indices for tped.

    cmd="./line_indices ${metric_type}_test.tped \
                        ${metric_type}_test_line_indices.bin"
    echo $cmd
    $cmd

    # Extract per-line allele labels from tped.

    if [ $metric_type != czekanowski ] ; then
      cmd="./allele_labels.sh ${metric_type}_test.tped \
                              ${metric_type}_test_allele_labels.txt"
      echo $cmd
      $cmd
    fi

    # Validate results against original tped file.

    cmd="./preprocess_validate $metric_type_prec \
                               ${metric_type_prec}_test.bin \
                               ${metric_type}_test.tped \
                               ${metric_type}_test_allele_labels.txt"
    echo $cmd
    $cmd

    for num_way in 2 3 ; do
    #for num_way in 3 ; do

      # CoMet run to generate metrics

      if [ 1 = 1 ] ; then

        if [ -e $PWD/../genomics_metric ] ; then
          exec=$PWD/../genomics_metric
        elif [ $prec = single ] ; then
          exec=$HOME/genomics/build_single_test_nompi_summit/genomics_metric
        else
          exec=$HOME/genomics/build_test_nompi_summit/genomics_metric
        fi

        if [ $metric_type_prec = czekanowski_single -a \
             $(basename $(dirname $(realpath $exec)) | grep single | wc -l) = 0 ] ; then
          continue
        fi

        if [ $metric_type_prec = czekanowski_double -a \
             $(basename $(dirname $(realpath $exec)) | grep single | wc -l) = 1 ] ; then
          continue
        fi

        launchcmd=""

        cmd="$launchcmd $exec \
          --num_way $num_way --metric_type $metric_type --sparse yes \
          --num_vector $num_vector --num_field $num_field \
          --all2all yes --compute_method CPU \
          --verbosity 1 --checksum yes \
          --input_file ${metric_type_prec}_test.bin \
          --output_file_stub ${metric_type_prec}_test_${num_way}way"
        echo $cmd
        $cmd

      fi

      cmd="./postprocess $metric_type $num_way \
        ${metric_type}_test_allele_labels.txt \
        ${metric_type}_test_line_labels.txt \
        ${metric_type_prec}_test_${num_way}way_0.bin \
        ${metric_type_prec}_test_${num_way}way_0.txt"
      echo $cmd
      $cmd

      cmd="tr ' ' '\t' < ${metric_type_prec}_test_${num_way}way_0.txt | \
        cut -f1-$(( 2 * $num_way )) | \
        ./validate $metric_type_prec $num_way ${metric_type}_test.tped \
          ${metric_type}_test_line_indices.bin \
          > ${metric_type_prec}_test_${num_way}way_0_validate.txt"
      echo $cmd
      bash -c "$cmd"

      cmd="./validate_all.sh $metric_type $num_way \
                             ${metric_type_prec}_test_${num_way}way_0.txt"
      echo $cmd
      set +e
      $cmd
      local errcode=$?
      set -e
      if [ $errcode = 0 ] ; then
        num_passed=$(( $num_passed + 1 ))
      else
        num_failed=$(( $num_failed + 1 ))
      fi

    done # num_way

    for num_way in 2 3 ; do

      if [ 1 = 1 ] ; then

        # Cleanup.

        rm -f ${metric_type}_test.tped \
              ${metric_type_prec}_test.bin \
              ${metric_type}_test_line_labels.txt \
              ${metric_type}_test_line_indices.bin \
              ${metric_type}_test_allele_labels.txt \
              ${metric_type_prec}_test_${num_way}way_0.bin \
              ${metric_type_prec}_test_${num_way}way_0.txt \
              ${metric_type_prec}_test_${num_way}way_0_validate.txt

      fi

    done

  done # metric_type_prec

  done # num_vector_num_field

  echo "TOTAL: $num_passed tests passed."
  echo "TOTAL: $num_failed tests failed."
} # main

#------------------------------------------------------------------------------

main $@

#==============================================================================

