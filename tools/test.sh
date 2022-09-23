#!/bin/bash
#==============================================================================

function main ()
{
  scratch_dir=/gpfs/alpine/stf006/scratch/$USER

  local metric_type_prec
  for metric_type_prec in ccc duo czekanowski_single czekanowski_double ; do
  #for metric_type_prec in ccc duo ; do
  #for metric_type_prec in duo ; do
  #for metric_type_prec in czekanowski_single ; do

    local metric_type=$(echo $metric_type_prec | cut -f1 -d_)
    local prec=$(echo $metric_type_prec | cut -f2 -d_)

    # Make tped file.
    num_vector=17 num_field=73

    cmd="./make_tped_example.sh $metric_type ${metric_type}_test.tped $num_vector $num_field"
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
      cut -f5- ${metric_type}_test.tped \
        | awk 'BEGIN{FS=OFS=ORS=""}{split($0,a); asort(a); for(i in a) 
            if (a[i] != a[i-1] && a[i] != "\t") printf a[i] ; print "\n"}' \
        | sed -e 's/^.$/&&/' \
        > ${metric_type}_test_allele_labels.txt
    fi

    # Validate results against original tped file.

    cmd="./preprocess_validate $metric_type_prec \
                               ${metric_type_prec}_test.bin \
                               ${metric_type}_test.tped \
                               ${metric_type}_test_allele_labels.txt"
    echo $cmd
    $cmd

    for num_way in 2 3 ; do

      # CoMet run to generate metrics

      if [ $prec = single ] ; then
        exec=$HOME/genomics/build_single_test_summit/genomics_metric
      else
        exec=$HOME/genomics/build_test_summit/genomics_metric
      fi

      if [ $metric_type = czekanowski ] ; then
        tc=0
      else
        tc=1
      fi

      if [ 1 = 1 ] ; then

        cp ${metric_type_prec}_test.bin $scratch_dir

        cmd="jsrun --nrs 1 --rs_per_host 1 --cpu_per_rs 1 -g 1 --tasks_per_rs 1 \
          -X 1 --smpiargs='-gpu' $exec \
          --num_way $num_way --metric_type $metric_type --sparse yes \
          --num_vector $num_vector --num_field $num_field --all2all yes --compute_method GPU \
          --verbosity 1 --checksum yes --tc $tc \
          --input_file $scratch_dir/${metric_type_prec}_test.bin \
          --output_file_stub $scratch_dir/${metric_type_prec}_test_${num_way}way"
        echo $cmd
        $cmd

        cp $scratch_dir/${metric_type_prec}_test_${num_way}way_0.bin .

      fi

      cmd="./postprocess $num_way ${metric_type}_test_allele_labels.txt \
        ${metric_type}_test_line_labels.txt \
        ${metric_type_prec}_test_${num_way}way_0.bin \
        ${metric_type_prec}_test_${num_way}way_0.txt"
      echo $cmd
      $cmd

      tr ' ' '\t' < ${metric_type_prec}_test_${num_way}way_0.txt | \
        cut -f1-$(( 2 * $num_way )) | \
        ./validate $metric_type_prec $num_way ${metric_type}_test.tped \
          ${metric_type}_test_line_indices.bin \
          > ${metric_type_prec}_test_${num_way}way_0_validate.txt

      cmd="./validate_all.sh $num_way ${metric_type_prec}_test_${num_way}way_0.txt"
      echo $cmd
      $cmd

    done # num_way

    if [ 0 = 1 ] ; then

      # Cleanup.

      rm -f ${metric_type_prec}_test.bin \
            ${metric_type}_test_line_labels.txt \
            ${metric_type}_test_line_indices.bin \
            ${metric_type}_test_allele_labels.txt

    fi

  done # metric_type_prec

} # main

#------------------------------------------------------------------------------

main $@

#==============================================================================

