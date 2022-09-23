#!/bin/bash
#==============================================================================

function random_max ()
{
  #local result=32767
  local result=714025
  echo $result
}

function random ()
{
  local r="$1"
  #local result=$RANDOM
  local result=$(( ( $r * 4096 + 150889 ) % 714025 ))
  echo $result
}

#------------------------------------------------------------------------------

function make_tped ()
{
  local metric_type="$1"
  local num_vector="$2"
  local num_field="$3"

  cp /dev/null "$output_file"
  local DELIM=$'\t'

  # random seed.
  local r=0

  local vector_num
  for vector_num in $(seq 1 $num_vector) ; do

    local c1=A
    local c2=G

    r=$(random $r)
    if [ $r -lt $(( $(random_max) / 2 )) ] ; then
      c1=C
      c2=T
    fi

    r=$(random $r)
    echo -n "col1"
    echo -n "${DELIM}"
    echo -n "label-${vector_num}-$r"
    echo -n "${DELIM}"
    echo -n "col3"
    echo -n "${DELIM}"
    echo -n "col4"

    local field_num
    for field_num in $(seq 1 $num_field) ; do

      echo -n "${DELIM}"

      #-----
      if [ "$metric_type" = czekanowski ] ; then 
      #-----

        r=$(random $r)
        perl -e "print $r * 5. / $(random_max) ;"

      #-----
      elif [ "$metric_type" = ccc ] ; then 
      #-----

        r=$(random $r)
        if [ $r -lt $(( $(random_max) / 10 )) ] ; then

          echo -n "0"
          echo -n "${DELIM}"
          echo -n "0"

        else

          r=$(random $r)
          if [ $r -lt $(( $(random_max) / 2 )) ] ; then
            echo -n "$c1"
          else
            echo -n "$c2"
          fi

          echo -n "${DELIM}"

          r=$(random $r)
          if [ $r -lt $(( $(random_max) / 2 )) ] ; then
            echo -n "$c1"
          else
            echo -n "$c2"
          fi

        fi

      #-----
      elif [ "$metric_type" = duo ] ; then 
      #-----

        r=$(random $r)
        if [ $r -lt $(( $(random_max) / 10 )) ] ; then
          echo -n "0"
        else
          r=$(random $r)
          if [ $r -lt $(( $(random_max) / 2 )) ] ; then
            echo -n "$c1"
          else
            echo -n "$c2"
          fi
        fi

      #-----
      fi
      #-----

    done # num_field

    echo

  done # num_vector


} # make_tped

#------------------------------------------------------------------------------

function main ()
{

  if [ -z "$1" ] ; then
cat <<EOF
${0##*/}: create a tped file for use in testing.
Usage: ${0##*/} <metric_type> <output_file> [<num_vector>] [<num_field>]
  where <metric_type> is czekanowski, ccc or duo.
EOF
    return
  fi

  local metric_type="$1"
  local output_file="$2"

  local num_vector="${3:-17}"
  local num_field="${4:-73}"

  if [ "$metric_type" != czekanowski -a "$metric_type" != ccc -a \
       "$metric_type" != duo ] ; then
    echo "Error: invalid metric_type. $1" 1>&2
    exit 1
  fi

  make_tped "$metric_type" "$num_vector" "$num_field" > "$output_file"

} # main

#------------------------------------------------------------------------------

main $@

#==============================================================================
