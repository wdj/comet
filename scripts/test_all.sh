#!/bin/bash
#==============================================================================
#
# Run code tests on ALL versions.
# This should be executed in the directory containing the genomics_gpu
# repo directory.
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

# qsub -Astf006 -lnodes=4 -lwalltime=2:0:0 -I
# bsub -P stf006 -Is -nnodes 2 -alloc_flags gpumps -W 120 $SHELL

#==============================================================================

function script_dir
{
  local TARGET_FILE="$0"

  if [ $(uname -s) != "Darwin" ] ; then
    echo $(dirname "$(readlink -f "$TARGET_FILE")")
    return
  fi

  cd $(dirname $TARGET_FILE)
  TARGET_FILE=$(basename $TARGET_FILE)
  # Iterate down a (possible) chain of symlinks
  while [ -L "$TARGET_FILE" ] ; do
    TARGET_FILE=$(readlink $TARGET_FILE)
    cd $(dirname $TARGET_FILE)
    TARGET_FILE=$(basename $TARGET_FILE)
  done
  # Compute the canonicalized name by finding the physical path
  # for the directory we're in and appending the target file.
  local PHYS_DIR=$(pwd -P)
  local RESULT=$PHYS_DIR/$TARGET_FILE
  echo $(dirname $RESULT)
}

#------------------------------------------------------------------------------

function repo_dir
{
  echo "$(script_dir)/.."
}

#==============================================================================

function main
{
  # Location of this script.
  #local SCRIPT_DIR=$(script_dir)
  local REPO_DIR="${COMET_REPO_DIR:-$(repo_dir)}"
  local SCRIPT_DIR="$REPO_DIR/scripts"
  # Perform initializations pertaining to platform of build.
  . $SCRIPT_DIR/_platform_init.sh

  # Pick up builds directory.
  local BUILDS_DIR="${COMET_BUILDS_DIR:-$PWD}"

  local DO_SINGLE=1
  local DO_DOUBLE=1

  while [ "${1:-}" != "" ] ; do
    case $1 in
      --nosingle) DO_SINGLE=0 ;;
      --nodouble) DO_DOUBLE=0 ;;
      *)          echo "${0##*/}: Unrecognized argumnt. $1" 1>&2 ; exit 1 ;;
    esac
    shift
  done

  local DIRS=""

  if [ $DO_SINGLE = 1 ] ; then
    if [ -d "$BUILDS_DIR/build_single_test_$COMET_PLATFORM_STUB" ] ; then
      DIRS+=" $BUILDS_DIR/build_single_test_$COMET_PLATFORM_STUB"
    else
      DIRS+=" $BUILDS_DIR/build_single_test_nompi_$COMET_PLATFORM_STUB"
    fi
  fi

  if [ $DO_DOUBLE = 1 ] ; then
    if [ -d "$BUILDS_DIR/build_test_$COMET_PLATFORM_STUB" ] ; then
      DIRS+=" $BUILDS_DIR/build_test_$COMET_PLATFORM_STUB"
    else
      DIRS+=" $BUILDS_DIR/build_test_nompi_$COMET_PLATFORM_STUB"
    fi
  fi

  # Do tests.

  local DIR
  for DIR in $DIRS ; do
    printf -- '-%.0s' {1..79}; echo ""
    pushd $DIR
    time make test ARGS=-V 2>&1 | tee out_test.txt
    if [ $? != 0 ] ; then
      exit $?
    fi
    popd
    printf -- '-%.0s' {1..79}; echo ""
  done

  # Final reporting.

  printf -- '-%.0s' {1..79}; echo ""

  for DIR in $DIRS ; do
    grep -H fail $DIR/out_test.txt
  done
  OUT_FILES="$(for DIR in $DIRS ; do echo $DIR/out_test.txt ; done)"
  if [ $(awk '/ tests fail/' $OUT_FILES | wc -l) = \
       $(awk '/ 0 tests fail/' $OUT_FILES | wc -l) ] ; then
    echo "!!! All tests PASSED !!!"
  else
    echo "!!! Some tests FAILED !!!"
  fi

  printf -- '-%.0s' {1..79}; echo ""
} # main

#==============================================================================

main "$@"

#==============================================================================
