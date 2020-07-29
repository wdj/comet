#!/bin/bash
#==============================================================================
#
# Build ALL versions.
# This should be executed in the directory containing the genomics_gpu
# repo directory.
#==============================================================================

# Exit immediately on error.
set -eu -o pipefail

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

#==============================================================================

function main
{
  # Location of this script.
  local SCRIPT_DIR=$(script_dir)
  # Perform initializations pertaining to platform of build.
  . $SCRIPT_DIR/_platform_init.sh

  local DO_TEST=1
  local DO_RELEASE=1
  local DO_MPI=1
  local DO_NOMPI=1
  local DO_SINGLE=1
  local DO_DOUBLE=1

  while [ "${1:-}" != "" ] ; do
    case $1 in
      --notest)     DO_TEST=0 ;;
      --norelease)  DO_RELEASE=0 ;;
      --nompi)      DO_MPI=0 ;;
      --nonompi)    DO_NOMPI=0 ;;
      --nosingle)   DO_SINGLE=0 ;;
      --nodouble)   DO_DOUBLE=0 ;;
      *)            echo "${0##*/}: Unrecognized argumnt. $1" 1>&2 ; exit 1 ;;
    esac
    shift
  done

  local DIRS=""

  if [ $DO_TEST = 1 ] ; then
    if [ $DO_MPI = 1 ] ; then
      if [ $DO_SINGLE = 1 ] ; then
        DIRS+=" build_single_test_$COMET_PLATFORM_STUB"
      fi
      if [ $DO_DOUBLE = 1 ] ; then
        DIRS+=" build_test_$COMET_PLATFORM_STUB"
      fi
    fi
    if [ $DO_NOMPI = 1 ] ; then
      if [ $DO_SINGLE = 1 ] ; then
        DIRS+=" build_single_test_nompi_$COMET_PLATFORM_STUB"
      fi
      if [ $DO_DOUBLE = 1 ] ; then
        DIRS+=" build_test_nompi_$COMET_PLATFORM_STUB"
      fi
    fi
  fi

  if [ $DO_RELEASE = 1 ] ; then
    if [ $DO_MPI = 1 ] ; then
      if [ $DO_SINGLE = 1 ] ; then
        DIRS+=" build_single_release_$COMET_PLATFORM_STUB"
      fi
      if [ $DO_DOUBLE = 1 ] ; then
        DIRS+=" build_release_$COMET_PLATFORM_STUB"
      fi
    fi
    if [ $DO_NOMPI = 1 ] ; then
      if [ $DO_SINGLE = 1 ] ; then
        DIRS+=" build_single_release_nompi_$COMET_PLATFORM_STUB"
      fi
      if [ $DO_DOUBLE = 1 ] ; then
        DIRS+=" build_release_nompi_$COMET_PLATFORM_STUB"
      fi
    fi
  fi

  local DIR
  for DIR in $DIRS ; do
    printf -- '-%.0s' {1..79}; echo ""
    pushd $DIR
    $SCRIPT_DIR/make.sh 2>&1 | tee out_make.txt
    local MYSTATUS=$?
    if [ $MYSTATUS != 0 ] ; then
      echo "Build failure." 1>&2
      exit $MYSTATUS
    fi
    popd
    printf -- '-%.0s' {1..79}; echo ""
  done
}

#==============================================================================

main "$@"

#==============================================================================
