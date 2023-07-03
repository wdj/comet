#!/bin/bash -l
#==============================================================================
#
# Build CoMet code.
#
# NOTE: it is recommended that this script not be called directly but
# instead that the script make_all.sh be used.
#
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

#------------------------------------------------------------------------------

function repo_dir
{
  echo "$(script_dir)/.."
}

#==============================================================================

function main
{
  # Initializations.

  # Location of this script.
  #local SCRIPT_DIR=$(script_dir)
  local REPO_DIR="${COMET_REPO_DIR:-$(repo_dir)}"
  local SCRIPT_DIR="$REPO_DIR/scripts"
  # Perform initializations pertaining to platform of build.
  . $SCRIPT_DIR/_platform_init.sh

  # make.

  time make -j4 VERBOSE=1
  #time make VERBOSE=1

  if [ $? != 0 ] ; then
    exit $?
  fi

  # make install.

  time make install

  if [ $? != 0 ] ; then
    exit $?
  fi

  # make and install the tools.

  (cd $REPO_DIR && tar cf - tools) | tar xf -
  pushd tools
  make
  for FILE in $(ls *.cc) ; do
    cp $(basename $FILE .cc) ../install_dir/bin
  done
  for FILE in $(ls *.sh) ; do
    cp $FILE ../install_dir/bin
  done
  popd

} # main

#==============================================================================

main "$@"

#==============================================================================
