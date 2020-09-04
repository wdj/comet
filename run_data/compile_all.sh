#!/bin/bash

echo ""
echo ""
echo "Compiling CoMet"
cd ${COMET_SRC}/comet_work_gcc8.3
${COMET_SRC}/scripts/make_all.sh --nompi --norelease --nosingle
echo "Done Compiling CoMet"

