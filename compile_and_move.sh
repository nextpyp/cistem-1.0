#!/bin/bash

make -j4
cd src
/bin/bash compile_refine3d.sh
/bin/bash compile_merge3d.sh
/bin/bash compile_reconstruct3d.sh
#/bin/bash compile_unblur.sh
#/bin/bash move_binaries.sh
