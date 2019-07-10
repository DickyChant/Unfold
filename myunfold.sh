#########################################################################
# File Name: run.sh
# Author: peng huo
# mail: penghuo12@gmail.com
# Created Time: Mon 21 Dec 2015 02:30:30 PM EST
#########################################################################
#!/bin/bash

root -b <<EOF
gSystem->Load("~/codes/RooUnfold/libRooUnfold");
.L myunfold.C+g
myunfold();
.q
EOF
