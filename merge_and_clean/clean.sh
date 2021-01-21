#!/bin/bash


# For merged_alldata.rds, need to set n_gene = 10
# For mesenchyme, need to set n_gene = 10

N_DROP=10
#N_DROP=$2

if [ "$#" -ne 2 ]
then
    echo "Needs two parameters: <input rds prefix> <n_workers>"
    exit
fi

echo "================================================="
echo "Starting step 1 - remove Gm42418"
echo "================================================="

if ! R --no-save -f ./pipeline/1_rm_gm42418.R --args $1 $N_DROP $2
then
    echo "Step 1 returned an error"
    exit
fi


echo "================================================="
echo "Starting step 2 - SCTransform"
echo "================================================="

if ! R --no-save -f ./pipeline/2_sct.R --args $1 $N_DROP $2
then
    echo "Step 2 returned an error"
    exit
fi

echo "================================================="
echo "Starting step 3 - split P7b"
echo "================================================="

if ! R --no-save -f ./pipeline/3_split_p7b.R --args $1 $N_DROP $2
then
    echo "Step 3 returned an error"
    exit
fi

echo "================================================="
echo "Starting step 4 - integrating P7b"
echo "================================================="

if ! R --no-save -f ./pipeline/4_integrate_p7b.R --args $1 $N_DROP $2
then
    echo "Step 4 returned an error"
    exit
fi

echo "================================================="
echo "Starting step 5 - rerun SCTransform after integration"
echo "================================================="

if ! R --no-save -f ./pipeline/5_reSCT.R --args $1 $N_DROP $2
then
    echo "Step 5 returned an error"
    exit
fi


echo "================================================="
echo "DONE"
echo "================================================="

