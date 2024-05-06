#!/bin/bash

baseDir="."

mkdir -p "$baseDir"

#cumulativeOutputFile="${baseDir}/cumulative_SegmentationCkov2_output.log"

# Function to handle the interrupt
cleanup() {
    echo "Keyboard interrupt (Ctrl+C) detected. Exiting..."
    # Any cleanup code goes here. For example:
    echo "Performing cleanup tasks..."
    exit 1  # Exit script with a status code
}

# Setting up trap
trap cleanup SIGINT

for i in {1315..2000}; do
    echo "Iteration: $i"

    if /home/eflatlan/alice/O2/prodtests/sim_challengeCust.sh -f sim -n 30 -j 30 -r 50 -s pbpb; then
        if root -b -q -l "SegmentationCkov2.C(1.75)" > "SegmentationCkov2_output_${i}.log"; then
            echo "run ok"

            if cp data.root "dataSimLeadSigma40_${i}.root"; then
                echo "data.root moved"
            else 
                echo "moving data.root failed"      
            fi

        else
            echo "ROOT commands failed"
        fi

    else
        echo "sim_challengeCustom.sh failed"
    fi
done

