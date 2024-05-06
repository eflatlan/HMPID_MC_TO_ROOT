#!/bin/bash
baseDir="."

mkdir -p "$baseDir"

#cumulativeOutputFile="${baseDir}/cumulative_SegmentationCkov2_output.log"

for i in {1..2000}; do
    echo "Iteration: $i"
    
    if ./sim_challengeCust.sh -f simk -s pbpb -n 10 -j 5 -r 50  > "sim${i}.log"; then
        if root -b -q -l "SegmentationCkov2.C(1.75)" > "SegmentationCkov2_output_${i}.log"; then
            echo "run ok"

          	if 	cp data.root "dataLead${i}.root"; then
          	echo "data.root moved"
		else 
      		echo "moving    data.root  failed"	
		fi
      else
            echo "ROOT commands failed"
        fi
    else
        echo "sim_challengeCustom.sh failed"
    fi
done

