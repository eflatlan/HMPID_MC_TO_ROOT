### Run with O2 Branch : [massHypPC](https://github.com/eflatlan/AliceO2/tree/massHypPC)

# Run files with run_checker.sh
Run the file [run_checker](https://github.com/eflatlan/HMPID_MC_TO_ROOT/blob/officialMC/run_checker.sh)
This creates the simulation files in the current working directory

**run_checker.sh** runs the [sim_challengeCust.sh](https://github.com/eflatlan/HMPID_MC_TO_ROOT/blob/officialMC/sim_challengeCust.sh) executable in a loop. 
It is ran in a loop because it sometimes crashes when running many runs consecutively.

## sim_challengeCust
**sim_challengeCust.sh** is a modified version of the official [https://github.com/AliceO2Group/AliceO2/blob/dev/prodtests/sim_challenge.sh](sim_challenge.sh).
Here we have removed the AOD, and all devices not necessary for the HMPID.
Also input parameters are changed. 
- r : interaction rate
- f
  -  simk/simp/simpr : particleGun using pion/kaon/proton
  -  sim : standard collision
 
- s : collision system, only valid if the **-f** is used with the **sim** option
  - pp : proton-proton
  - pbpb : lead-lead
 
- j : number of cores to use
- number of events

## SegmentationCkov2
[https://github.com/eflatlan/HMPID_MC_TO_ROOT/blob/officialMC/SegmentationCkov2.C](SegmentationCkov2) is a macro that runs the post-processing of the HMPID data-output from the previous section.
It reads the clustefile, matchfile, and also gathers the MC-truth information.


## Run example

`alienv enter O2/latest-<BRANC_TO_USE>-o2`

 make them executable, only needed once : 

`chmod +x run_checker.sh`

`chmod +x sim_challengeCust.sh`

Run the executable :  

`./run_checker.sh`
# Merging output files
When running **run_checker.sh**, we will have as many output-files as specified in the loop in the file. 
We then can merge them by reading the TTree content in each file and append to the vectors.
This is done using the [mergeTTrees.C](https://github.com/eflatlan/HMPID_MC_TO_ROOT/blob/massHyp/mergeTTrees.C) macro. 
This macro reads all files containing a specific string and ending with **".root"** and creates a single output file.
(The string is specified on this line : [mergeTTrees.C#L27](https://github.com/eflatlan/HMPID_MC_TO_ROOT/blob/ab385fdde49f0ece5ae3f60f20bbd838ee5b816f/mergeTTrees.C#L27))
