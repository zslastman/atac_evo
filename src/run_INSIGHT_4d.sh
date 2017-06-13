###A simple local script for running INSIGHT on TFBS data



screeningFile='/g/furlong/garfield/projects/INSIGHT_analyses/useful_annotations/ncRNA_and_exonic_masks/dmel-all-r5.57_genes_exons.bed'

inputBed=$1 #best to specify a full path
resultsFolder=$2 #best to specify a full path
deletePostSites=True #set to True at the moment...we almost never use these posterior sites
deleteLog=False ##set to False...the log file is the easiest place to read the results.
mainScript=/g/furlong/garfield/projects/INSIGHT_analyses/src/run_INSIGHT/run_INSIGHT.py
#freqCutOff=$3

#we can also code an alternative version for 4d sites
validSites=/g/furlong/garfield/projects/INSIGHT_analyses/valid_sites_file/validSites.txt
neutralFolder=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/50kb_withOverlaps_final/neutralSites
blockLocations=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/50kb_withOverlaps_final/blockLocations.bed
blockData=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/50kb_withOverlaps_final/blockData.pickle

#for freqCutOff in 5 15 25 35 50
#for freqCutOff in $(seq 1 50)
for freqCutOff in {5..50..5}
do
    bsub -o /dev/null python $mainScript $inputBed $validSites $neutralFolder $blockLocations $blockData $resultsFolder $deletePostSites $deleteLog False $freqCutOff #$screeningFile
done

#bsub -o /dev/null python $mainScript $inputBed $validSites $neutralFolder $blockLocations $blockData $resultsFolder $deletePostSites $deleteLog True $freqCutOff #$screeningFile
