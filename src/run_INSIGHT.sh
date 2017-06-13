###A simple local script for running INSIGHT on TFBS data



screeningFile='/g/furlong/garfield/projects/INSIGHT_analyses/useful_annotations/ncRNA_and_exonic_masks/dmel-all-r5.57_genes_exons.bed'

inputBed=$1 #best to specify a full path
resultsFolder=$2 #best to specify a full path
deletePostSites=True #set to True at the moment...we almost never use these posterior sites
deleteLog=False ##set to False...the log file is the easiest place to read the results.
mainScript=/g/furlong/garfield/projects/INSIGHT_analyses/src/run_INSIGHT/run_INSIGHT.py
#freqCutOff=15
freqCutOff=$3

# #use here for standard version - 5kb blocks of non-coding sequence
# validSites=/g/furlong/garfield/projects/INSIGHT_analyses/valid_sites_file/validSites.txt
# neutralFolder=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/5kb_withOverlaps_updated/neutralSites
# blockLocations=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/5kb_withOverlaps_updated/blockLocations.bed
# blockData=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/5kb_withOverlaps_updated/blockData.pickle

#we can also code an alternative version for 4d sites
validSites=/g/furlong/garfield/projects/INSIGHT_analyses/valid_sites_file/validSites.txt
neutralFolder=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/50kb_withOverlaps_final/neutralSites
blockLocations=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/50kb_withOverlaps_final/blockLocations.bed
blockData=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/50kb_withOverlaps_final/blockData.pickle

# #or for BIG files
# validSites=/g/furlong/garfield/projects/INSIGHT_analyses/valid_sites_file/validSites.txt
# neutralFolder=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/500kb_withOverlaps_4d/neutralSites
# blockLocations=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/500kb_withOverlaps_4d/blockLocations.bed
# blockData=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/500kb_withOverlaps_4d/blockData.pickle

# #we can also run using the alternative method for ancestral state reconstruction
# validSites=/g/furlong/garfield/projects/INSIGHT_analyses/valid_sites_file/validSites_ucscAnc.txt
# neutralFolder=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/50kb_withOverlaps_final_UCSC_ancState/neutralSites
# blockLocations=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/50kb_withOverlaps_final_UCSC_ancState/blockLocations.bed
# blockData=/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/50kb_withOverlaps_final_UCSC_ancState/blockData.pickle

python $mainScript $inputBed $validSites $neutralFolder $blockLocations $blockData $resultsFolder $deletePostSites $deleteLog True $freqCutOff #$screeningFile
