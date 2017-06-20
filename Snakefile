configfilename = 'atac_evo_config.yaml'
configfile: configfilename

workdir: config['snk_workdir']

# use a custom R profile for this project
shell.prefix()

rule all:
	input:
		config['snk_dhs_enh_reportdir'],
		config['snk_dhs_enh_folder']+"/all_dhs_enhancerstp68h.gff3",
		config["snk_transcripts"],
		config["snk_hsap_dnase_filt"],
		config["snk_fly_hum_directionality_plot"],

# ###############################################################################
# #########Tetsting snakmake
# ###############################################################################
rule testsnakeRscriptpre:
	output:
		config['snk_testinput']
	shell:
		"echo 'foobar' > {config[snk_testinput]}"

# this is just a test for using snakemake....
rule testsnakeRscript:
	input:
		config['snk_testinput']
	output:
		config['snk_testoutput']
	shell:
		"Rscript src/dermot/testsnakeRscript.R {output}"


# ###############################################################################
# #########Genome wide enhancer search.
# ###############################################################################

# this is the script to define enhancers based on the DHS data genome wide
rule define_enhancers_genomewide:
	input:
		config['snk_dnasepeakfolder'],
		config['snk_crm8008file'],
		config['snk_procapfolder'],
		config['snk_genesfile'],
		config['snk_cad4tissuerawboolfile'],
		config['snk_allpeaksfile'],
		config['snk_stageinffile'],
		"src/dermot/define_enhancers_genomewide.R"
	output:
		config['snk_dhs_enh_folder'],
		config['snk_dhs_enh_folder']+"/all_dhs_enhancerstp24h.gff3",
		config['snk_dhs_enh_folder']+"/all_dhs_enhancerstp68h.gff3"
	shell:
		"Rscript src/dermot/define_enhancers_genomewide.R"
# probs need a rule to create the reports folder....

# This rule creates a report on the dhs enhancers
rule define_enhancers_genomewide_report:
	input:
		'src/dermot/define_enhancers_report.R',
		{config['snk_dhs_enh_folder']},
		config['snk_dhs_enh_folder']+"/all_dhs_enhancerstp24h.gff3",
		config['snk_dhs_enh_folder']+"/all_dhs_enhancerstp68h.gff3"
	output:
		config['snk_dhs_enh_reportdir']
	shell:
		"mkdir "+config['snk_dhs_enh_reportdir']+" ; "
		"Rscript -e 'setwd(\"analysis/Reports\");knitr::spin(\"../../src/dermot/define_enhancers_report.R\")';"


# -------------------------------------------------------------------------------
# --------S2 analysis
# -------------------------------------------------------------------------------
# This rule loads the transcript and gene data for the s2 analysis
rule load_gene_data:
	input:	"src/dermot/load_genes_pol.R"
	output:
		config["snk_transcripts"],
		config["snk_human_transcripts"],
		config["snk_gene_models"],
	shell:
		"Rscript src/dermot/load_genes_pol.R"

rule assemble_human_dmel_s2_enhs:
	input:
		"src/dermot/assemble_human_dmel_s2_enhs.R",
		config["snk_transcripts"],
		config["snk_human_transcripts"],
		config["snk_s2_comp_enhancers"],
		config["snk_s2_dhs_folder"],
		config["snk_human_enh"],
		config["snk_hsap_dir"],
	output:
		config["snk_s2_tss"],
		config["snk_s2_dnase_filt"],
		config["snk_hg19_tss"],
		config["snk_hsap_dnase_filt"]
	shell:
		"Rscript src/dermot/assemble_human_dmel_s2_enhs.R"


# This rule uses a bash script to do liftovers and read counting with the
# Human and Fly
rule hsap_readcount:
	input:
		"src/dermot/hsap_readcount.sh",
		config["snk_hg19_chrom_sizes"],
		config["snk_hsap_dnase_filt"],
		config["snk_hg19_tss"],
		config["snk_hg18_19_lift"],
		config["snk_human_procap1"],
		config["snk_human_procap2"],
	output:
		config["snk_hsap_dnase_filt"]+".procap1.count",
		config["snk_hsap_dnase_filt"]+".procap2.count",
		config["snk_hg19_tss"]+".procap2.count",
		config["snk_hg19_tss"]+".procap2.count",
	shell:
		"bash src/dermot/hsap_readcount.sh config.json"
# now for fly...
rule dmel_readcount:
	input:
		"src/dermot/dmel_readcount.sh",
		config["snk_s2_tss"],
		config["snk_s2_dnase_filt"],
		config["snk_dmel_chr_sizes"],
		config["snk_dmel_groseq_plus"],
		config["snk_dmel_groseq_minus"],
	output:
		config["snk_s2_tss"]+".groseq.plus.bed",
		config["snk_s2_tss"]+".groseq.min.bed",
		config["snk_s2_dnase_filt"]+".groseq.plus.bed",
		config["snk_s2_dnase_filt"]+".groseq.min.bed",
	shell:
		"bash src/dermot/dmel_readcount.sh config.json"

# Finally produce the plots of readcounts
rule comp_readcounts:
	input:
		"src/dermot/comp_readcounts.R",
		config["snk_hsap_dnase_filt"]+".procap1.count",
		config["snk_hsap_dnase_filt"]+".procap2.count",
		config["snk_hg19_tss"]+".procap2.count",
		config["snk_hg19_tss"]+".procap2.count",
		config["snk_s2_tss"]+".groseq.plus.bed",
		config["snk_s2_tss"]+".groseq.min.bed",
		config["snk_s2_dnase_filt"]+".groseq.plus.bed",
		config["snk_s2_dnase_filt"]+".groseq.min.bed",
	output:
		config["snk_fly_hum_directionality_plot"],
		config["snk_readcount_comp_plot"],
	shell:
		"Rscript src/dermot/comp_readcounts.R"

# Now look at STARR-seq / DNase / GROseq relationshiop
rule starr_seq_groseq:
	input:
		"src/dermot/starr_seq_groseq.R",
		config["snk_s2_dhs_folder"],
		config["snk_s2_starr_folder"]+"*",
		config["snk_s2_starr_folder"],
		config["snk_s2_scRNA_folder"],
		config["snk_s2_cage_folder"],
		config["snk_s2_procap_folder"],
		config["snk_s2_dhs_cov_folder"],
	output:
		config["snk_readcount_comp_plot"],
		config["snk_s2_dhs_withexpr"],
		config["snk_dnase_eRNA_rocs_plot"],
	shell:
		"Rscript src/dermot/starr_seq_groseq.R"

# -------------------------------------------------------------------------------
# --------Calling CAGE clusters as per Phantom
# -------------------------------------------------------------------------------


# generate the ctss files for dpi, run it
rule run_dpi:
	input:
		"src/dermot/run_dpi",
		config["snk_s2_dhs_folder"],
		config["snk_s2_starr_folder"]+"*",
		config["snk_s2_starr_folder"],
		config["snk_s2_scRNA_folder"],
		config["snk_s2_cage_folder"],
		config["snk_s2_procap_folder"],
		config["snk_s2_dhs_cov_folder"],
	output:
		config["snk_readcount_comp_plot"],
		config["snk_s2_dhs_withexpr"],
		config["snk_dnase_eRNA_rocs_plot"],
	shell:
		"Rscript src/dermot/run_dpi.R"


		
rule run_dpi:
	input:
		"src/dermot/run_dpi",
		config["snk_s2_dhs_folder"],
		config["snk_s2_starr_folder"]+"*",
		config["snk_s2_starr_folder"],
		config["snk_s2_scRNA_folder"],
		config["snk_s2_cage_folder"],
		config["snk_s2_procap_folder"],
		config["snk_s2_dhs_cov_folder"],
	output:
		config["snk_readcount_comp_plot"],
		config["snk_s2_dhs_withexpr"],
		config["snk_dnase_eRNA_rocs_plot"],
	shell:
		"Rscript src/dermot/run_dpi.R"

# load up the dpi peaks, create some plots about them, filter them
rule process_dpi_peaks:
	input:
		"src/dermot/process_dpi_peaks.R",
		config["snk_s2_dhs_folder"],
		config["snk_s2_starr_folder"]+"*",
		config["snk_s2_starr_folder"],
		config["snk_s2_scRNA_folder"],
		config["snk_s2_cage_folder"],
		config["snk_s2_procap_folder"],
		config["snk_s2_dhs_cov_folder"],
	output:
		config["snk_readcount_comp_plot"],
		config["snk_s2_dhs_withexpr"],
		config["snk_dnase_eRNA_rocs_plot"],
	shell:
		"Rscript src/dermot/process_dpi_peaks.R"

# Use dpi peaks to predict enhancers
rule dpi_predict:
	input:
		"src/dermot/starr_seq_groseq.R",
		config["snk_s2_dhs_folder"],
		config["snk_s2_starr_folder"]+"*",
		config["snk_s2_starr_folder"],
		config["snk_s2_scRNA_folder"],
		config["snk_s2_cage_folder"],
		config["snk_s2_procap_folder"],
		config["snk_s2_dhs_cov_folder"],
	output:
		config["snk_readcount_comp_plot"],
		config["snk_s2_dhs_withexpr"],
		config["snk_dnase_eRNA_rocs_plot"],
	shell:
		"Rscript src/dermot/dpi_predict.R"


