## Workflow pipeline for Gene-Module Heritability Analysis

Step1: 
       Start with a csv or table format file, comprising of a data frame of size G x M where G is the number of genes and M is the number of modules.
       For example: start with the file: '/Users/kushaldey/Documents/singlecellLDSC/data/healthyvinflamed_score.csv'
       Then run the follwoing script: 
       **Code to run**: *workflow_process_modules.R*
       This will save M number of files saved by the names of the cell-types or module names that are column names of the file `healthyvinflamed_score.csv`
       to the desired output gene module directory: `/Users/kushaldey/Documents/singlecellLDSC/data/HealthyVInflamed_gene_score_Feb1` 
       We may also save all genes in a separate folder: here `/Users/kushaldey/Documents/singlecellLDSC/data/All_IBD_genes.txt`

Step2:
	In the next step, we carry out processing of the M gene module files created in Step 1. For this step, we need a cluster, so we may want to move the
	files to the cluster. In my case, I move these M files to `/n/groups/price/kushal/singlecellLDSC/data/Gene_Modules/HealthyVInflamed_gene_score_Feb1`.
	Next we build annotations from these gene score files using 100kb, 5kb or Roadmap union across tissues strategies. Batch jobs to submit here for each
	of the M files are:
	**Code to run**: *build_module_annotations.sh*
	This will generate in a user-specified `bed_dir` M many folders with each folder comprising of three bedfiles (bedgraph) corresponding to
	100kb (100kb.bed), 5kb (5kb.bed) and Roadmap union (Roadmap_Enhancer.bed). However these bed files require some cleaning as they may possess many
	overlapping intervals and may not be sorted. This leads us to Step3. The bedfiles directory for our example is:
	'/n/groups/price/kushal/singlecellLDSC/data/BEDFILES/HealthyVInflamed_gene_score_Feb1'. 

Step3: 
       In this Step, we clean the bedfiles created in Step2 to sort them and for overlapping intervals, use the maximum annotation. 
       **Code to run**: *clean_bedgraphs.sh*
       Get a compute node and run this straight sequentially. Better that way than submitting multiple jobs to the cluster.   

Step4:
	In this step, we use the cleaned bedgraph files from Step3 to generate SNP-level annotations. In this step, we create inside a particular annotation
	cell, M annotation folders, each folder containing 3 subfolder named 100kb, 5kb and Roadmap_Enhancer and each folder will contain 22 annotation  files 
	(chr*.annot.gz) after the 22 chromosomes. The annotation folders in case of our example are in:
	'/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/HealthyVInflamed_gene_score_Feb1'.
	**Code to run**: *create_annotation_from_bedgraph.sh*, *make_annot_combine_from_bedgraph.py*
	The input to this will the bedfiles directory from Step3 and the output will be these newly created annotation files in the annotation folders.

Step5:
	After creating the annotations, we run ldsc script to compute the LD-scores for the different annotations in these annotation files.
	**Code to run**: *ldsc_mega.sh*
	This will create the (*chr*.ldscore.gz) as well as other additional files. The number of files to run will be M*22 as batch submission.

Step6:
	We run the S-LDSC regression model for T traits where T can be taken to be 64 (traits with above 10% heritability and of particular interest from
	medical point of view). For each of the M modules and 3 annotation types (100kb, 5kb, Roadmap_Enhancer), we run the regression model for T traits.
	**Code to run**: *ldsc_reg.sh*
	The number of files run is T*M*3. The output of this run will create a number of (*.results) files in the user specified 'results_dir'. The results 
	directory for our example is `/n/groups/price/kushal/singlecellLDSC/data/LDSC_RESULTS/HVI_gene_score_Feb1_plus_Healthy_plus_All`.
	For this step, the user also has to provide a baseline against which to compare. For example, we use the baseline called 
	'/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Baselines/baselineLD_MAF'.

Step7: 
       Compute the standard deviation of the annotations that will be required for the post-processing of the S-LDSC output. 
       **Code to run**: *compute_sd_annot.R*

Step 8:
      Postprocess the results files in the results directory and using the annotation standard deviation etc to compute matrix of tau-star, p(tau-star),
      Enrichment, p(Enrichment) and save that information for each module across all traits in the results directory under 
      '/n/groups/price/kushal/singlecellLDSC/data/LDSC_RESULTS/HVI_gene_score_Feb1_plus_Healthy_plus_All/baselineLD_MAF/$modulename'.
