# Tony continue Arthurs pipeline for mapping ribo-profiling reads
#### Run R script to intersect reads (sam file) and predicted.ORF (repre.valid.ORF.gtf file), then generate a table with read count on codon frame 1,2,3
https://github.com/Arthurhe/Riboprofiling_pipeline_Arthur


### Set up enivronment
```
cd PATH/TO/THE/REPOSITORY
conda env create -f tony_r_env.yml   #create the appropriate environment to run the Rscript

```

# Activate the environment we created
```
conda activate tony_r

```

   
# To run intersect with merged ORF (ORF with the same set of CDSs) 
```
Rscript /PATH/TO/Intersect_mergeORF.R -f ${Read.sam} -g ${ORF.gtf} -o ${Output}
 
```

${Read.sam} = the path to your read.sam file  
${ORF.gtf} = the path to your repre.valid.gtf file. If run on multiple samples, can use a concatenated gtf from multiple sample. The script will remove duplicates. <br>
${Output} = the path to your output directory and filename 

  
# for example:  
```
Rscript /home/tonycheng/Tools/Ribo_Tony/Intersect_mergeORF.R -f /home/tonycheng/Data/Fleur_data/sam/GSC2_t1.sam -g /home/tonycheng/Data/Fleur_data/gtf/GSC_total_uniq.repre.valid.ORF.gtf -o /home/tonycheng/Data/Fleur_data/GSC2_t1_mergeORF.tsv

```
    

