Bulk alingment pipeline used in Itaconate study. I used `snakemake/snakemake:v7.30.0"` docker image to execute it. 

To repeat alignment you need fastqc files (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277689)  
1. Modify `config.yml` to adjust file locations.
2. Check that all snakemake files use proper file names as input.
3. Activate conda star_env  with `conda activate star_env`
4. Run pipeline with default snakemake run command in workflow directory `snakemake --cores  4`  (or change number of cores to available for you)


To run R analysis notebook specified in `scripts` directory you will need to create your own R environment and install all the packages mentioned in scripts ( I used R 4.3.0 ). 
