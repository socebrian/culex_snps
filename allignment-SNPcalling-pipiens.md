# Proyec: SNPS Culex - ECOGEN
# pipieline for fastp trimming, alignment, variant calling and filtering
## sonia cebrian camison scebrian27@gmail.com/sonia.cebrian@ebd.csic.es
### 10/09/2024 - 15/1/2025

### 1: fastp trimming (culex conda env on CESGA pipiens 123) (talapas conda env fastp pipiens4 y perex)
1. fastqc.sh: runs fastqc in all files and then do multiqc over the results. I do it before and after applying fastp to:
    a. check the filters i need to apply on fastp
    b. be able to compare if fastp is cleaning what i need to
2.  list-array-fasp.sh: it creates a list ready to run fastp as a job array on many files within a directory
3. fastp-array-sh: runs fastp as a job array using all the files contained in [1]

```bash
mkdir fastqc
#first run of fastqc
sbatch -n 1 --cpus-per-task=10 --mem-per-cpu=20 -t 12:00:00 fastqc.sh 4.enero

mkdir fastp
#Do fastp with polyG trimming and min 30pb length
sbatch list-array-fastp.sh dir/to/files name-output.txt
sbatch --array=1-20%4 fastp-array.sh directory list-fastp-output.txt #directory -> 2.mayo or other folder on LUSTRE. Note that the scripts usually have pre-defined directories, only missing a folder
```
## 1.2- fastqc and multiqc to ensure everything looks right
```bash
fastqc $LUSTRE/3.repetidas/data/*.fastp.fq.gz -o $LUSTRE/3.repetidas/fastqc
```
Programs:
- multiqc v1.26
- fastqc v0.12
- fastp v0.24

## Error checking through the workflow
During all the steps, a lot of output and error files are generated in each job and it is impossible to check all of them 1 by 1 to ensure everything run smoothly. To check all of them at once, I created the script <error-check.sh>, that checks all the .err and .out files inside the current directory and prints the erros of all of them together. If the script is empty, it means it has not find the word "error" nor exit code not equal to 0 in any of these files. 
```bash
 sbatch -n 1 --mem-per-cpu=1G -t 00:30:00 --job-name=check-files error-check.sh . error-cut-mitoc.txt
```

### 2. alining (culex conda env on talapas)
Reference genome weas cut to retain only data on 3 main chromosomes
and I renamed the chromosomes to 1, 2 y 3
```bash
samtools faidx idCulPipi1.1.primary.fa SUPER_2 > idCulPipi.chrom.fa
samtools faidx idCulPipi1.1.primary.fa SUPER_3 >> idCulPipi.chrom.fa
samtools faidx idCulPipi1.1.primary.fa SUPER_1 >> idCulPipi.chrom.fa

sed 's/^>SUPER_2$/>2/' idCulPipi.chrom.fa > idCulPipi.chrom2.fa
sed 's/^>SUPER_3$/>3/' idCulPipi.chrom2.fa > idCulPipi.chrom3.fa
sed 's/^>SUPER_1$/>1/' idCulPipi.chrom3.fa > idCulPipi.chrom4.fa

#rm all except chrom4 and rename it to the original idCulPipi.chrom.fa

grep '^>' data/idCulPipi.chrom.fa
>2
>3
>1
```

Then was indexed using  
```bash
# indice para BWA
bwa index data/idCulPipi.chrom.fa
# indice general
samtools faidx data/idCulPipi.chrom.fa
# diccionario para gatk
samtools dict data/idCulPipi.chrom.fa -o data/idCulPipi.chrom.dict
```

Scripts to run alignment:
* Barches 1 to 3
    * list-arrays-bams.sh: Do List with names and lanes to distinguish different runs  ???? SEGURO
    * bwa-align-array.sh : get paired ends files and align. We are not currently using unpaired alligns for snp calling so they are not included iun the allignment.
    * BUSCAR COMO SE HIZO FLAGSTAT

* Batch 4: 
    * list-align-index.sh: creates a list of samples with all the info needed to run the alignment and indexing
    * align-bams-index.sh: runs the alignment and indexing of the bams. It also generates flagstat reports for each sample. The methods ure are the same as in the previous script, bu in this case this step was unified with sample indexing to save time.
    * flagstat-summary.sh: generates a summary report with the flagstat results of all the samples in a single file.
    
    
Then I examined which samples had less than 60% of properly mapped reads using:
 ><awk -F'\t' 'NR==1 || ($6 < 60)' summary.txt > below60-properly-mapped.txt>
This samples were not excluded yet, but the imformation was taking into acount.

```bash
mkdir 3.repetidas/bams
sbatch list-array-align.sh 3.repetidas/data list-align-output.txt 
sbatch --array=1-5 -n 1 --cpus-per-task=15 --mem-per-cpu=10G -t 12:00:00 --job-name=allign-3repes bwa-align-array.sh 3.repetidas list-align-output.txt idCulPipi.chrom.fa 
#repeat for the other batch of samples
```
##### bwa-align-array.sh
```bash
#variables
dir=${1} #3.repetidas
input=${2} #3.list-fastp-output.txt
ref=${3} # data/idCulPipi.chrom.fa

# Specify the path to the config file
config="/home/scamison/kernlab/scamison/${dir}/to-map/${input}"

# need to put in txt nombre=${1} # PP2114-c60.L7
name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)


#paired alignment
bwa mem \
    /home/scamison/kernlab/scamison/data/${ref}\
    /home/scamison/kernlab/scamison/${dir}/to-map/${name}.R1.fastp.fq.gz \
    /home/scamison/kernlab/scamison/${dir}/to-map/${name}.R2.fastp.fq.gz \
    -t 15|
samtools view -hbS - -o /home/scamison/kernlab/scamison/${dir}/bams/${name}.PR.bam -@ 15
```

### 3. BAM indexing (culex conda env on talapas)
For samples from batch 4, this step was inlcuded in precious alignment step. Batches 1 to 3 were done like this:

* <list-array-index.sh>: creates a list of taskid with samples an all their info of reads (reasds 1 or R2) and lanes to do the indexing
* <index-bam-array.sh>: take the list we just created and index all the bams with info about batch, lane, etc.

```bash
sbatch list-array-index.sh 1.noviembre/bams list-index-output.txt 
sbatch --array=1-15%5 -n 1 --cpus-per-task=1 --mem-per-cpu=2G -t 02:00:00 --job-name=3rep.index  index-bam-array.sh 3.repetidas list-indez-bams.txt
#repeat for the other batchs of samples
```

#### index-bam-array.sh
```bash
#variables
dir=${1} #3.repetidas
input=${2} #list-index-output.txt

# Specify the path to the config file
config="/home/scamison/kernlab/scamison/${dir}/bams/${input}"


IFS=$'\t' read -r -a samples <<< $(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$config")
# Extract parameters
id=${samples[1]}
lane=${samples[2]}
r=${samples[3]}
ext=${samples[4]}
#debug parameters not shown for simplification of the workflow
# File paths
input_bam="/home/scamison/kernlab/scamison/${dir}/bams/${id}.${lane}.${r}.bam"
sorted_bam="/home/scamison/kernlab/scamison/${dir}/bams/${id}.${lane}.${r}.sorted.bam"
sorted_rg_bam="/home/scamison/kernlab/scamison/${dir}/bams/${id}.${lane}.${r}.sorted.rg.bam"

# Sort BAM file
echo "Sorting BAM file: $input_bam"
samtools sort "$input_bam" -o "$sorted_bam"

# Add or replace read groups (picard line may differ depending on whether you can run picard directly or not, in talapas I have to activate java first)
echo "Adding or replacing read groups in BAM file: $sorted_bam"
picard AddOrReplaceReadGroups \ 
    I="$sorted_bam" \
    O="$sorted_rg_bam" \
    RGID="${id}.${lane}.${r}-seq" \
    RGLB="$dir" \
    RGPL=Illumina \
    RGPU="${id}.${lane}.${r}" \
    RGSM="$id" \
    VALIDATION_STRINGENCY=SILENT

# Index the sorted BAM file with read groups
echo "Indexing BAM file: $sorted_rg_bam"
samtools index "$sorted_rg_bam"
```
NOTA: AQUI SE METE TAMBIEN UN CHECKEO POR FLAGSTAT PARA VER CUANTAS SEC ESTAN PROPERLY MAPPED
NOTA PICARD VERSION USED FOR MITOC DE 4.ENERO conda install bioconda::picard


### 4. Mark duplicates (vcf environment on talapas)

We use 2 scripts for that:
* <mark-duplicates-array.sh>: it takes the list of samples works as an array eliminating all the duplicates from each file. As we are working with different runs (lanes) from the same samples that we are later going to merge, is better to erase all the duplicates instead of just marking them. 
Need to do <mkdir mark-duplicates> to run the script

* <list-mark-duplicates.sh> was used to create a list of samples to work with in batches 1 to 3. For batch 4 we modified mark-duplicates-array.sh to use the sample list generated at the beginning of the workflow. 

In the case of batch 4 of samples, we run flagstat again to check if the percentage of properly mapped reads was still low. The results were similar to the ones obtained before marking duplicates, so we decided to keep the samples with low mapping quality for now. In batches 1-3 we did not check this again.


As a result of using picard MarkdDuplicates we obtain for each sample on txt file with metrics about the duplications <PP1271.L6_marked_dup_metrics.txt>. To quickly check the percentage of duplications all the samples had I created one summary report for each bash of samples using this:

```bash
for file in /home/scamison/kernlab/scamison/1.noviembre/mark-duplicates/*_marked_dup_metrics.txt; do
    sample=$(basename "$file" | sed 's/_marked_dup_metrics.txt//')
    percent_dup=$(awk '/PERCENT_DUPLICATION/ {getline; print $9}' "$file")
    echo -e "${sample}\t${percent_dup}" >> 1.noviembre/mark-duplicates/duplicates-summary-1nov.txt
done
```
According tot the report, the duplication levels for each batch of samples were:
* 1.noviembre: 10-20% approx
* 2.mayo: 10-18% approx (mas centrado en valores intermedios que 1.noviembre)
* 3.repetidas: 31-35%
* 4.enero: 8-22% approx (mas variable, algunos llegando hasta 25%)

After this we indexed all samples
```bash
for bam in ./*bam ; do echo "indexing $bam"; samtools index $bam; done
```
#### mark-duplicates-array.sh
```bash
#variables
dir=${1} #3.repetidas
input=${2} #list of bam files

# Specify the path to the config file
sample_list="/home/scamison/kernlab/scamison/${dir}/bams/${input}"
bam_dir="/home/scamison/kernlab/scamison/${dir}/bams"
out_dir="/home/scamison/kernlab/scamison/${dir}/mark-duplicates"

# Read the sample ID for the current array task
sample=$(awk -v taskid="$SLURM_ARRAY_TASK_ID" '$1 == taskid {print $2}' "$sample_list")
in_bam="$bam_dir/${sample}.PR.sorted.rg.bam"
# Define output filenames
output_bam="$out_dir/${sample}.md.bam"
output_metrics="$out_dir/${sample}_marked_dup_metrics.txt"

# Debug: Print the sample variable
echo "Sample: $sample"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=$in_bam \
    O=$output_bam \
    M=$output_metrics \
    REMOVE_DUPLICATES=true
```

### 5. SORT AND Merge alignments form same sample, differente reads/lanes (culex conda env on talapas)
FIRST THING SORT THE MERGEDDDD BAMFILESSSSSSSS (NOT NEEDED FOR THAT)

We use 2 scripts:
* <list-to-merge-bams.sh> : generates a list of samples with only PP**** codes, so can be used in future steps that only need to identify samples. This list need to have info separated for sample id and lane.
*<merge-bams-array.sh>: merge all the indexed bams correspondign to the same sample (paired and unpaired reads)
```bash
mkdir 3.repetidas/merged
sbatch list-to-merge-bams.sh 3.repetidas/mark-duplicates list-to-merge-bams.txt 
sbatch --array=1-50%5 -n 1 --cpus-per-task=1 --mem-per-cpu=2G -t 02:00:00 --job-name=2mayo.merge merge-bams-array.sh 2.mayo list-merge-bams.txt
```
Then merged files should be indexing
```bash
for bam_file in 3.repetidas/merged/*.bam; do
    echo "Indexing BAM file: $bam_file"
    samtools index "$bam_file"
done
#repeat for the other batch
```
#### merge-bams-aray.sh
```bash
# Directory containing BAM files
dir=${1} #3.repetidas
input=${2} #3.list-to-merge-bams.txt

sample_list="/home/scamison/kernlab/scamison/${dir}/mark-duplicates/${input}"
bam_dir="/home/scamison/kernlab/scamison/${dir}/mark-duplicates"
out_dir="/home/scamison/kernlab/scamison/${dir}/merged"

# Read the sample ID for the current array task
sample=$(awk -v taskid="$SLURM_ARRAY_TASK_ID" '$1 == taskid {print $2}' "$sample_list")

# Debug: Print the sample variable
echo "Sample: $sample"

# Find all BAM files for the current sample
bam_files=$(ls "$bam_dir"/"${sample}"*.md.bam)

# Check if BAM files exist
if [ -z "$bam_files" ]; then
    echo "Error: Input BAM files not found for sample ${sample} at ${bam_dir}"
    exit 1
fi

# Merge the BAM files
output_bam="$out_dir/${sample}.merged.bam"
samtools merge "$output_bam" $bam_files
```
### 4.5 validation of merged files 
```bash
bash validate-merged.sh > validation.out 2> validation.err
```
We did not sortd bams again after merging, but we used  <picard ValidateSamFile I="$bam" MODE=SUMMARY> to do checks on the bam files an none od them had error regarding lack of sorting.

### 5.checking of final bams: quality, duplicates and coverage (culex env on talapas)
1. To check the quality of the mapping qualimaps + multibamqc 
2. To check if there is new PCR duplicates: re-run MarkDuplicates
3. control coverage along the genome (to-do) 

1. qualimaps + multibamqc - QualiMap v.2.3

First I did qualimaps of everything using <qualimap.sh> and then did <qualimap multi-bamqc> for each batch of samples.
To run <qualimap multi-bamqc> I need to create a configuration file with info of all the samples I want to add to the report (id, path to bamqc results and "condition"-in my case, batch-). 

```bash
find /home/scamison/kernlab/scamison/3.repetidas/qualimaps -mindepth 1 -maxdepth 1 -type d | awk -F/ '{print $NF "\t" $0 "\t3.repetidas"}' > multibamqc_config_file.txt
find /home/scamison/kernlab/scamison/2.mayo/qualimaps -mindepth 1 -maxdepth 1 -type d | awk -F/ '{print $NF "\t" $0 "\t2.mayo"}' >> multibamqc_config_file.txt
find /home/scamison/kernlab/scamison/1.noviembre/qualimaps -mindepth 1 -maxdepth 1 -type d | awk -F/ '{print $NF "\t" $0 "\t2.mayo"}' >> multibamqc_config_file.txt

#then run multibamqc
qualimap multi-bamqc -d multibamqc/multibamqc_config_file.txt -outdir "./multibamqc" -outfile "pipiens_1234batch_multibamqc_report.pdf" -outformat PDF:HTML
```

SUMMARY
Number of samples 	230
Number of groups 	4 (4 batches of sequencing)
Total number of mapped reads 	11,270,506,831
Mean samples coverage 	12.18
Mean samples GC-content 	39.17
Mean samples mapping quality 	32.85
Mean samples insert size 	316

PCA shows some batch effects clustering samples from 2-3 separated bt PC1 from samples from batch 4. Biggest differences are shown between samples within batch 4.
Batch 4 clearly shows bigger probblems with GC content and mapping quality across the reverence. Some samples ave better coverage but sime others have way less than in precvious batches.

#### qualimap.sh
```bash
dir=${1} #3.repetidas

# Paths and variables
input_dir="/home/scamison/kernlab/scamison/${dir}/merged"	
output_dir="/home/scamison/kernlab/scamison/${dir}/qualimaps"

## Loop through each merged BAM file in the input directory
for bam_file in "$input_dir"/*merged.bam; do
    # Check if the BAM file exists to avoid errors
    if [[ -f "$bam_file" ]]; then
        # Extract sample ID from filename
        sample_id=$(basename "$bam_file" .merged.bam)

        echo "Running qualimap for sample: ${sample_id}"
        qualimap bamqc -bam "$bam_file" -outdir "${output_dir}/${sample_id}" -outfile "${sample_id}.merged.bamqc.pdf"  -outformat PDF:HTML
    else
        echo "No BAM files found in ${input_dir}"
    fi
done
```
#### 
2. To check if there is new PCR duplicates: re-run MarkDuplicates
I marked duplicates again with <mark-duplicates2.sh> yo check if there was any new duplicates created from merging. I redid a report with the results an all samples from 3 batches had <0.5% of duplicates.

Batches 1,2,3 show that there are not duplicates created by the merging process. Batch 4 was assumed to have not created duplicate reads either, so it was not checked again. 


3. check the coverage along the genome
* <coverage-mosdepth.sh>: loop thorught all the files within the actual directory to calculate depth statistics on each of them
* <mosdepth-coverage.R>: Uses the output <sample.mosdepth.global.dist.txt> from each sample to plot the distribution of depth along the genome (proportion of genome with different levels of depth). Then uses <sample.regions.bed.gz> to plot the mean coverage on 100K bases windows along each of the chormosomes. This script is run directly on an R session on the shell.

HACER RESUMEN DE CADA MUESTRA 
```bash
sbatch -n 1 --cpus-per-task=10 --mem-per-cpu=20G -t 23:00:00 coverage-mosdepth.sh
```
* Overview of Mosdepth outputs
    * sample.mosdepth.global.dist.txt: A global coverage distribution across all bases in the genome “How many bases have coverage 0, coverage 1, coverage 2, etc.?”
    * sample.mosdepth.region.dist.txt: Coverage distribution restricted to the intervals in your BED file (if you used --by somefile.bed or --by <size>). Similar to the global distribution but only for those “regions” as defined by the user or by Mosdepth’s default intervals.
    * sample.mosdepth.summary.txt A one-line summary for each “region set” (or entire genome) with average depth, total bases covered, etc. Quick overall stats, not detailed coverage along the genome.
    * sample.per-base.bed.gz Coverage for each individual base position in your sequence. The file will have lines like: chrom    start    end    coverage where start/end differ by 1 base. This can be very large if your genome is big.  Useful if you want to do your own binning later, but heavy on memory/disk.
    * sample.regions.bed.gz  Coverage summarized by each region. Where do these regions come from?
        If you used --by 100000, then these are 100 kb windows.
        If you provided a BED file of intervals, these regions match those intervals.
        If you did not provide any --by arguments, Mosdepth will produce coverage by “contigs” or by smaller default windows (depending on your version/config).
    Each line typically has:  chrom    start    end    mean_coverage. This is usually the most convenient file for plotting coverage in fixed windows.
    * Index files (.csi): These are just index files for the .bed.gz coverage data. Typically you can ignore them unless your downstream tool (e.g., samtools, bcftools, or some R packages) needs random-access to the coverage file.




### 6.Variant calling (snpbcf env on talapas)
* <snp-calling-bcftools.sh>: This script gets a list of bam files <bam_list.txt> and a list with information about sample's origin population <population-mpileup.txt> and runs bcftools mpileup and call over all the samples calling variants only for the chromosome indicated when running sbatch (1, 2 or 3).


The list of files must indicate the complet path to each one of them so they can be reached even if they are located in different directories.
```bash
find /home/scamison/kernlab/scamison/3.repetidas/merged -mindepth 1 -maxdepth 1 -type f -name "P*.merged.bam" > SNPcalling/bam_list.txt
find /home/scamison/kernlab/scamison/2.mayo/merged -mindepth 1 -maxdepth 1 -type f -name "P*.merged.bam" >> SNPcalling/bam_list.txt
find /home/scamison/kernlab/scamison/1.noviembre/merged -mindepth 1 -maxdepth 1 -type f -name "P*.merged.bam" >> SNPcalling/bam_list.txt

```
#### snp-calling-bcftools.sh
```bash
thr=${SLURM_CPUS_PER_TASK}
chromosome=${1} #1, 2 or 3

# Paths and variables
dir="/home/scamison/kernlab/scamison/SNPcalling"
reference="/home/scamison/kernlab/scamison/data/idCulPipi.chrom.fa"
population_file="/home/scamison/kernlab/scamison/SNPcalling/populations-mpileup.txt"  # File with population information
 
# Chromosome BAM lists. list must contain the COMPLETE path to the bam files
list="${dir}/bam-list.txt" 	

# Run bcftools mpileup and call
bcftools mpileup -Ou  --threads $thr -f $reference -q 20 -Q 15 -b $list -r "$chromosome" --annotate "FORMAT/AD,FORMAT/DP,INFO/AD" | \
bcftools call -m  --threads $thr -G $population_file -Ob --format-fields "GQ,GP" -o ${dir}/pipiens123_${chromosome}_pop.bcf

# sbatch -n 1 --cpus-per-task=10 --mem-per-cpu=20G -t 12-00:00:00 --job-name=SNPcall2 snp-calling-bcftools.sh 2
```

#### How many sites are in each bcf?
Using <bcftools view -H pipiens123_3_pop.bcf | wc -l> 
<bcftools +counts vcf> is way faster
* chrom 1: 123,274,950 sites
* chrom 2: 206,781,587 sites
* chrom 3: 184,554,457 sites

### 7. Variant filtering (vcf env en poppy)
Se esta usando version 1.19
NOTA 2: TAMBIEN SI SE PUEDE FILTRAR LOS SNPS QUE QUEDEN A 5BP O MENOS DE INDELS 
NOTA 3: PLOTTEAR TODAS LAS VARIABLES A FILTRAR PARA DECIDIR TRESHOLDS (MIRAR PAPER POPPIP)

I do a frist filter of the bcf in which I filter out repetitions/TE that were detected by EarlGrey (FALTA PARTE DE COMO HICE LO DE EARL GRAY) 
```bash
bash filter-repeats.sh 1 4.enero > pipiens4_maskrep.out 2> pipiens4_maskrep.err
```
I did this for pipiens123, then for pipiens4, and the i merged the results 

```bash
bcftools merge pipiens123_1_norep.bcf pipiens4_1_norep.bcf -Ob --threads 30 -o pipiens1234_1_norep.bcf ; bcftools merge pipiens123_2_norep.bcf pipiens4_2_norep.bcf -Ob --threads 30 -o pipiens1234_2_norep.bcf ; bcftools merge pipiens123_3_norep.bcf pipiens4_3_norep.bcf -Ob --threads 30 -o pipiens1234_3_norep.bcf
```

Then I filtered and separated variants and invariants sites
```bash
bash filter-snps-invariants.sh "1" "pipiens1234_1" pipiens1234_1_norep.bcf > pip1234_1_invariants.out 2> pip1234_1_invariants.err
bash filter-snps-invariants.sh "2" "pipiens1234_2" pipiens1234_2_norep.bcf > pip1234_2_invariants.out 2> pip1234_2_invariants.err
bash filter-snps-invariants.sh "3" "pipiens1234_3" pipiens1234_3_norep.bcf > pip1234_3_invariants.out 2> pip1234_3_invariants.err
```
Pipiens1234_1_invariants.bcf:  16.488.511  variant sites
                              33.948.196 invariant sites
Pipiens1234_2_invariants.bcf:  34.868.347  variant sites
                             53.817.046 invariant sites
Pipiens1234_3_invariants.bcf:  23.700.404  variant sites    
                             48.657.554 invariant sites

This same filter generates reports about missing sites per individual, which I check
Is easy to spot without plotting the individuals that stand out, which match the ones that presented very low percentage of mapped reads. I filter those samples out

```bash
for bcf in *variant.bcf; do
echo "erasing samples with miss sites >60% $bcf"
name="${bcf%.bcf}"
output="${name}_clean.bcf"
echo "Output file: $output"
bcftools view -s ^PP1173,PP1615,PP1440,PP1589,PP1599,PP1867,PP1759,PP1803,PP1850,PP1621,PP1666 $bcf -Ob -o $output
echo "done with $bcf"
done 
```

I also generate some statistics about depth of coverage and plot the results using <plot-depth.R>
```bash
for bcf in *ldepth; do
name="${bcf%.ldepth}"
number=$(echo "$bcf" | cut -d'_' -f2)
    echo "Running Rscript for $bcf with chrom=$number and base name=$name"
Rscript depth_plot.R \
       "${number}" \
        "${name}.ldepth" \
        "${name}.ldepth.mean" \
       "${name}" 
done

#Rscript depth_plot.R 1 pipiens1234_1_variant_clean.ldepth pipiens1234_1_variant_clean.ldepth.mean pipiens1234_1
```

Filters I did apply after that:
* retain only biallellic SNPs
* exclude sirtes with INFO/DP> meanDP*2
* DP per site & individual >=5
* GQ per site & individual >=30
* Missinf GT per site < 50%
```bash
bash filter-variants.sh  pipiens1234_1_variant_clean.bcf   53023pipiens1234_1  5  30 |& tee 53023pipiens1234_1.log
```
This filter gives me the number of variants retained after each step (or the porportion of hwat would enter now the filter of missingness in case of DP and GQ, which don't filter out sites, but mark GT as ./.). It also gives metrics on missing site per individuals after each step. 

### hoe many sites are left after filtering?
dp>5 gq>30 miss<20 maf>=3%
chromosoma 1:
* variant sites 16488511
* biallelic snps 12722393
* variant sites miss<20% 10190233
* variant sites dp 12602478 (sin quitar missingness)
* sites g.depth >5 7.949.907
* sites gq>30 7170688 (al filtrar del todo missingness queda esto)
(MISS 0.5: 9477819)
* sites maf>=3%  1.025.647

chromosme 2:
* variant sites 34868347
* biallelic snps 26666156
* variant sites miss<20%  20810742
* variant sites dp 26537901 (sin quitar missingness)
* sites g.depth >5  16447977
* sites gq>30  15154646 (al filtrar del todo missingness queda esto)
(MISS 0.5:  19500000)
* sites maf>=3%  1.846.719

chromosoma 3:
* variant sites 23700404
* biallelic snps  19009093
* variant sites miss<20%  14994998
* variant sites dp 18870451 (sin quitar missingness)
* sites g.depth >5  11651718
* sites gq>30  10455113 (al filtrar del todo missingness queda esto)
* sites maf>=3%  1.565.303

The numbers of DP and GQ are counted as how many would be excluded If i applied a 20% missing site filter.
The MAF remaininf sites are calculated because I originally filtered that filter, but I used an intermediate file obtained after GQ filter for the next steps, as I identified some individuals that should be removed from the amalyses before further filtering.

After the 50% miss gt filter, I plot the missingness per individual using <plot-missingness.R> to see how many individuals have more than 50% of missing genotypes. FOr that, I first sorted, indexed and merged the 3 bcf files.

```bash
for bcf in *BDQM5.bcf; do
echo "sorting $bcf"
name="${bcf%.bcf}"
output="${name}_sorted.bcf"
echo "Output file: $output"
bcftools sort --temp-dir /sietch_colab/scamison/talapas/variant_filtering/temp-bcftools -Ob -o $output $bcf
echo "done with $bcf"
done
 #also indexed them (show just 1 as example)
bxftools index 530pipiens1234_1_BDQM5_sorted.bcf

bcftools concat -Ob --threads 30 -o 530pipiens123_all.bcf 530pipiens1234_1_BDQM5_sorted.bcf 530pipiens1234_2_BDQM5_sorted.bcf 530pipiens1234_3_BDQM5_sorted.bcf 
#index concat
bcftools index 530pipiens123_all.bcf
#the recalculate the missingness per individual
vcftools --bcf  530pipiens123_all.bcf --out 530pipiens123_all --missing-indv
bcftools +counts 
```
I plotted missing using plot_missingness.R in Rstudio.
Given the results, I decided to remove two more individuals that always show as outliers in missing sites.
```bash
bcftools view -s ^PP1271,PP1255 530pipiens123_all.bcf -Ob -o 530pipiens1234_all_clean.bcf
```
Then I applied the last filters:
* retain snps with <20% missing genotypes
* retain alleles with MAF >= 0.03
```bash
bash filter_variants_MISyMAF.sh 530pipiens1234_all_clean.bcf 53023pipiens1234_all_clean |& tee 53023pipiens1234_all_clean_BDQMF.log
#count final snps
bcftools view -H 53023pipiens1234_all_clean_BDQM.bcf  | cut -f1 | sort | uniq -c
```

After filter misiing <20% y maf >3% y depth >5 y gq>30
chrom1 1.049.728 
chrom2 1.884.757
chrom3 1.600.874 

```bash



