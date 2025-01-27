# Proyec: SNPS Culex - ECOGEN
# pipieline for fastp trimming, alignment, variant calling and filtering
## sonia cebrian camison scebrian27@gmail.com/sonia.cebrian@ebd.csic.es
### 10/09/2024 - 15/1/2025

### 1: fastp trimming (culex conda env on CESGA)
1. list-array-fasp.sh: it creates a list ready to run fastp as a job array on many files within a directory
2. fastp-array-sh: runs fastp as a job array using all the files contained in [1]

```bash
#Do fastp with polyG trimming and min 30pb length
sbatch list-array-fastp.sh dir/to/files name-output.txt
sbatch --array=1-20%4 fastp-array.sh directory list-fastp-output.txt #directory -> 2.mayo or other folder on LUSTRE. Note that the scripts usually have pre-defined directories, only missing a folder
```
## 1.2- fastqc and multiqc to ensure everything looks right
```bash
fastqc $LUSTRE/3.repetidas/data/*.fastp.fq.gz -o $LUSTRE/3.repetidas/fastqc

```
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
* list-arrays-bams.sh: Do List with names and lanes to distinguish different runs
* bwa-align-array.sh : get paired ends files and align. Then align unpaired reads files. Generates three outputs. (R1 and R2 unpaired, and paired reads)
As they are too much samples to run all the array at once, we do it in groups of 50 tasks (when doing in CESGA, not necessary on talapas), so we run the script several times for eacg batch of samples.

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

We use 2 scripts for that:
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

### 4. Mark duplicates (vcf environment on talapas)
We use 2 scripts for that:
* <list-mark-duplicates.sh>: it creates a list with taskid and samples names (i.e. PP1255.L2) to be used as a config file in the next script
* <mark-duplicates-array.sh>: it takes the list of samples created and works as an array eliminating all the duplicates from each file. As we are working with different runs (lanes) from the same samples that we are later going to merge, is better to erase all the duplicates instead of just marking them. 

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

### 5. Merge alignments form same sample, differente reads/lanes (culex conda env on talapas)

We use 2 scripts:
* <list-to-merge-bams.sh> : generates a list of samples with only PP**** codes, so can be used in future steps that only need to identify samples
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


### 5.checking of final bams: quality, duplicates and coverage (culex env on talapas)
1. To check the quality of the mapping qualimaps + multibamqc 
2. To check if there is new PCR duplicates: re-run MarkDuplicates
3. control coverage along the genome (to-do) 

1. qualimaps + multibamqc 
First I did qualimaps of everything using <qualimap.sh> and then did <qualimap multi-bamqc> for each batch of samples.
To run <qualimap multi-bamqc> I need to create a configuration file with info of all the samples I want to add to the report (id, path to bamqc results and "condition"-in my case, batch-). 

```bash
find /home/scamison/kernlab/scamison/3.repetidas/qualimaps -mindepth 1 -maxdepth 1 -type d | awk -F/ '{print $NF "\t" $0 "\t3.repetidas"}' > multibamqc_config_file.txt
find /home/scamison/kernlab/scamison/2.mayo/qualimaps -mindepth 1 -maxdepth 1 -type d | awk -F/ '{print $NF "\t" $0 "\t2.mayo"}' >> multibamqc_config_file.txt
find /home/scamison/kernlab/scamison/1.noviembre/qualimaps -mindepth 1 -maxdepth 1 -type d | awk -F/ '{print $NF "\t" $0 "\t2.mayo"}' >> multibamqc_config_file.txt

#then run multibamqc
qualimap multi-bamqc -d multibamqc/multibamqc_config_file.txt -outdir "./multibamqc" -outfile "123batch_multibamqc_report.pdf" -outformat PDF:HTML
```
When analysing 3 first batches (1.noviembre, 2.mayo, 3.repetidas), the PCA plot in the multi-bamqc report shows the variance in the data across the first two principal components (PC1 and PC2). The samples corresponding to different batches (3.repetidas, 2.mayo, and 1.noviembre) form distinct clusters on the PCA1 vs PC2 graph.
- The clustering of samples by batch indicates the presence of batch effects.
- Batch effects are non-biological variations introduced during sample preparation, sequencing, or other technical processes.
It is neccessary to use batch effects corrections in downstream analyses. 

SUMMARY
Number of samples 	96
Number of groups 	3
Total number of mapped reads 	3,449,399,688
Mean samples coverage 	9.15
Mean samples GC-content 	39.33
Mean samples mapping quality 	33.49
Mean samples insert size 	354

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

3. check the coverage along the genome
FALTA COVERAGE CHECK

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
* chrom 1: 123,274,950 sites
* chrom 2: 206,781,587 sites
* chrom 3: 184,554,457 sites

### 7. Variant filtering (vcf env en poppy)

NOTA: UPDATE VERSION DE BCFTOOLS A 1.21 ANTES DE REHACER ESTE PASO PARA FUTURAS MUESTARS

To do a first check of the vcf I use the python notebook <explore_vcf_examples.ipynb>
I do a frist filter of the bcf in which I filter out repetitions/TE that were detected by EarlGrey (FALTA PARTE DE COMO HICE LO DE EARL GRAY) and retain only sites that have a minimum of 1 read for one of the minor alleles.
To do so I use <variant-filter1n.sh> (it removes repetitions and filter by  --min-ac 1).

#### How many CALLABLE sites are in each chromosome?
Using <bcftools view -H pipiens123_1_filter1n_qd.vcf.gz | wc -l>
* chrom 1: 12,673,699
* chrom 2: 23,742,523
* chrom 3: 20,459,567


#### variantfilter1n.sh
```bash
chromosome=${1}
bcf=${2} #pipiens123

# Define file paths based on the chromosome number
repeat_bed="variant_filtering/repeats/s${chromosome}.cxpipiens.repeats.bed"
input_bcf="SNPcalling/${bcf}_${chromosome}_pop.bcf"
repeat_masked_vcf="variant_filtering/${bcf}_${chromosome}_repeatmasked.vcf.gz"
filtered_vcf="variant_filtering/${bcf}_${chromosome}_filter1n.vcf.gz"

# Print the file paths
echo "Repeat BED file: ${repeat_bed}"
echo "Input BCF file: ${input_bcf}"
echo "Repeat masked VCF file: ${repeat_masked_vcf}"
echo "Filtered VCF file: ${filtered_vcf}"

# Run bcftools commands
bcftools filter -M ^$repeat_bed --soft-filter repeat_masked_super${chromosome} --threads 30 -Oz -o $repeat_masked_vcf $input_bcf
bcftools view -f PASS,. --min-ac 1 --threads 30 -Oz -o $filtered_vcf $repeat_masked_vcf
```
#### 
These vcf do not include information on quality by depth (QD), which is necessary for the filtering so I add it using <qd.sh>: <qd.sh bash qd.sh pipiens123_1_filter1n.vcf.gz pipiens123_1_filter1n_qd.vcf.gz> for each chromosome.
#### qd.sh
```bash
#!/usr/bin/env bash
zcat $1 | gawk '
BEGIN{OFS="\t"}
($1 ~ /^##/) { print }
($1 == "#CHROM") {
print "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Quality by depth\">"
print $0
}
!($1 ~ /^#/) {
QUAL=$6;
match($8, /DP=([0-9]+)/, DP);
if (DP[1] > 0) {
QD = QUAL/DP[1]
} else {
QD = 0.0
}
$8=$8";QD="QD
print
}
' | bgzip > $2
```
####

Then the vcf must be indexed using tabix
```bash
tabix -p vcf pipiens123_1_filter1n_qd.vcf.gz
```

Filters I did apply after that:

1. Using <variant-filter-BDQF.sh>
* retain only biallellic SNPs
* DP per site & individual >=5
* GQ pero site & individual >=20  (Using 30 I lost a few thousand more, here I am just trying to have something to star working with)
* MAC (minor allele count) >= 3   (If I did filter by MAF >=0.05 I ended up with )

```bash
sbatch -n 1 --cpus-per-task=15 --mem-per-cpu=20G -t 12:00:00 --job-name=DPQF_1 variant-filter-BDQF.sh pipiens123_1_filter1n_qd.vcf.gz pipiens123_1_f1n
```
#### variant-filter-BDQF.sh
```bash
input=${1} #
output=${2}

module load miniconda3
conda activate vcf

bcftools view -v snps -m2 -M2 "$input" | bcftools filter -e 'FORMAT/DP<5' --threads 30 -Ob -o "${output}_BD3.bcf"
bcftools filter -e 'FORMAT/GQ<20' --threads 30 -Ob -o "${output}_BDQ.bcf" "${output}_BD3.bcf"
#bcftools filter -e 'MAF<0.05' --threads 30 -Ob -o "${output}_BDQF.bcf" "${output}_BDQ.bcf"
bcftools filter -e 'MAC<3' --threads 30 -Ob -o "${output}_BDQF.bcf" "${output}_BDQ.bcf"
# sbatch -n 1 --cpus-per-task=15 --mem-per-cpu=20G -t 12:00:00 --job-name=DPQF_1 variant-filter-BDQF.sh pipiens123_1_filter1n_qd.vcf.gz pipiens123_f1n
```

### How many SNPs I have?
* after biallelic, DP and GQ filter
    * chrom 1: 92.586
    * chrom 2: 111.795
    * chrom 3: 128.944

* afer MAC filter:
    * chrom1: 13.116
    * chrom 2: 15.372
    * chrom 3: 17.709


2. Then I wanted to filter by QD>2 or even QD>1  using <variant-filtering-QD.sh>,  but three chromosomes retainedonly 20 to 50 SNPs.
```bash
 sbatch -n 1 --cpus-per-task=15 --mem-per-cpu=20G -t 12:00:00 --job-name=DPQF_1 variant-filter-QD.sh pipiens123_1_f1n_BDQF.bcf pipiens123_1_1nBDQF
```
#### variant-filter-QD.sh
```bash
input=${1} #
output=${2}
module load miniconda3
conda activate vcf

bcftools filter -e 'QD<1' -Ob -o "${output}_QD.bcf" "$input"
# sbatch -n 1 --cpus-per-task=15 --mem-per-cpu=20G -t 12:00:00 --job-name=DPQF_1 variant-filter-BDQF.sh pipiens123_1_filter1n_qd.vcf.gz pipiens123_f1n
```
####

3. Filter by Missing Data <20%
I did <bcftools filter -e 'F_MISSING<0.20' pipiens123_1_1nBDQF_QD.bcf | bcftools view -H | wc -l> and had no output, so none of the chromosomes have misisng data after BQDF filter.
To check I still did <vcftools --bcf pipiens123_1_f1n_BDQF.bcf --missing-indv --out pipiens123_1_f1n_BDQF> and the <.imiss> reports for three chromosomes reported 0 missing genotypes
I still run <variant-filter-M.sh> and got the same number of SNPs as before.
#### SNPS after filtering Miising genotypes (not QD filter applied)
    * chrom1: 13.116
    * chrom 2: 15.372
    * chrom 3: 17.709

That are the files that I will be working with for now. NOTE: 4 of my samples are males, and I did not took this into account when working with chormosome 1.
```bash
 vcftools --bcf pipiens123_1_f1n_BDQF.bcf --missing-indv --out pipiens123_1_f1n_BDQF
```
#### variant-filter-M.sh
```bash
input=${1} #
output=${2}

module load miniconda3
conda activate vcf

echo 'F_MISSING<0.2'

bcftools filter -e 'F_MISSING>=0.2'  -Ob -o "${output}_M.bcf" "$input"
```


