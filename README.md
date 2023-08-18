# FLAQ-AMR (Florida Assembly Quality - AMR Detection)
FL BPHL's standard bacterial assembly pipeline for taxonomic characterization and AMR detection

## About
FLAQ-AMR was developed to analyze Illumina paired-end, whole-genome sequencing data generated from cultured bacterial isolates. The pipeline generates annotated, de novo assemblies along with reports including quality metrics, taxonomic classification, ST, and AMR/virulence detection. The current version will run only on [HiPerGator] (https://www.rc.ufl.edu/about/hipergator/) (HPG) using local Singularity containers for each pipeline process.

Stay tuned for FLAQ-AMR's upgrade to [Sanibel](https://github.com/BPHL-Molecular/Sanibel), a platform agnostic [Nextflow](https://www.nextflow.io/) workflow with added tools for expanded species-specific subtyping. Sanibel is currently under active development.

## Dependencies
- Python3
- Singularity/Apptainer
- Git

To load python3 into your current environment on HiPerGator, either use `module load python` to get the lastest version of python or activate your base conda environment. For more information on how to set up your base conda environment on HPG, see the [HiPerGator Analysis Reference Guide](https://github.com/StaPH-B/southeast-region/tree/master/hipergator)).

Singularity/Apptainer will be loaded as a module during your job execution on HPG using the sbatch job script in this repository. 

Git is already installed in your HPG environment upon login.

## Usage

For first time use, clone this repository to a directory in blue on HPG, such as in /blue/bphl-\<state\>/\<user\>/repos/bphl-molecular/.
```
cd /blue/bphl-<state>/<user>/repos/bphl-molecular/
git clone https://github.com/BPHL-Molecular/flaq_amr.git
```
For future use, update any changes to your local repository on HPG by navigating to the flaq_amr repository and pulling any changes.
```
cd flaq_amr/
git pull
```
To run the FLAQ-AMR pipeline, copy all files from the flaq_amr local repository to your analysis folder. Make an input directory and copy your fastqs.
```
mkdir <analysis_dir>
cd <analysis_dir>
cp /blue/bphl-<state>/<user>/repos/bphl-molecular/flaq_amr/* .
mkdir fastqs/
cp /path/to/fastqs/*.fastq.gz fastqs/
```
Rename your fastq files to the following format: sample_1.fastq.gz, sample_2.fastq.gz. See below for a helpful command to rename your R1 and R2 files.
```
cd fastqs/
for i in *_R1_001.fastq.gz; do mv -- "$i" "${i%-F*}_1.fastq.gz"; done
for i in *_R2_001.fastq.gz; do mv -- "$i" "${i%-F*}_2.fastq.gz"; done
```
Edit your sbatch job submission script to include your email to receive an email notification upon job END or FAIL. Replace ENTER EMAIL in `#SBATCH --mail-user=ENTER EMAIL` with your email address. Make sure there is no space between = and your email address. Edit additional sbatch parameters as needed to run your job succesfully, such as the length of time the job will run.

Submit your job.
```
sbatch sbatch_flaq_amr.sh
```

## Main processes
- [Fastqc](https://github.com/s-andrews/FastQC)
- [Trimmomatic](https://github.com/usadellab/Trimmomatic)
- [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
- [Mash](https://github.com/marbl/Mash)
- [MLST](https://github.com/tseemann/mlst)
- [Unicycler](https://github.com/rrwick/Unicycler)
- [Quast](https://github.com/ablab/quast)
- [Prokka](https://github.com/tseemann/prokka)
- [Kraken2](https://github.com/DerrickWood/kraken2)
- [AMRFinderPlus](https://github.com/ncbi/amr)

## Primary outputs

Outputs from each process for each individual sample can be found in a sample-specific subdirectory within the FLAQ-AMR analysis directory. Report.txt contains the main summary report with quality metrics, taxonomic id, and ST. Additional details can be found in the report outputs from each process. Final assemblies (.fasta), annotated features (.gff), and amrfinder results are copied into the run directory for easier access for use in downstream analyses.

```
analysis_dir/
|__ <date>_flaq_run/
     |__ report.txt
     |__ sample1/
     |__ sample2/
|__ amrfinder_results/
|__ annotations/
|__ assemblies/
```

## Developed by:
[@SESchmedes](https://www.github.com/SESchmedes)<br />

Please email bphl16BioInformatics@flhealth.gov for questions and support.
