#!/usr/bin/env python

#Author: Sarah Schmedes
#Email: sarah.schmedes@flhealth.gov

'''
This program runs the FLaq pipeline in addition to AMR detection.
Assumes installation of python3.7+, Singularity.
Assumes paired-end data sample_[R]1[2].fastq[.gz]
'''

import os
import sys
import subprocess
import argparse
import datetime
import fnmatch
import re
import pandas as pd

#Parse arguments, get path for fastqs
parser = argparse.ArgumentParser(usage='flaq_amr.py <input_dir> [options]')
parser.add_argument('input', help='path to dir with fastqs')
parser.add_argument('--threads', default=8, dest='threads', help='specify number of threads, (default: %(default)s)')

if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

input_dir = args.input
threads = str(args.threads)

if input_dir[-1] != '/':
    input_dir = input_dir + '/'

out_path = os.getcwd()
output_dir = out_path + '/' + datetime.date.today().strftime('%Y-%m-%d') + '_flaq_run'
subprocess.run('mkdir -p ' + output_dir, shell=True, check=True) #make output directory date_flaq_run

#Get sample names
samples = []
fastqs = []

for f in os.listdir(input_dir):
    if fnmatch.fnmatch(f, '*.fastq*'):
        fastqs.append(f)
        sn = f.split("_")
        sn = sn[0]
        samples.append(sn)
unique = set(samples)
samples = list(unique)
samples.sort()

#Create output file
report = open(output_dir + '/report.txt', 'w')
header = ['sampleID', 'speciesID_mash', 'nearest_neighbor_mash', 'mash_distance', 'speciesID_kraken', 'kraken_percent', 'mlst_scheme', 'mlst_st', 'num_clean_reads', 'avg_readlength', 'avg_read_qual', 'est_coverage', 'num_contigs', 'longest_contig', 'N50', 'L50', 'total_length', 'gc_content', 'annotated_cds'] 
report.write('\t'.join(map(str,header)) + '\n')

#Make directories for output files
subprocess.run('mkdir -p assemblies', shell=True, check=True)
subprocess.run('mkdir -p annotations', shell=True, check=True)
subprocess.run('mkdir -p amrfinder_results', shell=True, check=True)

#Run analysis pipeline for each sample
for s in samples:
    sample_dir = output_dir + '/' + s + '/'
    subprocess.run('mkdir -p ' + sample_dir, shell=True, check=True) #mkdir for each sample name
    subprocess.run('cp ' + input_dir + s +'*.fastq* ' + sample_dir, shell=True, check=True) #cp fastqs to new dir
    
    out_log = open(sample_dir + s + '.out', 'w')
    err_log = open(sample_dir + s + '.err', 'w')
    #Run fastqc on original reads
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/fastqc_0.11.9.sif fastqc ' + sample_dir + '*.fastq* --threads ' + threads, shell=True, stdout=out_log, stderr=err_log, check=True)
    #Rename fastqc output files
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.html ' + sample_dir + s + '_1_original_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.zip ' + sample_dir + s + '_1_original_fastqc.zip', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.html ' + sample_dir + s + '_2_original_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.zip ' + sample_dir + s + '_2_original_fastqc.zip', shell=True, check=True)

    #Run trimmomatic
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/trimmomatic_0.39.sif trimmomatic PE -threads ' + threads + ' -phred33 -trimlog ' + sample_dir + s + '.log ' + sample_dir + s + '_*1.fastq* ' + sample_dir + s + '_*2.fastq* ' + sample_dir + s + '_1.trim.fq.gz ' + sample_dir + s + '_unpaired_1.trim.fq.gz ' + sample_dir + s + '_2.trim.fq.gz ' + sample_dir + s + '_unpaired_2.trim.fq.gz ILLUMINACLIP:/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:71 TRAILING:20', shell=True, stdout=out_log, stderr=err_log, check=True)
    #rm unpaired reads
    subprocess.run('rm ' + sample_dir + s + '_unpaired*.trim.fq.gz', shell=True, check=True)
    #rm fastq files copied from previous dir
    subprocess.run('rm ' + sample_dir + s + '*.fastq.gz', shell=True, check=True)

    #Run bbduk to remove Illumina adapter sequences and any PhiX contamination
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/bbtools_38.76.sif bbduk.sh in1=' + sample_dir + s + '_1.trim.fq.gz in2=' + sample_dir + s + '_2.trim.fq.gz out1=' + sample_dir + s + '_1.rmadpt.fq.gz out2=' + sample_dir + s + '_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=' + sample_dir + s + '.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/bbtools_38.76.sif bbduk.sh in1=' + sample_dir + s + '_1.rmadpt.fq.gz in2=' + sample_dir + s + '_2.rmadpt.fq.gz out1=' + sample_dir + s + '_1.fq.gz out2=' + sample_dir + s + '_2.fq.gz outm=' + sample_dir + s + '_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=' + sample_dir + s + '_phixstats.txt', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('rm ' + sample_dir + '*.trim.fq.gz', shell=True, check=True)
    subprocess.run('rm ' + sample_dir + '*.rmadpt.fq.gz', shell=True, check=True)

    #Run fastqc on trimmed forward and reverse reads
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/fastqc_0.11.9.sif fastqc ' + sample_dir + '*.fq.gz --threads ' + threads, shell=True, stdout=out_log, stderr=err_log, check=True)
    #Rename fastqc output files
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.html ' + sample_dir + s + '_1_clean_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_1_fastqc.zip ' + sample_dir + s + '_1_clean_fastqc.zip', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.html ' + sample_dir + s + '_2_clean_fastqc.html', shell=True, check=True)
    subprocess.run('mv ' + sample_dir + s + '_2_fastqc.zip ' + sample_dir + s + '_2_clean_fastqc.zip', shell=True, check=True)

    #Run multiqc
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/multiqc_1.8.sif multiqc ' + sample_dir + '*_fastqc.zip -o ' + sample_dir, shell=True, stdout=out_log, stderr=err_log, check=True)

    #mkdir for mash_output
    subprocess.run('mkdir -p ' + sample_dir + 'mash_output', shell=True, check=True)
    #Run mash for species ID
    subprocess.run('zcat ' + sample_dir + s + '_1.fq.gz ' + sample_dir + s + '_2.fq.gz | singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/mash_2.2.sif mash sketch -r -m 2 - -o ' + sample_dir + 'mash_output/' + s + '_sketch', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/mash_2.2.sif mash dist /db/RefSeqSketchesDefaults.msh ' + sample_dir + 'mash_output/' + s + '_sketch.msh > ' + sample_dir + 'mash_output/' + s + '_distances.tab', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('sort -gk3 ' + sample_dir + 'mash_output/' + s + '_distances.tab | head > ' + sample_dir + 'mash_output/' + s + '_distances_top10.tab', shell=True, stdout=out_log, stderr=err_log, check=True)
    #Get mash result (Genus and species ID)
    with open(sample_dir + 'mash_output/' + s + '_distances_top10.tab', 'r') as mash:
        top_hit = mash.readline()
        gn = re.sub(r'.*-\.-', '', top_hit)
        gn = gn.split()
        gn = gn[0]
        gn = re.split(r'^([^_]*_[^_]*)(_|\.).*$', gn)[1]
        genus = gn.split('_')[0]
        species = gn.split('_')[1]
        distance = top_hit.split()[2]
        accession = top_hit.split("-")[5]
    
    #Run unicycler
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/unicycler_0.4.7.sif unicycler -1 ' + sample_dir + s + '_1.fq.gz -2 ' + sample_dir + s + '_2.fq.gz -o ' + sample_dir + s + '_assembly --min_fasta_length 300 --keep 1 --min_kmer_frac 0.3 --max_kmer_frac 0.9 --verbosity 2', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('mv ' + sample_dir + s + '_assembly/assembly.fasta ' + sample_dir + s + '_assembly/' + s + '.fasta', shell=True, check=True)
    #Run Quast
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/quast_5.0.2.sif quast.py -o ' + sample_dir + s + '_assembly/quast_results/ ' + sample_dir + s +'_assembly/' + s + '.fasta', shell=True, stdout=out_log, stderr=err_log, check=True)
    df = pd.read_table(sample_dir + s + '_assembly/quast_results/report.tsv', sep="\t")
    assem = list(df.columns)[1]
    contigs = df[assem][12].astype(int)
    long_contig = df[assem][13].astype(int)
    n50 = df[assem][16].astype(int)
    l50 = df[assem][18].astype(int)
    genome = df[assem][14].astype(int)
    gc = df[assem][15].astype(int)
    #Run cg-pipeline scripts
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/lyveset_1.1.4f.sif run_assembly_shuffleReads.pl -gz ' + sample_dir + s + '_1.fq.gz ' + sample_dir + s + '_2.fq.gz > ' + sample_dir + s + '_clean_shuffled.fq.gz', shell=True, stdout=out_log, stderr=err_log, check=True)
    subprocess.run('singularity exec --bind ' + sample_dir + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/lyveset_1.1.4f.sif run_assembly_readMetrics.pl ' + sample_dir + s + '_clean_shuffled.fq.gz -e ' + genome.astype(str) + ' > ' + sample_dir + s + '_readMetrics.txt', shell=True, stdout=out_log, stderr=err_log, check=True)
    #Parse readMetrics
    with open(sample_dir + s + '_readMetrics.txt', 'r') as metrics:
        firstline = metrics.readline()
        firstline = firstline.rstrip()
        firstline = firstline.split()
        secondline = metrics.readline()
        secondline = secondline.rstrip()
        secondline = secondline.split()
        rm = dict(zip(firstline, secondline))
        avg_read_len = rm['avgReadLength']
        avg_qual = rm['avgQuality']
        num_read = rm['numReads']
        cov = rm['coverage']

    #Run prokka
    subprocess.run('singularity exec --bind ' + sample_dir + s + '_assembly:/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/prokka_1.14.5.sif prokka --cpus ' + threads + ' --genus ' + genus + ' --species ' + species + ' --strain ' + s + ' --outdir ' + sample_dir + s + '_assembly/prokka --prefix ' + s + ' --force --compliant --locustag ' + genus + ' ' + sample_dir + s + '_assembly/' + s + '.fasta', shell=True, stdout=out_log, stderr=err_log, check=True)
    #Parse prokka summary output
    cds = ''
    with open(sample_dir + s + '_assembly/prokka/' + s + '.txt', 'r') as genes:
        for line in genes:
            line = line.rstrip()
            content = line.split()
            if content[0] == 'CDS:':
                cds = content[1]

    #Run ncbi's amrfinderplus
    subprocess.run('singularity exec --cleanenv -B /tmp:/tmp -B ' + sample_dir + s + '_assembly/:/data /apps/staphb-toolkit/containers/ncbi-amrfinderplus_3.10.1.sif amrfinder -n /data/' + s + '.fasta --plus -o /data/' + s + '_amrfinderplus_report.tsv', shell=True, stdout=out_log, stderr=err_log, check=True)

    #Run mlst for all samples
    subprocess.run('singularity exec --bind ' + sample_dir + s + '_assembly:/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/mlst_2.19.0.sif mlst ' + sample_dir + s + '_assembly/' + s + '.fasta --nopath > ' + sample_dir + s + '_assembly/' + s + '.mlst', shell=True, stdout=out_log, stderr=err_log, check=True)
    #Parse mlst output
    with open(sample_dir + s + '_assembly/' + s +  '.mlst') as mlst:
        for line in mlst:
            out = line.rstrip().split()
            scheme = out[1]
            st = out[2]

    #Run Kraken2 with the mini db to confirm top species hit
    kraken_out = sample_dir + 'kraken_out/'
    subprocess.run('mkdir -p ' + kraken_out, shell=True, check=True)
    subprocess.run('singularity exec --bind ' + out_path + ':/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/kraken2_2.0.8-beta.sif kraken2 --db /kraken2-db/minikraken2_v1_8GB/ --threads ' + threads + ' --use-names --report ' + kraken_out + s + '.report --output ' + kraken_out + s + '_kraken.out --paired ' + input_dir + s + '_1.fastq.gz ' + input_dir + s + '_2.fastq.gz', shell=True, stdout=out_log, stderr=err_log, check=True)
    #Parse Kraken2 output
    with open(kraken_out + s + '.report', 'r') as kreport:
        lines = kreport.readlines()
        for l in lines:
            l_parse = l.lstrip().rstrip().split("\t")
            percent = l_parse[0]
            tax_level = l_parse[3]
            tax = l_parse[5].lstrip()
            if tax_level == 'S':
                break

    #Copy results files to respective directories
    subprocess.run('cp ' + sample_dir + s + '_assembly/' + s + '.fasta assemblies/', shell=True, check=True)
    subprocess.run('cp ' + sample_dir + s + '_assembly/prokka/' + s + '.gff annotations/', shell=True, check=True)
    subprocess.run('cp ' + sample_dir + s + '_assembly/' + s + '_amrfinderplus_report.tsv amrfinder_results/', shell=True, check=True) 

    #Write to output file
    results = [s, genus + '_' + species, accession, distance, tax, percent, scheme, st, num_read, avg_read_len, avg_qual, cov, contigs, long_contig, n50, l50, genome, gc, cds]
    report.write('\t'.join(map(str,results)) + '\n')
    out_log.close()
    err_log.close()

report.close()
