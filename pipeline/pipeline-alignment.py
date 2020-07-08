#################################################################
#################################################################
############### Genome Alignment ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Guccione Laboratory
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import ruffus.cmdline as cmdline
import sys
import os
import glob

##### 2. LSF #####
# 2.1 Import
sys.path.append('/hpc/users/torred23/pipelines/support')
import lsf

# 2.2 Default parameters
r_source = 'pipeline/scripts/alignment.R'
py_source = 'pipeline/scripts/Alignment.py'
P = 'acc_GuccioneLab'
q = 'express'
W = '00:30'
GB = 5
n = 1
mkdir_val = True

# 2.3 Wrappers
# CMD
def run_job(cmd_str, outfile, W = W, GB = GB, n = n, **kwargs):
	lsf.run_job(cmd_str, outfile, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

# R
def run_r_job(func_name, func_input, outfile, W = W, GB = GB, n = n, **kwargs):
	lsf.run_r_job(func_name, func_input, outfile, r_source=r_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

# Py
def run_py_job(func_name, func_input, outfile, W = W, GB = GB, n = n, **kwargs):
	lsf.run_py_job(func_name, func_input, outfile, py_source=py_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

##### 3. Custom script imports #####
# 3.1 Python
#sys.path.append('pipeline/scripts')
#import Alignment as P

# 3.2 R
# r.source(r_source)

#############################################
########## 2. General Setup
#############################################
##### 1. Pipeline running #####
# Pipeline args
options = cmdline.get_argparse().parse_args()

##### 2. Variables #####
# Dictionary of genomes and URLs
genome_dict = {
    "mouse": [
        {
            "assembly": "mm10",
            "synonym": "GRCm38",
            "primary_assembly_url": "ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz",
			"cdna_url": "ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz",
            "gtf_url": "ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz",
        }
    ],
    "human": [
        {
            "assembly": "hg38",
            "synonym": "GRCh38",
            "primary_assembly_url": "ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
			"cdna_url": "ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
            "gtf_url": "ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz"
        }
    ]
}

#######################################################
#######################################################
########## S1. Download data
#######################################################
#######################################################

#############################################
########## 1. Download
#############################################

def downloadJobs():
	for organism, genomes in genome_dict.items():
		if organism in ['human']:
			for genome in genomes:
				outdir = 'hydra/{assembly}/ensembl/{assembly}.log'.format(**locals(), **genome)
				yield [None, outdir, organism, genome]

@files(downloadJobs)

def downloadGenomes(infile, outfile, organism, genome):

	# Create directory
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Download
	os.system('wget -P {outdir} {primary_assembly_url} && wget -P {outdir} {gtf_url} && wget -P {outdir} {cdna_url} &> {outfile}'.format(**locals(), **genome))

#############################################
########## 2. Unzip genomes
#############################################

@follows(downloadGenomes)

@transform('hydra/*/ensembl/*.gz',
		   suffix('.gz'),
		   '')

def unzipGenomes(infile, outfile):

	# Unzip
	os.system('gunzip --verbose {infile}'.format(**locals()))

#############################################
########## 3. Zip genomes
#############################################

@follows(downloadGenomes)

@transform(glob.glob('hydra/*/ensembl/*.gtf') + glob.glob('hydra/*/ensembl/*.fa'),
		   suffix(''),
		   '.gz')

def zipGenomes(infile, outfile):

	# Zip
	os.system('gzip --verbose {infile}'.format(**locals()))

#######################################################
#######################################################
########## S2. Build indices
#######################################################
#######################################################

#############################################
########## 1. STAR
#############################################

@follows(downloadGenomes, unzipGenomes)

@transform('hydra/*/ensembl/*.gtf',
		   regex(r'(.*hg.*)/ensembl/.*.gtf'),
		   add_inputs(r'\1/ensembl/*primary_assembly.fa'),
		   r'\1/STAR/')

def buildStarIndex(infiles, outfile):
	
	# Command
	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {infiles[1]} --sjdbGTFfile {infiles[0]} --runThreadN 12 --outFileNamePrefix {outfile}'''.format(**locals())

	# Create log dir
	logdir = os.path.join(outfile, 'job')
	if not os.path.exists(logdir):
		os.makedirs(logdir)

	# Get log files
	log_files = {x: os.path.join(logdir, 'job.')+x for x in ['stdout', 'stderr', 'lsf']}

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.3a'], W='02:00', GB=3, n=12, ow=True, **log_files)

#############################################
########## 2. kallisto
#############################################

@follows(downloadGenomes, zipGenomes)

@transform('hydra/*/ensembl/*cdna.all.fa.gz',
		   regex(r'(.*)/ensembl/(.*).fa.gz'),
		   r'\1/kallisto/\2.idx')

def buildKallistoIndex(infile, outfile):

	# Command
	cmd_str = '''kallisto index -i {outfile} {infile}'''.format(**locals())

	# Create log dir
	logdir = os.path.join(os.path.dirname(outfile), 'job')
	if not os.path.exists(logdir):
		os.makedirs(logdir)

	# Get log files
	log_files = {x: os.path.join(logdir, 'job.')+x for x in ['stdout', 'stderr', 'lsf']}

	# Run
	run_job(cmd_str, outfile, modules=['kallisto/0.46.1'], W='01:00', GB=24, n=1, ow=True, **log_files)


# STAR_INDEX="/Users/maayanlab/data/starindex/mouse_ensemble_90"

# GENOME="/Users/maayanlab/data/hydra/Mus_musculus.GRCm38.dna.chromosome.fa"
# ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/

# GTF="/Users/maayanlab/data/hydra/Mus_musculus.GRCm38.90.gtf"
# ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/

# ~/OneDrive/star/STAR \
# --runMode genomeGenerate \
# --genomeDir $STAR_INDEX \
# --genomeFastaFiles $GENOME \
# --sjdbGTFfile $GTF \
# --runThreadN 8 \
# --sjdbOverhang 100


# # Command
# cmd_str = '''STAR \
# 	--quantMode TranscriptomeSAM GeneCounts \
# 	--genomeDir {infiles[1]} \
# 	--readFilesIn {infiles[0]} \
# 	--readFilesCommand zcat \
# 	--outFileNamePrefix {prefix} \
# 	--runThreadN 8 \
# 	--outSAMstrandField intronMotif \
# 	--outFilterIntronMotifs RemoveNoncanonical \
# 	--outSAMtype BAM SortedByCoordinate \
# 	--limitIObufferSize 45000000
# '''.format(**locals())

# cmd_str = '''kallisto quant -i {infiles[0][1]} -o {outfile} {infiles[0][0]} {infiles[1][0]} '''.format(**locals())

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
if __name__ == '__main__':
	cmdline.run(options)
print('Done!')