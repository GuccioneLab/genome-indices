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
def run_job(cmd_str, outfile, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs):
	lsf.run_job(cmd_str, outfile, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir, **kwargs)

# R
def run_r_job(func_name, func_input, outfile, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs):
	lsf.run_r_job(func_name, func_input, outfile, r_source=r_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir, **kwargs)

# Py
def run_py_job(func_name, func_input, outfile, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs):
	lsf.run_py_job(func_name, func_input, outfile, py_source=py_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir, **kwargs)

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
            "gff3_url": "ftp://ftp.ensembl.org/pub/release-99/gff3/mus_musculus/Mus_musculus.GRCm38.99.gff3.gz"
        }
    ],
    "human": [
        {
            "assembly": "hg38",
            "synonym": "GRCh38",
            "primary_assembly_url": "ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
			"cdna_url": "ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
            "gtf_url": "ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz",
            "gff3_url": "ftp://ftp.ensembl.org/pub/release-99/gff3/homo_sapiens/Homo_sapiens.GRCh38.99.gff3.gz"
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
				outdir = 'arion/{assembly}/ensembl/{assembly}.log'.format(**locals(), **genome)
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

# @follows(downloadGenomes)

@transform('arion/*/ensembl/*gtf.gz',
		   suffix('.gz'),
		   '')

def unzipGenomes(infile, outfile):

	# Unzip
	os.system('gunzip --verbose {infile}'.format(**locals()))

#############################################
########## 3. Zip genomes
#############################################

@follows(unzipGenomes)

@transform(glob.glob('arion/*/ensembl/*.gtf') + glob.glob('arion/*/ensembl/*.fa'),
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

# @follows(downloadGenomes)

@transform('arion/*/ensembl/*.gtf',
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
	run_job(cmd_str, outfile, modules=['star/2.7.3a'], W='02:00', GB=3, n=12, ow=False, **log_files)

# Alignment
'''
# Single-end
STAR \
	--quantMode TranscriptomeSAM GeneCounts \
	--genomeDir genomeDir \
	--readFilesIn read.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix outfile_prefix \
	--runThreadN 8 \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical \
	--outSAMtype BAM SortedByCoordinate \
	--limitIObufferSize 45000000

# Paired-end
STAR \
	--quantMode TranscriptomeSAM GeneCounts \
	--genomeDir genomeDir \
	--readFilesIn read1.fastq.gz read2.fastq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix outfile_prefix \
	--runThreadN 8 \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical \
	--outSAMtype BAM SortedByCoordinate \
	--limitIObufferSize 45000000
'''

#############################################
########## 2. kallisto
#############################################

# @follows(downloadGenomes)

@transform('arion/*/ensembl/*cdna.all.fa.gz',
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
	run_job(cmd_str, outfile, modules=['kallisto/0.46.1'], W='01:00', GB=24, n=1, ow=False, **log_files)

# Alignment
'''
# Single-end
kallisto quant -i index.idx -o outdir/ --single -l 200 -s 20 -t 8 input.fastq.gz &> output.log

# Paired-end
kallisto quant -i index.idx -o outdir/ read1.fastq.gz read2.fastq.gz &> output.log
'''

#############################################
########## 3. Salmon
#############################################

# @follows(downloadGenomes)

@transform('arion/*/ensembl/*cdna.all.fa.gz',
		   regex(r'(.*)/ensembl/(.*).fa.gz'),
		   r'\1/salmon/index')

def buildSalmonIndex(infile, outfile):

	# Command
	cmd_str = '''salmon index -t {infile} -i {outfile}'''.format(**locals())

	# Create log dir
	logdir = os.path.join(os.path.dirname(outfile), 'job')
	if not os.path.exists(logdir):
		os.makedirs(logdir)

	# Get log files
	log_files = {x: os.path.join(logdir, 'job.')+x for x in ['stdout', 'stderr', 'lsf']}

	# Run
	run_job(cmd_str, outfile, modules=['salmon/1.2.1'], W='01:00', GB=8, n=4, ow=False, **log_files)

# Alignment 
'''
# Single-end
salmon quant -i transcripts_index -l A -r reads.fq.gz -p 8 --validateMappings -o transcripts_quant

# Paired-end
salmon quant -i transcripts_index -l A -1 reads1.fq.gz -2 reads2.fq.gz -p 8 --validateMappings -o transcripts_quant
'''
# See https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype for info on -l flag

#######################################################
#######################################################
########## S3. Splicing
#######################################################
#######################################################

#############################################
########## 1. SUPPA
#############################################

@transform('arion/*/ensembl/*.gtf',
		   regex(r'(.*)/ensembl/(.*).gtf'),
		   r'\1/SUPPA/io./\2*.io.',
		   r'\1/SUPPA/{file_format}/\2')

def buildSuppaIndex(infile, outfiles, outfileRoot):

	# Loop through formats
	for file_format in ['ioi', 'ioe']:

		# Get outfile
		outfile_basename = outfileRoot.format(**locals())

		# Create log dir
		logdir = os.path.join(os.path.dirname(outfile_basename), 'job')
		if not os.path.exists(logdir):
			os.makedirs(logdir)

		# Get log files
		log_files = {x: os.path.join(logdir, 'job.')+x for x in ['stdout', 'stderr', 'lsf']}
		
		# Command
		if file_format == 'ioe':
			cmd_str = '''python \$SUPPA_HOME/suppa.py generateEvents -i {infile} -o {outfile_basename} -f {file_format} -e SE SS MX RI FL'''.format(**locals())
		elif file_format == 'ioi':
			cmd_str = '''python \$SUPPA_HOME/suppa.py generateEvents -i {infile} -o {outfile_basename} -f {file_format}'''.format(**locals())

		# Run
		run_job(cmd_str, outfile_basename, W='00:10', modules=['suppa/2.3'], GB=3, n=1, ow=False, mkdir=False, **log_files)

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