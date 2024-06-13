#!/usr/bin/env python3
# Version: 1.3

""""
This script processes paired-end raw reads into cleaning, decontamination, alignment, strand detection, assembly.

Author: Roger Silva
Date: 30/03/2024
"""

import argparse, os
import subprocess
import tempfile as tmp

# Global variables
# Prefix for creating the conda environment
ENV = 'lncrna'
# Defining directories
WORK_DIR=os.getcwd() + '/'
TOOL_DIR = WORK_DIR + 'tool/'
OUT_DIR = WORK_DIR + 'output/'

# Resource directories
R_ADAP = WORK_DIR + 'resources/adapters/'
R_GENO = WORK_DIR + 'resources/genome/GCF_902167145.1/'

# Hisat Index directory
R_HIDX = WORK_DIR + 'output/index/'

R_RAWR = WORK_DIR + 'resources/raw/'
R_RNAS = WORK_DIR + 'resources/rRNA/'
R_TRAN = WORK_DIR + 'resources/transcriptome/'

# Output directories
O_CLEAN = WORK_DIR + 'output/clean/'
O_DECON = WORK_DIR + 'output/decont/'
O_ALIGN = WORK_DIR + 'output/hisat/'
O_ASSEM = WORK_DIR + 'output/stringTie/'
O_ASMEM = WORK_DIR + 'output/stringTie/merge'
O_ASMAB = WORK_DIR + 'output/stringTie/abundance'
O_ASMCO = WORK_DIR + 'output/stringTie/coverage'
O_ASMTR = WORK_DIR + 'output/stringTie/track'
O_SALM  = WORK_DIR + 'output/salmon/'
O_QUAS  = WORK_DIR + 'output/rnaquast/'
O_BLAST = WORK_DIR + 'output/blast/'

# Quality trimming sides - r: right only and l: left only or rl = both sides
QTRIM='rl'
# Trim phred quality
TRIMQ=25
# Read must be at least as long as
MINLEN=40
# Intron size
INTRON=124
# Minimum lncRNA size
MINLNC=200
# Threads 
WORKER=15

def call_variables(dir_path, types):
    """
    Defining the raw read files.
    """
    exps = []
    files = []

    # list fastq.gz files only
    for file in os.listdir(dir_path):
        if file.endswith(types):
            files.append(os.path.join(dir_path, file))

    exps = list(set([x.split('_')[0] for x in files if 'ALL' not in x]))

    return exps
    
# ISF = FR, ISR = RF

def run(env, cmd):
    print('CMD', cmd)
    result = subprocess.call('conda run -n ' + env + ' ' + cmd, shell=True)

def clean(param):
    print('Running the cleaning process...')
    cmd = f"bbduk.sh threads={param['workers']} in1={param['r1']} in2={param['r2']} out1={param['out_r1']} out2={param['out_r2']} ref={param['adapters']} ktrim=r hdist=1 tpe tbo qtrim={param['qtrim']} trimq={param['trimq']} minlen={param['minlen']} outm={param['discarted']} stats={param['stats']} bhist={param['summary']}"
    run('hot', cmd)

def decont(param):
    print('Running decontamination process...')
    cmd = f"bbduk.sh threads={param['workers']} in1={param['r1']} in2={param['r2']} out1={param['out_r1']} out2={param['out_r2']} ref={param['rnas']} k=27 stats={param['stats']} &>> {param['output']}"
    run('hot', cmd)

def hisat(param):
    print('Running HISAT2 alignment process...')
    c1 = f"hisat2 -p {param['workers']} --dta --rna-strandness {param['strand']} -a --secondary --no-mixed --no-discordant --novel-splicesite-outfile {param['splice']} -x {param['index']} "
    c2 = f"-1 {param['r1']} -2 {param['r2']} -S {param['sam']} --un {param['unpair']} --al {param['unalig']} "
    c3 = f"--summary-file {param['summary']}"
    cmd = c1 + c2 + c3
    run('hot', cmd)
    
    print('Sorting BAM files...')
    cmd = f"samtools sort -@ 6 -O bam {param['sam']} -o {param['sort']}"
    run('hot', cmd)

    print('Creating Index files...')
    cmd = f"samtools index {param['sort']}"
    run('hot', cmd)

    print('Removing SAM file...')
    cmd = f"rm {param['sam']}"
    run('hot', cmd)

def stringtie(param):
    """
    -j 5 = there should be at least 5 spliced reads that align across a junction
    -c 5 = minimum 5 read coverage allowed to predict transcripts
    -p process / threads
    """
    print('Running StringTie assembly process...')
    # https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
    cmd = f"stringtie -j 5 -c 5 --{param['strand']} -C {param['cover']} -p {param['workers']} -v -m {param['lncrna']} -o {param['gtf']} -G {param['anot']} {param['sort']}"
    run('hot', cmd)

def merge(param):
    """ 
    -g minimum locus gap separation. Reads closer than 50 bp will be merged (default). 
    """
    print('Running StringTie Merge ALL transcripts from ALL experiments...')
    print(param['files'])
    f_tmp = tmp.NamedTemporaryFile(mode='w', delete=False)
    for i in param['files']:
        f_tmp.write(i + '\n')
    f_tmp.seek(0)
    f_tmp.close()

    print('Creating directories...')
    cmd = f"mkdir {O_ASSEM}abundance {O_ASSEM}merge {O_ASSEM}track"
    run('hot', cmd)

    cmd = f"stringtie --merge -m 200 -p {param['workers']} -G {param['anot']} -o {param['gtf_all']} {f_tmp.name}"
    run('hot', cmd)

    print('Evaluating transcriptome assembly process...')
    cmd = f"gffcompare -V --debug -r {param['anot']} -o {param['out']} {param['gtf_all']}"
    run('hot', cmd)

    print('Tracking identical transcripts...')
    cmd = f"gffcompare -V --debug -r {param['anot']} -o {param['out_trk']} -i {f_tmp.name}"
    run('hot', cmd)

    print('Extracting fasta files from tracking transcripts process...')
    cmd = f"gffread -W -F -w {param['fasta_trk']} -g {param['geno']} {param['gtf_trk']}"
    run('hot', cmd)

    print('Extracting fasta files from transcripts process...')
    cmd = f"gffread -W -F -w {param['fasta_mer']} -g {param['geno']} {param['gtf_all']}"
    run('hot', cmd)

    print('Joining transcripts names - easier identification...')
    cmd = f"{TOOL_DIR}/mstrg_prep.pl {param['gtf_all']} > {param['gtf_out']}"
    run('hot', cmd)


def abundace(param):
    print('Measuring transcript abundance...')
    cmd = f"stringtie -e -B -p {param['workers']} -G {param['gtf_out']} -A {param['abund']} -o {param['trans']} {param['bam']}"
    run('hot', cmd)

def salmon(param):
    print('Running salmon to find strandedness...')
    cmd = f"salmon quant -i {param['idx']} -l A -1 {param['r1']} -2 {param['r2']} -g {param['anot']} -p 4 --gcBias --validateMappings --numGibbsSamples 200 --seqBias -o {param['output']}"
    run('hot', cmd)

def rnaquast(param):
    print('Running rnaQUAST quality...')
    cmd = f"{TOOL_DIR}/rnaQUAST-2.2.2/rnaQUAST.py --transcripts {param['fasta']} --reference {param['geno']} --gtf {param['gff']} -o {param['output']}"
    run('hot', cmd)

def check_size(param):
    print('Checking sizes')
    cmd = '#!/bin/bash; for file in *sorted.bam; do echo $(samtools view -F 4 $file |head -n 1000000|cut -f 10|sort -n| wc -c)/1000000|bc done'
    run('hot', cmd)

def index(param):
    """
    Generates index to be used to hisat, star and salmon
    """
    print('Generating salmon index...')
    cmd = f"salmon index -t {param['cdna']} -i {param['out_l']}"
    run('hot', cmd)

    print('Generating HISAT2 index...')
    cmd = f"hisat2_extract_splice_sites.py {param['ano_h']} > {param['out_h']}ss.exon.txt"
    run('hot', cmd)

    cmd = f"hisat2_extract_exons.py {param['ano_h']} > {param['out_h']}exon.exon.txt"
    run('hot', cmd)

    cmd = f"hisat2-build --ss {param['out_h']}ss.exon.txt --exon {param['out_h']}exon.exon.txt -f {param['geno']}  {param['out_h']}T.thermophilus"
    run('hot', cmd)


def call_directory():
    import os
    '''
    If the folder does not exit, create it!
    '''
    dirs = [TOOL_DIR, R_ADAP, R_GENO, R_HIDX, R_RAWR, R_RNAS, O_CLEAN,O_DECON, O_SALM, O_QUAS, O_ASMEM, O_ASMAB, O_ASMCO, O_ASMTR]
    for folder in dirs:
        if not os.path.exists(folder):
            print(f'Creating directory {folder}')
            os.makedirs(folder)

    return None

def call_clean():
    typ = '.fastq.gz'
    dir = R_RAWR
    files = call_variables(dir, typ)

    for file in files:
        param=dict()

        file_name = os.path.basename(file)
        
        param['adapters']   = R_ADAP + 'adapters.fa'
        param['r1']         = file + '_1.fastq.gz'
        param['r2']         = file + '_2.fastq.gz'
        param['out_r1']     = O_CLEAN + file_name + '_1_cleaned.fastq.gz' 
        param['out_r2']     = O_CLEAN + file_name + '_2_cleaned.fastq.gz'
        param['qtrim']      = QTRIM
        param['trimq']      = TRIMQ
        param['minlen']     = MINLEN
        param['discarted']  = O_CLEAN + file_name + '_discarted.fastq'
        param['stats']      = O_CLEAN + file_name + '_stats.out'
        param['summary']    = O_CLEAN + file_name + '_summary.out'
        param['output']     = O_CLEAN + file_name + '.out'
        param['workers']    = WORKER

        clean(param)


def call_decont():
    typ = '.fastq.gz'
    dir = O_CLEAN
    files = call_variables(dir, typ)

    for file in files:
        param=dict()

        file_name = os.path.basename(file)
            
        param['rnas']       = R_RNAS + 'RF00001.fa,' + R_RNAS + 'RF00002.fa,' + R_RNAS + 'RF00005.fa,' + R_RNAS + 'RF01960.fa,' + R_RNAS + 'RF02541.fa,' + R_RNAS + 'RF02543.fa'
        param['r1']         = file + '_1_cleaned.fastq.gz'
        param['r2']         = file + '_2_cleaned.fastq.gz'
        param['out_r1']     = O_DECON + file_name + '_1_cleaned.fastq.gz' 
        param['out_r2']     = O_DECON + file_name + '_2_cleaned.fastq.gz'
        param['stats']      = O_DECON + file_name + '_stats.out'
        param['output']     = O_DECON + file_name + '.out'
        param['workers']    = WORKER

        decont(param)

def call_hisat():
    typ = '.fastq.gz'
    dir = O_DECON
    files = call_variables(dir, typ)

    for file in files:
        param=dict()

        file_name = os.path.basename(file)
            
        param['r1']         = file + '_1_cleaned.fastq.gz'
        param['r2']         = file + '_2_cleaned.fastq.gz'
        param['splice']     = O_ALIGN + file_name + '_splice_sites.out' 
        param['unpair']     = O_ALIGN + file_name + '_unpair.fastq'
        param['unalig']     = O_ALIGN + file_name + '_unpair_aligned.fastq'
        param['summary']    = O_ALIGN + file_name + '_summary.out'
        param['sam']        = O_ALIGN + file_name + '.sam'
        param['bam']        = O_ALIGN + file_name + '.bam'
        param['sort']       = O_ALIGN + file_name + '_sorted.bam'
        param['index']      = R_HIDX + 'Zea'
        param['strand']     = 'FR'
        param['workers']    = WORKER

        hisat(param)

def call_stringtie():
    typ = '_sorted.bam'
    dir = O_ALIGN
    files = call_variables(dir, typ)

    for file in files:
        param=dict()

        file_name = os.path.basename(file).split('_')[0]

        param['cover']      = O_ASSEM + 'coverage/' + file_name + '_covered.gtf'
        param['lncrna']     = MINLNC
        param['gtf']        = O_ASSEM + file_name + '.gtf'
        param['file']       = file_name
        param['anot']       = R_GENO + 'genomic.gff'
        param['sort']       = O_ALIGN + file_name + '_sorted.bam'
        param['abund']      = O_ASSEM + 'abundance/' + file_name +'/' + file_name + '_abund.tab'
        param['strand']     = 'fr'
        param['workers']    = WORKER

        stringtie(param)

def call_merge():
    typ = '.gtf'
    dir = O_ASSEM
    files = call_variables(dir, typ)
    # sort them
    files.sort()
    files = [x for x in files]
    print(files)
    param=dict()

    param['files']      = files
    param['anot']       = R_GENO + 'genomic.gff'
    param['geno']       = R_GENO + 'GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna'
    param['gtf_all']    = O_ASSEM + 'merge/ALL_merged.gtf'
    param['gtf_trk']    = O_ASSEM + 'track/ALL_merged_track.combined.gtf'
    param['gtf_out']    = O_ASSEM + 'merge/ALL_merged_prep.gtf'
    param['out']        = O_ASSEM + 'merge/ALL_merged'
    param['out_trk']    = O_ASSEM + 'track/ALL_merged_track'
    param['fasta_mer']  = O_ASSEM + 'merge/ALL_merged.fasta'
    param['fasta_trk']  = O_ASSEM + 'track/ALL_merged_track.fasta'
    param['workers']    = WORKER
    
    merge(param)


def call_abundance():
    typ = '.gtf'
    dir = O_ASSEM
    files = call_variables(dir, typ)

    for file in files:
        param=dict()

        file_name = os.path.splitext(os.path.basename(file))[0]

        param['abund']         = O_ASSEM + 'abundance/' + file_name +'/' + file_name + '_abund.tab'
        param['gtf_out']       = O_ASSEM + 'track/ALL_merged_track.combined.gtf'
        param['trans']         = O_ASSEM + 'abundance/' + file_name +'/' + file_name + '_combined.gtf'
        param['dir']           = O_ASSEM + 'abundance/' + file_name
        param['bam']           = O_ALIGN + file_name + '_sorted.bam'
        param['workers']       = WORKER

        abundace(param)


def call_rnaquast():
    
    param=dict()
        
    param['gff']        = R_GENO + 'T.thermophilus.gff'
    param['geno']       = R_GENO + 'T.thermophilus.fna'
    param['fasta']      = O_ASSEM + 'ALL_merged.fasta'
    param['output']     = O_QUAS + 'ALL'

    rnaquast(param)


if __name__ == "__main__":

    def dir_path(path):
        if os.path.isdir(path):
            return path
        else:
            raise argparse.ArgumentTypeError(f'Path: {path} is not a valid path')

    parser = argparse.ArgumentParser(description="Options for finding lncRNA pipeline")
    parser.add_argument('-f','--folders', action='store_true', help='Create folder tree for the pipeline to save its processed files')
    parser.add_argument('-c','--clean', action='store_true', help='Clean all raw reads using bbmap tool')
    parser.add_argument('-d','--decont', action='store_true', help='Decontamination of all raw reads')
    parser.add_argument('-i','--hisat', action='store_true', help='Run Hisat for paired-end reads')
    parser.add_argument('-s','--stringtie', action='store_true', help='Run StringTie for all BAM files generated by Hisat2')
    parser.add_argument('-m','--merge', action='store_true', help='Merge transcripts from all samples')
    parser.add_argument('-a','--abundance', action='store_true', help='Count all transcripts from all samples')
    parser.add_argument('-all','--all', action='store_true', help='Run all stages together (clean, decont, hisat, stringtie, merge, abundance)')
    parser.add_argument('-ass','--ass', action='store_true', help='Run only the assembling stages ( stringtie, merge, abundance)')
    parser.add_argument('-q','--quast', action='store_true', help='Run RNAQuast on the assembled meta-transcriptome')
    parser.add_argument('-b','--about', action='store_true', help='Execute the pipeline for the first time')


    args = parser.parse_args()

    if args.about:
        print('Transcriptome assembly pipeline - version: 1.0')

    if args.folders:
        call_directory()

    if args.clean:
        call_clean()

    if args.decont:
        call_decont()

    if args.hisat:
        call_hisat()

    if args.stringtie:
        call_stringtie()

    if args.merge:
        call_merge()

    if args.abundance:
        call_abundance()

    if args.quast:
        call_rnaquast()

    if args.debug:
        typ = '.fastq.gz'
        dir = R_RAWR
        files = call_variables(dir, typ)

    if args.all:
        call_clean()
        call_decont()
        call_hisat()
        call_stringtie()
        call_merge()
        call_abundance()

    if args.ass:
        call_stringtie()
        call_merge()
        call_abundance()

    if args.tools:
        call_rnaquast()
