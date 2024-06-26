#!/usr/bin/env python3
"""
Python script used to identify structuraly identical transcripts.

Description:
This script selects structuraly identical transcripts from RNA sequencing libraries,
outputs a fasta file with all identical mRNA and lncRNAs, and creates a CSV file with information about 
exons, fpkm, tpm, coverage and transcript length.

version: 1.1
"""

import os
import statistics
from Bio import SeqIO
import argparse

# Constant variables
WORK_DIR=os.getcwd() + '/'

def get_info(_group=[], _count=3):
    # If any transcripts were expressed -skip
    if (len([g for g in _group if '-' in g]) <= _count):
        # Exon - mode, fpkm & tpm - median, cov - mean, len - mode 
        exons, fpkms, tpms, coverages, lengths = [],[],[],[],[]
        
        for each in _group:
            # if the transcript was expressed in the experiment
            if '-' not in each:
                # exons
                exons.append(int(each.split('|')[2]))
                fpkms.append(float(each.split('|')[3]))
                tpms.append(float(each.split('|')[4]))
                coverages.append(float(each.split('|')[5]))
                lengths.append(int(each.split('|')[6]))
        # Exon, FPKM, FPKM STD, TPM, TPM STD, coverage, coverage STD, length, length STD
        if _count < 2:
            return [statistics.mode(exons), statistics.mean(fpkms), statistics.stdev(fpkms), statistics.mean(tpms), statistics.stdev(tpms), statistics.mean(coverages), statistics.stdev(coverages), statistics.mean(lengths), statistics.stdev(lengths)]
        else:
            return [statistics.mode(exons), statistics.mean(fpkms), None, statistics.mean(tpms), None, statistics.mean(coverages), None, statistics.mean(lengths), None]
    
    return None

# Transcripts without validation
def get_tracking(_tracking_file, _code=['u','i','x','='], _group={}, _exp=3):
    # Keys:
    #   tracking file - it is the stringtie file
    #   info - exons, TPM, length. total of isoforms
    #   groups: [ex1, ex2, ex3, ...]
    dic_ret = {}
    exp = 0 if _exp == 3 else 1 if _exp == 2 else 2

    dic_group = dict()
    # Correct base 0 in Python
    for k,v in _group.items():
        dic_group[k] = [x + 3 for x in v]
    
    with open(_tracking_file) as out:
        # Keep all the info
        dic_tmp = dict()
        # Keep temporary experiments and groups
        
        for line in out.read().splitlines():
            # frag, locus, gene, code, exp*
            #    0,     1,    2,    3, ...
            param = []
            for each in line.split('\t'):
                    param.append(each)
            
            if param[3] in _code:
                for k, v in dic_group.items():
                    group = []
                    # selecting groups from the dict group: q1,q2,q12 / q3,q4,q5 ... by position
                    for each in v:
                        group.append(param[each])

                    dic_tmp[k] = group

                # Processing the experiments 35: [q3,q4,q5], 40: [q6,q7,q8]...
                dic_info = {}
                
                for k, group in dic_tmp.items():
                    # If the group doesn't have any transcript, skip
                    dic_info[k] = get_info(group, exp)
                
                # Check if all values in the dict is None
                tmp_dict = {k: v for k, v in dic_info.items() if v is not None}
                
                if len(tmp_dict) > 0:
                    dic_ret[(param[0],param[1],param[3],param[2])] = dic_info

    return dic_ret

# Count all transcripts from tracking
def count_transcript(tmp, typ):
    dic_exp = {'35': 0, '40': 0, '45': 0, '50': 0}
    dict_u = {key: tmp[key] for key in tmp.keys() if typ in key}

    for _,v in dict_u.items():
        for key, value in v.items():
            if value != None:
                dic_exp[key] = dic_exp[key] + 1

    return dic_exp

# Get only one transcript type in a dict
def get_only_transcripts(tmp, typ):
    # Get only specific transcripts
    dic = {key: tmp[key] for key in tmp.keys() if (typ in key)}
    return dic

# return locus from a transfrag
def get_locus(trans, all_transcripts):
    for k, _ in all_transcripts.items():
        if trans in k:
            return k[1]

def get_validated_symbols():
    tmp = {'pcg':[], 'lncrna':[]}
    info = get_validated()
    for k in info['='].keys():
        tmp['pcg'].append(k[3].split('|')[0].split('-')[1])
    for k in info['u'].keys():
        tmp['lncrna'].append(k[1])
    for k in info['x'].keys():
        tmp['lncrna'].append(k[1])
    for k in info['i'].keys():
        tmp['lncrna'].append(k[1])
        
    return tmp

def get_validated_symbols_temp():
    tmp = {'35': [], '40': [], '45': [], '50': []}
    for _,v in get_validated().items():
        for key, val in v.items():
            for each in ['35','40','45','50']:
                if val[each] is not None:
                    if key[2] == '=':
                        tmp[each].append(key[3].split('|')[0].split('-')[1])
                    else:
                        tmp[each].append(key[1])
        
    return tmp


def get_validated():
    for each in ['u','i','x','=']:
        print(each, count_transcript(get_validated()[each], each))

# Read the tracking transcritps into a dict of biopython sequences - valid transcripts
def get_fasta():
    with open(WORK_DIR + '/stringTie/track/ALL_merged_track.fasta') as tmp:
        seq_dic = SeqIO.to_dict(SeqIO.parse(tmp, "fasta"))
    return seq_dic

# Get validated transcripts after Blastx, InterProScan and CPC2 - only unique transcripts
def get_validated(all_transcripts, ids):
    # Check in the file for x,i,u and = transcripts
    put_lnc_anti = [x.strip() for x in open(WORK_DIR + '/lncrna/antisense_putative_lncrnas_ids')]
    put_lnc_inter = [x.strip() for x in open(WORK_DIR + '/lncrna/intergenic_putative_lncrnas_ids')]
    put_lnc_intra = [x.strip() for x in open(WORK_DIR + '/lncrna/intragenic_putative_lncrnas_ids')]

    dic_inter = get_only_transcripts(all_transcripts, 'u')
    dic_antis = get_only_transcripts(all_transcripts, 'x')
    dic_intra = get_only_transcripts(all_transcripts, 'i')
    dic_mrna = get_only_transcripts(all_transcripts, '=')

    dic_inter_tmp = dict(dic_inter)
    dic_antis_tmp = dict(dic_antis)
    dic_intra_tmp = dict(dic_intra)

    for k,v in dic_inter.items():
        if k[0] not in put_lnc_inter:
            del dic_inter_tmp[k]

    for k,v in dic_antis.items():
        if k[0] not in put_lnc_anti:
            del dic_antis_tmp[k]

    for k,v in dic_intra.items():
        if k[0] not in put_lnc_intra:
            del dic_intra_tmp[k]

    return {'=':dic_mrna, 'u': dic_inter_tmp, 'x': dic_antis_tmp, 'i': dic_intra_tmp}

def save_output():
    None


def file_exists(filepath):
    """Check if the file exists and return the file path if it does."""
    if not os.path.isfile(filepath):
        raise argparse.ArgumentTypeError(f"The file {filepath} does not exist. Please, verify!")
    return filepath

def check_molecules(lst):
    """Check the transcript classification codes"""
    stringtie_elements = ['=','c','k','m','n','j','e','o','s','x','i','y','p','r','u']
    elements = lst.strip('[]').split(',')
    result = [element.strip() for element in elements]
    if not all(item in stringtie_elements for item in result):
        raise argparse.ArgumentTypeError(f"The following element(s) were not recognized: {set(result) - set(stringtie_elements)}. Please, verify!")
    return result

def parse_groups(grs):
    import json
    try:
        return json.loads(grs)
    except json.JSONDecodeError as e:
        raise argparse.ArgumentTypeError(f"Invalid element(s): {e}")
    
def check_groups(grs):
    try:
        items = grs.split('/')
        result = {}
        for item in items:
            key, value = item.split(':')
            print(key,value)
            key = key.strip()
            value = list(value.strip('[]').replace(',', ''))
            result[key] = value
        return result
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Error parsing input string: {e}")    
    

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Identify structurally identical transcripts.")
    # putative lncRNA file - containing the IDs of all putative lncRNAs
    parser.add_argument('-i', '--ids', type=file_exists, help='File containing only transfrag IDs that were previsouly checked - putative transcripts')
    
    # Transcripts default [u,i,x,=]
    parser.add_argument('-m', '--molecules', type=check_molecules, required=True, default=['u','x','i','='], help='Which molecules will be selected to be analyzed')
    
    # Groups [q1,q2,q3][q4,q5,q6]
    # it will name ex1, ex2, ex3 according to the total of brackets
    parser.add_argument('-g', '--groups', type=check_groups, required=True, help='Group your samples in the same order they were processed by StringTie')

    # Tracking File argument
    parser.add_argument('-t', '--tracking', type=file_exists, help='Tracking file from StringTie2 output to be read')
    # Fasta File argument
    parser.add_argument('-f', '--fasta', type=file_exists, help='Fasta file from StringTie2 output to be used to save identical transcripts fasta file')
    # Output csv file
    parser.add_argument('-o', '--output', type=str, default='output.csv', help='Output file (default: output.csv)')
    # it will return all identical mRNAs
    args = parser.parse_args()
    return args


def main():
    """ Main function """
    # Parse command-line arguments
    args = parse_arguments()
    print(args.molecules)
    print(args.groups)

    group = {'g1':[1,2,3],'g2':[4,5,6]}
    # _tracking_file, _code=['u','i','x','='], _group={}, _exp=3
    #all_transcripts = get_tracking(args.tracking, _group=group, _exp=3)
    #print(all_transcripts)


    
if __name__ == "__main__":
    main()