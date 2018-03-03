#!/usr/bin/env python2.7

from collections import OrderedDict
import os
import argparse
import sys

GROUPHOME = os.environ['GROUPHOME']


def get_spec_sen(output, R, S, R_tot, S_tot, gene, isolate,sn_th,sp_th):
    sens = R / float(R_tot) * 100
    spec = 100 - (S / float(S_tot) * 100)
    if sp_th <= spec and sn_th <= sens:
        d = get_lineage(isolate)
        lineages, counts=[], []
        for k, v in d.iteritems():
            lineages.append(str(v))
        output.write('\t'.join(('%.3f' % sens, '%.3f' % spec, gene, '\t'.join([i for i in lineages]))) + '\n')
    else:
        pass


def get_lineage(l):
    l = l.split(',')
    d = OrderedDict()
    lineages = ['Lineage', 'Lineage 7', 'Euro-American', 'M. bovis', 'Cannot be determined',
                'M. Africanum, West African 1', 'Indo-Oceanic', 'M. Africanum, West African 2',
                'East-Asian', 'East-African-Indian']
    for lineage in lineages:
        d[lineage] = 0
    for line in open(GROUPHOME + '/metadata/lineage.txt'):
            li = line.strip().split('\t')
            isolate = li[0]
            lineage = li[1]
            if lineage == 'Lineage':
                lineage = 'Cannot be determined'
            for iso in l:
                if iso == isolate and iso in get_R_or_S():
                    d[lineage] = d.get(lineage, 0) + 1
    return d


def get_drug():
    dst_file = open(GROUPHOME + '/metadata/dst.txt')
    header = dst_file.readline()
    dst_header = header.split()
    col = -1
    for column in dst_header:
        col += 1
        if column == drug:
            return column, col


def get_R_or_S():
    r_isolates=[]
    drug_col = get_drug()
    DRUG = drug_col[0]
    COL = drug_col[1]
    for line in open(GROUPHOME + '/metadata/dst.txt'):
        cols = line.split()
        if drug=='resistance':
            if cols[COL]=='MDR' or cols[COL]=='preXDR' or cols[COL]=='XDR':
                r_isolates.append(cols[0])
        else:
            if cols[COL] == 'R':
                r_isolates.append(cols[0])
    return r_isolates


def total():
    dst_dict = {}
    with open(GROUPHOME + '/metadata/dst.txt', 'r') as dst_handle:
        for line in dst_handle:
            column = line.rstrip('\n').split('\t')
            if 'isolate' in column[0]:
                header = column
                continue
            dst_dict[column[0]] = column[header.index(drug)]
    return dst_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--drug', required=True,
                        help='Input drug name, lowercase, abbreviation as seen in the dst.txt file')
    parser.add_argument('-g', '--gene-file', dest='gene_file', required=True,
                        help='File name for the mutation file that was created for the drug of interest.'
                             'Include the full file path')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'),
                        help='Output filename. Writes to standard out by default')
    parser.add_argument('-sn', '--sensitivity_threshold', default=0.0, type=float,
                        help='specify the sensitivity threshold (0-100%) ')
    parser.add_argument('-sp', '--specificity_threshold', default=0.0, type=float,
                        help='specify the specificity threshold (0-100%) ')
    args = parser.parse_args()
    output = args.output
    drug = args.drug
    # Creating drug as a global variable
    global drug
    dst = total()

    if drug=="resistance":
        R_total=len([r for r in dst.values() if r == 'MDR' or r == 'pre-XDR' or r == 'XDR' ])
        S_total = len([s for s in dst.values() if s == 'panS'])
    else:
        R_total = len([r for r in dst.values() if r == 'R'])
        S_total = len([s for s in dst.values() if s == 'S'])
    
    # genes count GWA output file name 
    gene_file = args.gene_file

    with open(gene_file) as f:
        lines = f.readlines()[1:]
    output.write('\t'.join(('Sensitivity', 'Specificity','gene', 'Lineage', 'Lineage 7', 'Euro-American', 'M. bovis', 
                            'Cannot be determined', 'M. Africanum, West African 1', 'Indo-Oceanic', 
                            'M. Africanum, West African 2', 'East-Asian', 'East-African-Indian')) + '\n')
    for line in lines:
        line=line.strip()
        var=line.split()
        if len(var) > 6:
            # use var[1] and[2] for the 'all' columns and var[3] and [4] for the 'unique' columns 
            get_spec_sen(output=output,
                         R=int(var[1]),
                         S=int(var[2]),
                         R_tot=R_total,
                         S_tot=S_total,
                         gene=var[0],
                         isolate=var[8], 
                        sn_th=args.sensitivity_threshold,
                        sp_th=args.specificity_threshold)

if __name__ == '__main__':
    main()
