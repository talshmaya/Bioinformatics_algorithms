#!/usr/bin/env python2.7

import os
import argparse
import sys
from collections import OrderedDict

GROUPHOME = os.environ['GROUPHOME']


def get_spec_sen(input_drug, R, S, R_tot, S_tot, gene, mutation, position, isolate,sn_th,sp_th):

    sens = R / float(R_tot) * 100
    spec = 100 - (S / float(S_tot) * 100)
    if sp_th <= spec and sn_th <= sens:
        d = get_lineage(input_drug, isolate)
        lineages, counts = [], []
        for k,v in d.iteritems():
            lineages.append(str(v))

        return '\t'.join(('%.3f' % sens, '%.3f' % spec, gene,mutation, position, '\t'.join([i for i in lineages])))
    else:
        pass


def get_lineage(input_drug, l):
    def get_R_or_S():
        def get_drug():
            dst_file = open(GROUPHOME + '/metadata/dst.txt')
            header = dst_file.readline()
            drugs = header.split()
            col = -1
            for drug in drugs:
                col += 1
                if drug == input_drug:
                    return drug, col
        r_isolates=[]
        drug_col = get_drug()
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

    l = l.split(',')
    lineages = ['Lineage', 'Lineage 7', 'Euro-American', 'M. bovis', 'Cannot be determined',
                'M. Africanum, West African 1', 'Indo-Oceanic', 'M. Africanum, West African 2',
                'East-Asian', 'East-African-Indian']
    d = OrderedDict()
    for lineage in lineages:
        d[lineage] = 0

    for line in open(GROUPHOME + '/metadata/lineage.txt'):
        li = line.strip().split()
        isolate = li[0]
        lineqage = li[1]
        if lineqage == 'Lineage':
            lineqage = 'Cannot be determined'
        for iso in l:
            if iso == isolate and iso in get_R_or_S():
                d[lineqage] = d.get(lineqage, 0) + 1
    return d


def total(drug):
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
    parser.add_argument('-m', '--mutation-file', dest='mutation_file', required=True,
                        help='File name for the mutation file that was created for the drug of interest.'
                             'Include the full file path')
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'),
                        help='Output filename. Writes to standard out by default')
    parser.add_argument('-sn', '--sensitivity_threshold', default=0.0, type=float,
                        help='specify the sensitivity threshold (0-100%) ')
    parser.add_argument('-sp', '--specificity_threshold', default=0.0, type=float,
                        help='specify the specificity threshold (0-100%) ')
    args = parser.parse_args()
    dst = total(args.drug)
    output = args.output
    global drug
    drug=args.drug
    if drug=="resistance":
        R_total=len([r for r in dst.values() if r == 'MDR' or r == 'preXDR' or r == 'XDR' ])
        S_total = len([s for s in dst.values() if s == 'panS'])
    else:
        R_total = len([r for r in dst.values() if r == 'R'])
        S_total = len([s for s in dst.values() if s == 'S'])
    #R_total=37 #mono R PZA
    # genes count GWA output file name
    
    mut_file = args.mutation_file

    with open(mut_file) as f:
        lines=f.readlines()[1:]
    output.write('\t'.join(('Sensitivity', 'Specificity', 'gene', 'mutation', 'position', 'Lineage', 'Lineage 7', 'Euro-American', 'M. bovis', 'Cannot be determined',
                'M. Africanum, West African 1', 'Indo-Oceanic', 'M. Africanum, West African 2',
                'East-Asian', 'East-African-Indian')) + '\n')
    out_list = []
    for line in lines:
        line=line.strip()
        var=line.split()
        if len(var) > 6:
            out = get_spec_sen(input_drug=args.drug,
                               R=int(var[4]),
                               S=int(var[5]),
                               R_tot=R_total,
                               S_tot=S_total,
                               gene=var[0],
                               mutation=var[1],
                               position=var[3],
                               isolate=var[6], 
                               sn_th=args.sensitivity_threshold,
                               sp_th=args.specificity_threshold)
            out_list.append(out)
    
    for line in out_list:
        if line:
            output.write(line + '\n')


if __name__ == '__main__':
    main()
