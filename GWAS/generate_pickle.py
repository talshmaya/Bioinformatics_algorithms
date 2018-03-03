#!/usr/bin/env python
# $ ./generate_pickle.py.py -h

from collections import Counter, defaultdict
import sys
import os
import fileinput
import argparse
import gzip
GROUPHOME = os.environ['GROUPHOME']
import pickle
import itertools


def count_variation(dst_file, bed_range):
    gene_AA = tuple()  # tuple that has the gene and amino acid
    info_per_line = []
    AA_isolate_gene_list = []
    clonality_dictionary = dict()

    MUMMER_PATH = '/data/variants/mummer-denovo/'
    PBHOOVER_PATH = '/data/variants/pbhoover/pbhooverV1.0.0a5/'
    OTHER_PATH = '/data/variants/public-genomes/'  
    extension = '.vcf.gz'
    
    vep_full_path = ""  # this variable will hold the full path to the vep file

    # Iterate through all isolate.high-impact.vcf.gz files
    for text in isolate_file:  # need to add a try block, will fail if file doesn't exist
        if "1-0028" in text:
            continue
        isolate = text.strip()
        if isolate + \
                extension in os.listdir(GROUPHOME + PBHOOVER_PATH):
            vep_full_path = GROUPHOME + PBHOOVER_PATH + isolate + extension
        elif isolate + extension in os.listdir(GROUPHOME + MUMMER_PATH):
            vep_full_path = GROUPHOME + MUMMER_PATH + isolate + extension
        elif isolate + extension in os.listdir(GROUPHOME + OTHER_PATH):
            vep_full_path = GROUPHOME + OTHER_PATH + isolate + extension
        else:
            print "couldn't find vep file for isolate: " + str(isolate)
            continue

        for line in gzip.open(vep_full_path):
            columns = line.split()
            if "#" in line:
                pass
            else:  # read each line that is not header
                AA_isolate_gene = ('', '', '', '', '')
                bb = ()
                isolate_gene = tuple()
                combined_info_lines, combined_info_lines_mut = [], []
                del combined_info_lines[:]
                del combined_info_lines_mut[:]
                del dup_mut[:]
                position = columns[1]
                ref = columns[3]
                alt = columns[4]
                infos = columns[7]
                CSQs = infos.split('CSQ=')

                if "|" in CSQs[1] and (int(position) in bed_range):
                    CSQ = str(CSQs[1]).split(',')
                    num_of_csqs = len(CSQ)
                    prev_cnq, prev_gene = "", ""

                    for dd in CSQ:
                        mutation = ""
                        # parse data from.high-impact.vcf.gz consequence line
                        info_per_line = dd.split('|')
                        AA = info_per_line[15]
                        gene = info_per_line[4]
                        gene_name = info_per_line[3]
                        if gene_name != "":  # if the gene name field is populated, use that, if not, get the gene
                            gene = gene_name
                        Consequence = info_per_line[1]
                        IMPACT = info_per_line[2]
                        DISTANCE = info_per_line[18]
                        cDNA_position = info_per_line[12]
                        PICK = info_per_line[19]
                        STRAND = info_per_line[17]
                        VARIANT_CLASS = info_per_line[22]
                        CODON = info_per_line[14]
                        BIOTYPE = info_per_line[7]
                        if "/" in AA:  # split mutation to ref and alt
                            AA_split = AA.split('/')
                            AA_ref = AA_split[0]
                            AA_alt = AA_split[1]

                        elif AA != "":
                            AA_ref = AA
                            AA_alt = ""
                        elif AA == "":
                            AA_ref = columns[3]
                            AA_alt = columns[4]

                        # parse mutation according to whether its in a gene or
                        # not and if its a snp/indel
                        isolate_gene = (gene, isolate)
                        if "upstream" not in Consequence and "downstream" not in Consequence:  # variant in a gene
                            if "SNV" in VARIANT_CLASS:  # its a snp
                                if "RNA" in BIOTYPE or "pseudogene" in BIOTYPE:
                                    mutation = AA_ref + cDNA_position + AA_alt
                                    AA_isolate_gene = (
                                        mutation, gene, isolate, position, Consequence)
                                    AA_isolate_gene_list.append(
                                        AA_isolate_gene)
                                else:
                                    mutation = AA_ref + CODON + AA_alt
                                    AA_isolate_gene = (
                                        mutation, gene, isolate, position, Consequence)
                                    AA_isolate_gene_list.append(
                                        AA_isolate_gene)
                            elif "substitution" in VARIANT_CLASS:
                                if "pseudogene" in BIOTYPE or "RNA" in BIOTYPE:
                                    mutation = AA_ref + cDNA_position + AA_alt
                                    AA_isolate_gene = (
                                        mutation, gene, isolate, position, Consequence)
                                    AA_isolate_gene_list.append(
                                        AA_isolate_gene)
                                else:
                                    if AA == "" and "upstream" not in prev_cnq:  # same variants, different genes
                                        dnap = cDNA_position.split("-")
                                        mutation = str(ref) + \
                                            str(dnap[-1]) + str(alt)
                                        AA_isolate_gene = (mutation + "/" + prev_mutation, gene + "/" + prev_gene, isolate, position,
                                                           Consequence)
                                        AA_isolate_gene_list.append(
                                            AA_isolate_gene)
                                        AA_isolate_gene = (
                                            mutation, gene, isolate, position, Consequence)
                                        AA_isolate_gene_list.append(
                                            AA_isolate_gene)
                                        dup_mut.append((gene, mutation))
                                        dup_mutation.append(
                                            [x for x in dup_mut])

                                    spl_cod = ""
                                    if "-" in CODON:
                                        spl_cod = CODON.split("-")
                                        start_cod = int(spl_cod[0])
                                    else:
                                        start_cod = CODON
                                    i = 0
                                    for rf, al in zip(AA_ref, AA_alt):

                                        if rf != "" and al != "":
                                            cod_pos = int(start_cod) + i
                                            i += 1
                                            mutation = str(
                                                rf) + str(cod_pos) + str(al)
                                            dup_mut.append((gene, mutation))
                                    mutation = (
                                        "+".join([x[1] for x in dup_mut]))
                                    AA_isolate_gene = (
                                        mutation, gene, isolate, position, Consequence)

                                    AA_isolate_gene_list.append(
                                        AA_isolate_gene)
                            elif "insertion" in VARIANT_CLASS or "deletion" in VARIANT_CLASS or "indel" in VARIANT_CLASS or "protein_altering" in VARIANT_CLASS:  # its an indel
                                if "-" in cDNA_position:
                                    # if it's a range, we only want the start
                                    # position
                                    ll = cDNA_position.split("-")
                                    cdna_pos = ll[0]
                                else:
                                    cdna_pos = cDNA_position

                                if "," in alt:
                                    alt_all = []
                                    sp = alt.split(",")
                                    for S in sp:
                                        alt_all.append(S)
                                    alt = ";".join(alt_all)

                                mutation = str(ref) + \
                                    str(cDNA_position) + str(alt)

                                AA_isolate_gene = (
                                    mutation, gene, isolate, position, Consequence)
                                AA_isolate_gene_list.append(AA_isolate_gene)

                        # elif int(DISTANCE) <= MAX_DIST:# and not args.bed: #
                        # not in a gene but within 200 BP of a gene

                        # comment this block out and uncomment the line before
                        # to go back to original function

                        # only consider upstream conseq if it's the only one,
                        # pick the one with smallest distance
                        elif AA_isolate_gene[0] == "":
                            dist_per_csqs = []
                            # append all distances per consq
                            c = 0
                            min_dist = 10000000
                            for csq_dist in CSQ:
                                c += 1  # keep count of which consq this is
                                dist_per_line = csq_dist.split('|')
                                if dist_per_line[18] != "" and "upstream" in dist_per_line[1]:
                                    # if upstream distance is less than prev
                                    if dist_per_line[18] < min_dist:
                                        min_dist = dist_per_line[18]
                                        min_index = c
                                        gene = dist_per_line[4]
                                        gene_name = dist_per_line[3]
                                        if gene_name != "":  # if the gene name field is populated, use that, if not, get the gene
                                            gene = gene_name
                                        DISTANCE = dist_per_line[18]

                        ###########

                            mutation = str(ref) + "-" + \
                                str(DISTANCE) + str(alt)
                            gene = str(gene) + "_prom"
                            AA_isolate_gene = (
                                mutation, gene, isolate, position, Consequence)
                            AA_isolate_gene_list.append(AA_isolate_gene)
                            break

                        prev_cnq = Consequence
                        prev_gene = gene
                        prev_mutation = mutation
        pickle.dump(AA_isolate_gene_list, open(
            '1.' + str(bed_range[0]) + '.' + str(bed_range[-1]) + '.pickle.int', 'wb'), -1)


if __name__ == "__main__":

    MAX_DIST = 200  # max bp distance from a gene to output
    dst_file = open(GROUPHOME + '/metadata/dst.txt')
    header = dst_file.readline()
    # get the list of all the isolates from the dst file
    all_islts_dst, dup_mutation, dup_mut = [], [], []
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Counts mutations/genes for multiple isolates. Input a file with the '
                                                 'list of isolates, a drug name and the what you wish to count by '
                                                 '(mutation/gene/both) and an output file with counts will be generated')
    parser.add_argument("-i", "--input", dest="filename", required=False, action='store',
                        help="input file with isolate names, one isolate per line", metavar="FILE")
    parser.add_argument("-b", "--bed-file", dest="bed", required=False, action='store',
                        help="input a bed file with 3 columns, tab seperated, CHR START STOP", metavar="FILE")

    args = parser.parse_args()

    if args.filename:
        isolate_file = open(args.filename)
    else:
        for line in dst_file:
            t = line.split('\t')
            if "1-0028" not in t:
                all_islts_dst.append(t[0])
            isolate_file = all_islts_dst

    bed_range = []
    # if there's a bed file, get the range of all the positions to check,
    # otherwise use all the positions, need to think of a more efficient way
    # of doing this
    if args.bed:

        bed_file = args.bed
        spos = bed_file.split(".")
        # for POS in bed_file:
        start = int(spos[1])
        stop = int(spos[2])
        """
        bed_file=open(args.bed)
        for POS in bed_file:
            spos = POS.split()
            start=int(spos[1])
            stop=int(spos[2])
        """
        for i in range(start, stop+1, 1):
            bed_range.append(i)
    else:
        for i in range(1, 5000000, 1):
            bed_range.append(i)

    count_variation(dst_file, bed_range)
