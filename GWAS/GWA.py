#!/usr/bin/env python
# $ ./GWA.py -h

from collections import Counter, defaultdict
import sys
import os
import errno
#import fileinput
import argparse
#import gzip
GROUPHOME = os.environ['GROUPHOME']
import pickle
import glob
from os import listdir
from os.path import isfile
import itertools
import time
import datetime

def count_variation(dst_file, bed_range,isolate_file,arg_drug, arg_all, arg_type, arg_output):
    """Looks for a gene_list.pickle file, if it doesn't exist
    Looks for intermidiate pickle files, then merges them into gene_list.pickle
    Passing the dictionary from the final pickle file to the next method 
    """
    AA_isolate_gene_list = []

    final_pickle_file = 'gene_list.pickle'
    vep_full_path = ""  # this variable will hold the full path to the vep file

    if os.path.isfile(final_pickle_file):
        print "Using a pickle file  ! ! ! ->    " + str(final_pickle_file)
        AA_isolate_gene_list = pickle.load(open(final_pickle_file, 'rb'))

    elif glob.glob('*.pickle.int') != "":
        print "Merging and removing interm pickle files into one file ->  " + str(final_pickle_file)

        pickle_files = []
        onlyfiles = glob.glob('*.pickle.int')
        for f in onlyfiles:
            with open(f, 'rb') as f_pickle:
                pickle_files.append(pickle.load(f_pickle))
                os.remove(f.strip())

        AA_isolate_gene_list = list(
            itertools.chain.from_iterable(pickle_files))
        pickle.dump(AA_isolate_gene_list, open('gene_list.pickle', 'wb'))

    return pop_dict(AA_isolate_gene_list, dst_file, isolate_file, bed_range,arg_drug, arg_all, arg_type, arg_output)


def pop_dict(AA_isolate_gene_list, dst_file, isolate_file, bed_range, arg_drug, arg_all, arg_type, arg_output):
    """Populates a list with the mutation, gene, position, isolate and consequence
    and R or S for the drug inputted 
    If --all is used when executing the script, will use all positions and isolates
    otherwise, only the ones passed using the -i and -b flags will be used
    """

    AA_isolates = []
    dict_isolt_gene_mut_all_R, dict_isolt_gene_mut_all_S = defaultdict(
        list), defaultdict(list)
    AA_set_R, AA_set_S = [], []
    ts= time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    #print "time stamp 2: "
    #print st
    for item in AA_isolate_gene_list:
        mutation = item[0]
        gene = item[1]
        position = item[3]
        isolate = item[2]
        consequence_col = item[4]

        if arg_all:
            if arg_drug == "MDR":
                s_or_r = get_MDR(isolate, dst_file)
            elif arg_drug == "MDR+":
                s_or_r = get_MDR_plus(isolate, dst_file)
            else:
                s_or_r = get_R_or_S(isolate, dst_file,arg_drug)
            bb = [mutation, gene, consequence_col, position]
            if s_or_r and 'R' in s_or_r and bb[0]:
                dict_isolt_gene_mut_all_R[isolate].append(bb)
                AA_set_R.append(','.join([mutation, gene, position]))
                AA_isolates.append(
                    (mutation, gene, position, consequence_col, isolate))

            if s_or_r and 'S' in s_or_r and bb[0]:
                dict_isolt_gene_mut_all_S[isolate].append(bb)
                AA_set_S.append(','.join([mutation, gene, position]))
                AA_isolates.append(
                    (mutation, gene, position, consequence_col, isolate))

        else:
            if isolate in isolate_file and int(position) in bed_range:
                if arg_drug == "MDR":
                    s_or_r = get_MDR(isolate, dst_file)
                elif arg_drug == "MDR+":
                    s_or_r = get_MDR_plus(isolate, dst_file)
                else:
                    s_or_r = get_R_or_S(isolate, dst_file, arg_drug)

                bb = [mutation, gene, consequence_col, position]
                if s_or_r and  'R' in s_or_r and bb[0]:
                    dict_isolt_gene_mut_all_R[isolate].append(bb)
                    AA_set_R.append(','.join([mutation, gene, position]))
                    AA_isolates.append(
                        (mutation, gene, position, consequence_col, isolate))

                if s_or_r and 'S' in s_or_r and bb[0]:
                    dict_isolt_gene_mut_all_S[isolate].append(bb)
                    AA_set_S.append(','.join([mutation, gene, position]))
                    AA_isolates.append(
                        (mutation, gene, position, consequence_col, isolate))

        mut_gene_R_all, mut_gene_S_all = [('', '')], [('', '')]

    AA_mut_list_R = get_count_lists(AA_set_R)
    AA_mut_list_S = get_count_lists(AA_set_S)
    
    ts= time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    #print "time stamp 3: "
    #print st
    return output_lines(
        AA_mut_list_S,
        AA_mut_list_R,
        dict_isolt_gene_mut_all_R,
        dict_isolt_gene_mut_all_S,
        AA_isolates, arg_all, arg_type, arg_output)

def get_drug(arg_drug):
    """Returns which column in the dst file the drug is."""
    dst_file = open('dst.txt')
    header = dst_file.readline()
    drugs = header.split()
    col = -1
    for drug in drugs:
        col += 1
        if drug in arg_drug:
            return drug, col
 

def get_R_or_S(isolate, dst_file, arg_drug):
    """Checks the dst file and returns if that isolate is R or S for 
    the given drug.
    """

    drug_col = get_drug(arg_drug)
    DRUG = drug_col[0]
    COL = drug_col[1]
    #for line in open('dst.txt'):
    for line in open(GROUPHOME + '/metadata/dst.txt'):
        cols = line.split()
        if isolate == cols[0]:
            return cols[COL]

def get_MDR_plus(isolate, dst_file):
    """Checks the dst file, resistance column and returns if that 
    isolate is MDR, XDR or pre-XDR or panS for the given drug.
    """

    for line in open(GROUPHOME + '/metadata/dst.txt'):
        cols = line.split()
        if isolate == cols[0]:
            profile=cols[10]
            if profile.lower()=="MDR".lower() or profile.lower()=="XDR".lower() or profile.lower()=="preXDR".lower(): 
                return "R"
            elif profile.lower()=="panS".lower():
                return "S"

def get_MDR(isolate, dst_file):
    """Checks the dst file, resistance column and returns if that 
    isolate is MDR or panS for the given drug.
    """

    for line in open(GROUPHOME + '/metadata/dst.txt'):
        cols = line.split()
        if isolate == cols[0]:
            profile=cols[10]
            if profile.lower()=="MDR".lower(): 
                return "R"
            elif profile.lower()=="panS".lower():
                return "S"


def get_count_lists(count_set):
    """Get counts for genes and isolates and gene,mutation isolates."""

    gene_mut_list = {}
    for mutation in count_set:
        j = 0
        if len(mutation) > 0:
            gene_mut = mutation
            if gene_mut not in gene_mut_list.keys():
                gene_mut_list[gene_mut] = 1
            else:
                gene_mut_list[gene_mut] += 1

    return gene_mut_list


def output_lines(aa_S, aa_R, all_R, all_S, aa_isolates, arg_all, arg_type, arg_output):
    """Get all the relevant counts for mutation/gene/both count
    print to standard out or write the output to a file
    """
    
    ts= time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    #print "time stamp 4: "
    #print st

    mut_w_pos = list(set([','.join(x[0:3]) for x in aa_isolates]))

    gene_mut_count = [("", 0)]
    ii_genes = [('')]
    R, S = [], []
    R_mut, S_mut = [], []
    
    # set of isolates/consequences for isolates/consq columns in output
    output_isolates, output_conseqs = set(), ""  
    
    # gene, mutation, # Rs, # Ss, genomic position, consequence
    count_all = []  

    # gene, unique # of R per gene, unique # of S per gene, # of isolates per 
    # gene that no S has, # of isolates per gene that no R has number of unique 
    # mutations in the gene,  # of mutations per gene that no S has,  # of mutations per gene that no R has
    # this is a list of tuples in the form of (gene,mutation) for all isolates
    # (both R and S)
    count_gene_output = []  
    gene_mutation_position = {}
    
    for mutation_gene_pos in aa_S:
        gene_mutation_position[mutation_gene_pos] = ''
    for mutation_gene_pos in aa_R:
        gene_mutation_position[mutation_gene_pos] = ''

    # iterate through all gene and mutation pairs
    for mutation_gene_pos in gene_mutation_position.keys(
    ):  
        mutation = mutation_gene_pos.split(',')[0]
        gene = mutation_gene_pos.split(',')[1]
        position = mutation_gene_pos.split(',')[2]
        if mutation and gene and position:
            z = all_R.copy()
            z.update(all_S)
            output_isolates = set()
            for k, v in z.iteritems():
                #GI=gene isolate pair 
                for GI in v:
                    if GI[0] == mutation and GI[1] == gene and GI[3] == position:
                        output_conseqs = GI[2]
                        output_isolates.add(k)

            output_isolates = ','.join(output_isolates)

            if mutation_gene_pos in aa_S and mutation_gene_pos in aa_R:
                count_all.append([gene,
                                  mutation,
                                  aa_R[mutation_gene_pos],
                                  aa_S[mutation_gene_pos],
                                  position,
                                  output_conseqs,
                                  output_isolates])
            elif mutation_gene_pos in aa_S and mutation_gene_pos not in aa_R.keys():  # no R's for that gene+mutation
                count_all.append([gene,
                                  mutation,
                                  0,
                                  aa_S[mutation_gene_pos],
                                  position,
                                  output_conseqs,
                                  output_isolates])
            elif mutation_gene_pos in aa_R and mutation_gene_pos not in aa_S.keys():  # no S's for that gene+mutation
                count_all.append([gene,
                                  mutation,
                                  aa_R[mutation_gene_pos],
                                  0,
                                  position,
                                  output_conseqs,
                                  output_isolates])

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
    
    ts= time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    #print "time stamp 5: "
    #print st

    output_list=[]
    # for the mutation argument (or both mutation and gene)
    if (arg_type == "mutation" or arg_type == "both"):
        output_list.append(
                "gene" +
                "\t" +
                "mutation" +
                "\t" +
                "consequence" +
                "\t" +
                "position" +
                "\t" +
                "R" +
                "\t" +
                "S" + '\n')
        for gene_col, mut_col, R_count, S_count, pos_col, consq, isolts_col in count_all:
                if gene_col != "" and mut_col != "":
                    output_list.append(
                        "\t".join(
                            (str(gene_col),
                             str(mut_col),
                                str(consq),
                                str(pos_col),
                                str(R_count),
                                str(S_count),
                                str(isolts_col))) + '\n')
        if arg_output:
            with open('../data/GWAS/' + arg_output + ".mutation_count", 'w') as output_file:
                for i in output_list:
                    output_file.write(i)
        else:
            for i in output_list:
                print i.strip()
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
    
    ts= time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    #print "time stamp 6: "
    #print st
    
    # for the gene argument (or both gene and mutation)
    if (arg_type == "gene" or arg_type == "both"):
        # creating a set of all the genes, to iterate through
        for k, v in all_R.iteritems():
            for i in v:
                if "/" not in i[1]:
                    ii_genes.append(i[1])
                else:
                    [ii_genes.append(x) for x in i[1].split("/")]
        for k, v in all_S.iteritems():
            for i in v:
                if "/" not in i[1]:
                    ii_genes.append(i[1])
                else:
                    [ii_genes.append(x) for x in i[1].split("/")]
        
        # a set of all the unique genes for all the isolates (R and S)
        all_uniq_genes = set(ii_genes)
        ts= time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        #print "time stamp 7: "
        #print st
        
        for i in all_uniq_genes:  # iterate through all the genes and get count for each gene
            if i != "":
                num_of_isolates_R, num_of_isolates_S, num_of_mutations, num_of_mutations_R, num_of_mutations_S, isolt_R_all_count, isolt_R_all_count = 0, 0, 0, 0, 0, 0, 0
                del R[:]
                del S[:]
                num_of_mutations = [x[0] for x in count_all].count(
                    i)  # count how many unique mutations total
                num_of_mutations_R = ([x[0] for x in count_all if x[3] == 0]).count(
                    i)  # if number of S's is 0, count how many mutations in that gene
                # ^ same but for S.
                num_of_mutations_S = (
                    [x[0] for x in count_all if x[2] == 0]).count(i)

                isolates_R_in_gene = set([pairs[0] for k, v in all_R.iteritems(
                ) for pairs in v if pairs[1] == str(i)])  # all S mutations for that gene
                isolates_S_in_gene = set(
                    [pairs[0] for k, v in all_S.iteritems() for pairs in v if pairs[1] == str(i)])
                ts= time.time()
                st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

                R, S = [], []
                for k, v in all_R.iteritems():  # k is the isolate and v is a pair (mutation, gene)
                    for pairs in v:
                        if pairs[1] == i:  # for the gene of interest
                            # append mutation and isolate
                            R.append((pairs[0], k))
                        elif "/" in pairs[1] and i in pairs[1]:
                            R.append((pairs[0], k))
                for k, v in all_S.iteritems():
                    for pairs in v:
                        if pairs[1] == i:
                            S.append((pairs[0], k))
                        elif "/" in pairs[1] and i in pairs[1]:
                            R.append((pairs[0], k))

                # count how many unique isolates for that gene
                isolt_R_all_count = len(set(([r[1] for r in R])))
                isolt_S_all_count = len(set(([s[1] for s in S])))

                # set of all the isolates that has (gene,mutation)
                isolate_list_col = set(([r[1] for r in R])) | set(
                    ([s[1] for s in S]))
                isolate_list_col = ','.join(isolate_list_col)

                bloo = set([k for k, v in all_R.iteritems()
                            for pairs in v if pairs[1] == str(i)])
                Sloo = set([k for k, v in all_S.iteritems()
                            for pairs in v if pairs[1] == str(i)])

                num_of_isolates_R = len(set([vloo for vloo in bloo for r in R if (
                    str(r[1]) == str(vloo) and r[0] not in isolates_S_in_gene)]))
                num_of_isolates_S = len(set([svloo for svloo in Sloo for s in S if (
                    str(s[1]) == str(svloo) and s[0] not in isolates_R_in_gene)]))

                count_gene_output.append(
                    (i,
                     isolt_R_all_count,
                     isolt_S_all_count,
                     num_of_isolates_R,
                     num_of_isolates_S,
                     num_of_mutations,
                     num_of_mutations_R,
                     num_of_mutations_S,
                     isolate_list_col))
        ts= time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        #print "time stamp 9: "
        #print st
        output_list=[]
        output_list.append(
            "Gene" +
            "\t" +
            " Total Number R (unique) [isolate count]" +
            "\t" +
            " Total Number S (unique) [isolate count]" +
            "\t" +
            "# R with mutations that No S has [isolate count, subset of Col2]" +
            "\t" +
            "# S with mutations that No R has [isolate count, subset of Col3]" +
            "\t" +
            "# of unique mutations" +
            "\t" +
            " # Mutations that only R has [mutation count]" +
            "\t" +
            " # Mutations that only S has [mutation count]"+'\n')
        for a, b, c, d, e, f, g, h, isolts_col in count_gene_output:  # unpack and print all the counts
            if a != "" and b != "" and c != "" and e != "" and f != "" and g != "" and h != "":
                output_list.append(
                    "\t".join(
                        (str(a),
                         str(b),
                            str(c),
                            str(d),
                            str(e),
                            str(f),
                            str(g),
                            str(h),
                            str(isolts_col)))+'\n')
        if arg_output:
            with open('../data/GWAS/' + arg_output + ".genes_count", 'w') as output_file:
                for i in output_list:
                    output_file.write(i)
        else:
            for i in output_list:
                print i.strip()

    ts= time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    #print "time stamp 10: "
    #print st\
    return output_list

def check_4_dups(a, b):
    """This is for the case where we have a gene inside a gene and we express 
    them as gene1/gene2 but don't want to also print gene1 and gene2
    """
    dup = []
    #if len(dup_mutation) > 0:  
    #for d_mut in dup_mutation:
    for j_mut in d_mut:
        if j_mut[0] == a and j_mut[1] == b:
            dup.append(j_mut)
    return dup


if __name__ == "__main__":
    """Paeses the command line flags, creats output directory if needed
    and parses the data in the input files, for example the list of isolates
    or list of positions
    """
    # max bp distance from a gene to output
    MAX_DIST = 200  
    dst_file = open('dst.txt')
    header = dst_file.readline()
    # get the list of all the isolates from the dst file
    dup_mutation= [] 
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Counts mutations/genes for multiple isolates. Input a file with the '
                                                 'list of isolates, a drug name and the what you wish to count by '
                                                 '(mutation/gene/both) and an output file with counts will be generated')
    parser.add_argument("-i", "--input", dest="filename", required=False, action='store',
                        help="input file with isolate names, one isolate per line", metavar="FILE")
    parser.add_argument("-o", "--output", help="Directs the output to a name of your choice. "
                                               "Will generate an OUTPUT.genes_count or OUTPUT.mutation_count (or both) "
                                               "based on the type of analysis inside a data/GWAS/ dir. Default is to write to "
                                               "standatd output", action='store', dest='output')
    parser.add_argument("-t", "--analysis-type", help="Either 'gene', 'mutation' or 'both'. Default is both ",
                        action='store', dest='type_analysis', required=False, default='both')
    parser.add_argument("-r", "--resistance", help="Choose one of the drug names as it appears on the dst table's header ",
                        action='store', dest='resistance',choices=['inh','rif','amk','cap','kan','mox','ofx','fq','pza','MDR','MDR+'] ,required=True)
    parser.add_argument("-b", "--bed-file", dest="bed", required=False, action='store',
                        help="input a bed file with 3 columns, tab seperated, CHR START STOP", metavar="FILE")
    parser.add_argument('-A', '--all', required=False, action='store_true',
                        help='Using this flag will run the GWA on all isolates and all positions')
    args = parser.parse_args()

    # if the output flag is used, create a data and a GWAS directory to put
    # the data inside
    arg_all=args.all
    arg_output=args.output
    arg_type=args.type_analysis
    arg_drug=args.resistance
    if args.output:
        try:
            os.makedirs('../data/GWAS/')
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise exc
            pass
    isolates=[]
    if args.filename:
        isolate_file = open(args.filename)
        for text in isolate_file:  # need to add a try block, will fail if file doesn't exist
            if "1-0028" in text:
                continue
            isolates.append(text.strip())
    else:
        for line in dst_file:
            t=line.split('\t')
            if "1-0028" not in t:
                isolates.append(t[0])

    bed_range = []
    # if there's a bed file, get the range of all the positions to check,
    # otherwise use all the positions, need to think of a more efficient way
    # of doing this
    if args.bed:

        bed_file = open(args.bed)
        for POS in bed_file:
            spos = POS.split()
            start = int(spos[1])
            stop = int(spos[2])

            for i in range(start, stop, 1):
                bed_range.append(i)
    else:
        for i in range(1, 5000000, 1):
            bed_range.append(i)
    
    ts= time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    #print "time stamp 1: "
    #print st
    
    print count_variation(dst_file, bed_range, isolates,arg_drug, arg_all, arg_type, arg_output)