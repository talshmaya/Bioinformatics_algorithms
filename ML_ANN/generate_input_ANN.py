#!/usr/bin/env python
# Tal Shmaya May 12, 2016
# $ ./GWA.py -h

from collections import Counter, defaultdict
from difflib import Differ
import sys, os, errno
import fileinput
import argparse
import gzip
GROUPHOME = os.environ['GROUPHOME']


def count_variation(dst_file, bed_range, isolate_file):

    gene_AA=tuple() #tuple that has the gene and amino acid
    info_per_line=[]
    gene_set_S=set()
    gene_set_R=set()
    AA_set_S = []
    AA_set_R = []
    AA_isolates=[]
    dict_isolt_gene_mut_all_R,dict_isolt_gene_mut_all_S = defaultdict(list),defaultdict(list)
    clonality_dictionary = dict()

    MUMMER_PATH = '/data/variants/mummer-denovo/'
    PBHOOVER_PATH = '/data/variants/pbhoover/pbhooverV1.0.0a5/'
    OTHER_PATH = '/data/variants/public-genomes/' # ?
    mummer_extension='.annotated.high-impact.vcf.gz'
    pbhoover_extension= '.vcf.gz'
    other_extension = '.vcf.gz' # ?

    vep_full_path="" # this variable will hold the full path to the vep file
    #Iterate through all isolate.high-impact.vcf.gz files
    for text in isolate_file: # need to add a try block, will fail if file doesn't exist
        isolate=text.strip()
        s_or_r= get_R_or_S(isolate,dst_file)
        if isolate+pbhoover_extension in os.listdir(GROUPHOME+PBHOOVER_PATH):
            vep_full_path = GROUPHOME+PBHOOVER_PATH+isolate+pbhoover_extension
        elif isolate+mummer_extension in os.listdir(GROUPHOME+MUMMER_PATH):
            vep_full_path = GROUPHOME+MUMMER_PATH+isolate+mummer_extension
        elif isolate+other_extension in os.listdir(GROUPHOME+OTHER_PATH):
            vep_full_path = GROUPHOME+OTHER_PATH+isolate+other_extension
        else:
            print "couldn't find vep file for isolate: " + str(isolate)
            continue

        for line in gzip.open(vep_full_path):
            columns = line.split()
            if "#" in line:
                pass
            else: #read each line that is not header
                AA_isolate_gene=('','','','','')
                bb=()
                isolate_gene=tuple()
                combined_info_lines,combined_info_lines_mut=[],[]
                AA_isolate_gene_list=[]
                del combined_info_lines[:]
                del combined_info_lines_mut[:]
                del AA_isolate_gene_list[:]
                del dup_mut[:]
                position=columns[1]
                ref=columns[3]
                alt=columns[4]
                infos= columns[7]
                CSQs=infos.split('CSQ=')

                if "|" in CSQs[1] and (int(position) in bed_range):

                    CSQ=str(CSQs[1]).split(',')
                    num_of_csqs=len(CSQ)
                    prev_cnq, prev_gene="",""
                    for dd in CSQ:
                        mutation=""
                        #parse data from.high-impact.vcf.gz consequence line
                        info_per_line=dd.split('|')
                        AA = info_per_line[15]
                        gene = info_per_line[4]
                        gene_name=info_per_line[3]
                        if gene_name!="": # if the gene name field is populated, use that, if not, get the gene
                            gene=gene_name
                        Consequence=info_per_line[1]
                        if Consequence == "synonymous_variant":
                            continue
                        IMPACT=info_per_line[2]
                        DISTANCE=info_per_line[18]
                        cDNA_position=info_per_line[12]
                        PICK=info_per_line[19]
                        STRAND=info_per_line[17]
                        VARIANT_CLASS=info_per_line[22]
                        CODON=info_per_line[14]
                        BIOTYPE=info_per_line[7]
                        if "/" in AA: #split mutation to ref and alt
                            AA_split=AA.split('/')
                            AA_ref = AA_split[0]
                            AA_alt = AA_split[1]

                        elif AA != "":
                            AA_ref = AA
                            AA_alt = ""
                        elif AA == "":
                            AA_ref = columns[3]
                            AA_alt = columns[4]

                        # parse mutation according to whether its in a gene or not and if its a snp/indel
                        isolate_gene=(gene,isolate)
                        if "upstream" not in Consequence and "downstream" not in Consequence: # variant in a gene
                            if "SNV" in VARIANT_CLASS: # its a snp
                                if "RNA" in BIOTYPE or "pseudogene" in BIOTYPE:
                                    mutation=AA_ref+cDNA_position+AA_alt
                                    #print mutation
                                    AA_isolate_gene = (mutation, gene, isolate, position, Consequence)
                                    AA_isolate_gene_list.append(AA_isolate_gene)
                                else:
                                    #continue
                                    mutation=AA_ref+CODON+AA_alt
                                    AA_isolate_gene=(mutation,gene,isolate,position,Consequence)
                                    AA_isolate_gene_list.append(AA_isolate_gene)
                            """
                            elif "substitution" in VARIANT_CLASS:
                                if "pseudogene" in BIOTYPE or "RNA" in BIOTYPE:
                                    mutation=AA_ref+cDNA_position+AA_alt
                                    AA_isolate_gene=(mutation,gene,isolate,position,Consequence)
                                    AA_isolate_gene_list.append(AA_isolate_gene)
                                else:
                                    if AA == "" and "upstream" not in prev_cnq: # same variants, different genes
                                        dnap= cDNA_position.split("-")
                                        mutation=str(ref)+str(dnap[-1])+str(alt)
                                        AA_isolate_gene=(mutation+"/"+prev_mutation,gene+"/"+prev_gene,isolate,position,
                                                         Consequence)
                                        AA_isolate_gene_list.append(AA_isolate_gene)
                                        AA_isolate_gene=(mutation,gene,isolate,position,Consequence)
                                        AA_isolate_gene_list.append(AA_isolate_gene)
                                        dup_mut.append((gene,mutation))
                                        dup_mutation.append([x for x in dup_mut])

                                    spl_cod=""
                                    if "-" in CODON:
                                        spl_cod=CODON.split("-")
                                        start_cod=int(spl_cod[0])
                                    else:
                                        start_cod=CODON
                                    i=0
                                    for rf,al in zip(AA_ref,AA_alt):

                                        if rf!="" and al!="":
                                            cod_pos=int(start_cod)+i
                                            i+=1
                                            mutation=str(rf)+str(cod_pos)+str(al)
                                            dup_mut.append((gene,mutation))
                                    mutation=("+".join([x[1] for x in dup_mut]))
                                    AA_isolate_gene=(mutation,gene,isolate,position,Consequence)

                                    AA_isolate_gene_list.append(AA_isolate_gene)
                            elif "insertion" in VARIANT_CLASS or "deletion" in VARIANT_CLASS or "indel"  in VARIANT_CLASS: #its an indel

                                if "-" in cDNA_position:
                                    ll=cDNA_position.split("-") # if it's a range, we only want the start position
                                    cdna_pos=ll[0]
                                else:
                                    cdna_pos=cDNA_position

                                # combined = []
                                # difference = list(Differ().compare(ref, alt))
                                # edited = [base.replace(' ', '') for base in difference[1:]]
                                # for i in edited[1:]:
                                #     if '+' in i:
                                #         combined.append(i.replace('+', ''))
                                #     elif '-' in i:
                                #         combined.append(i.replace('-', ''))
                                #     else:
                                #         combined.append(i)
                                # edited[1:] = combined
                                # reff = ''.join(edited)

                                # for multiple alt's seperated by a comma, use all w/ a "+" b/w each one
                                if "," in alt:
                                    alt_all=[]
                                    sp=alt.split(",")
                                    for S in sp:
                                        alt_all.append(S)
                                    alt=";".join(alt_all)
                                mutation=str(ref)+str(cDNA_position)+str(alt)

                                AA_isolate_gene=(mutation,gene,isolate,position,Consequence)
                                AA_isolate_gene_list.append(AA_isolate_gene)
                            """
                        elif int(DISTANCE) <= 100:# and not args.bed: # not in a gene but within 200 BP of a gene

                        ####### comment this block out and uncomment the line before to go back to original function
                        # else:
                        #     dist_per_csqs=[]
                        #     #append all distances per consq
                        #     for csq_dist in CSQ:
                        #         dist_per_line=csq_dist.split('|')
                        #         dist_per_csqs.append(dist_per_line[18])
                        #     DISTANCE=max(dist_per_csqs)
                        #     # get index of consq w max dist
                        #     max_dist_index=dist_per_csqs.index(DISTANCE)

                        #     for count, csq_dist in enumerate(CSQ):
                        #         if count==max_dist_index:
                        #             dist_per_line=csq_dist.split('|')
                        #             gene = dist_per_line[4]
                        #             gene_name=dist_per_line[3]
                        #             if gene_name!="": # if the gene name field is populated, use that, if not, get the gene
                        #              gene=gene_name

                        ###########

                            mutation=str(ref)+"-"+str(DISTANCE)+str(alt)
                            gene=str(gene)+"_prom"
                            AA_isolate_gene=(mutation,gene,isolate,position,Consequence)
                            AA_isolate_gene_list.append(AA_isolate_gene)
                        prev_cnq=Consequence
                        prev_gene=gene
                        prev_mutation=mutation
                for item in AA_isolate_gene_list:
                    mutation = item[0]
                    gene = item[1]
                    position = item[3]
                    isolate = item[2]
                    consequence_col=item[4]
                    bb = [mutation, gene, consequence_col]
                    if 'R' in s_or_r and bb[0]:
                        dict_isolt_gene_mut_all_R[isolate].append(bb)
                        gene_set_R.add(isolate_gene)
                        AA_set_R.append(','.join([mutation, gene, position]))
                        AA_isolates.append((mutation,gene,position,Consequence,isolate))
                        if ','.join([mutation, gene, position]) in clonality_dictionary.keys():
                            clonality_dictionary[','.join([mutation, gene, position])].append(isolate)
                        else:
                            clonality_dictionary[','.join([mutation, gene, position])] = [isolate]
                    if 'S' in s_or_r and bb[0]:
                        dict_isolt_gene_mut_all_S[isolate].append(bb)
                        gene_set_S.add(isolate_gene)
                        AA_set_S.append(','.join([mutation, gene, position]))
                        AA_isolates.append((mutation,gene,position,Consequence,isolate))
                        if ','.join([mutation, gene, position]) in clonality_dictionary.keys():
                            clonality_dictionary[','.join([mutation, gene, position])].append(isolate)
                        else:
                            clonality_dictionary[','.join([mutation, gene, position])] = [isolate]

        mut_gene_R_all, mut_gene_S_all=[('','')],[('','')]

    AA_mut_list_R=get_count_lists(AA_set_R,"AA")
    AA_mut_list_S=get_count_lists(AA_set_S,"AA")

    output_lines(AA_mut_list_S,AA_mut_list_R,dict_isolt_gene_mut_all_R,dict_isolt_gene_mut_all_S,AA_set_S,AA_set_R,AA_isolates)

def rm_low_dp(dst_file,bed_range, isolate_file):
    MUMMER_PATH = '/data/variants/mummer-denovo/'
    PBHOOVER_PATH = '/data/variants/pbhoover/pbhooverV1.0.0a5/'
    OTHER_PATH = '/data/variants/public-genomes/' # ?
    mummer_extension='.annotated.high-impact.vcf.gz'
    pbhoover_extension= '.vcf.gz'
    other_extension = '.vcf.gz' # ?

    vep_full_path="" # this variable will hold the full path to the vep file

    #Iterate through all isolate.high-impact.vcf.gz files
    for text in isolate_file: # need to add a try block, will fail if file doesn't exist
        isolate=text.strip()
        s_or_r= get_R_or_S(isolate,dst_file)
        if isolate+pbhoover_extension in os.listdir(GROUPHOME+PBHOOVER_PATH):
            vep_full_path = GROUPHOME+PBHOOVER_PATH+isolate+pbhoover_extension
        elif isolate+mummer_extension in os.listdir(GROUPHOME+MUMMER_PATH):
            vep_full_path = GROUPHOME+MUMMER_PATH+isolate+mummer_extension
        elif isolate+other_extension in os.listdir(GROUPHOME+OTHER_PATH):
            vep_full_path = GROUPHOME+OTHER_PATH+isolate+other_extension
        else:
            print "couldn't find vep file for isolate: " + str(isolate)
            continue

        for line in gzip.open(vep_full_path):
            columns = line.split()
            if "#" in line:
                pass
            else: #read each line that is not header
                position=columns[1]
                dps=line.split("DP=")
                DP=dps[1].split(";")
                #l= [''.join(map(str, x)) for x in bed_range if int(position)==x]
                #if len(l)>0:
                #    print l
                pos_strip= int(''.join(map(str,position)))
                if ( int(DP[0]) < 10 ) and pos_strip in bed_range  :
                    bed_range.remove(int(pos_strip))
    return bed_range


#returns which column in the dst file the drug is
def get_drug():
    drugs= header.split()
    col= -1
    for drug in drugs:
        col+=1
        if drug in args.drug:
            return drug, col

# checks the dst file and returns if that isolate is R or S for given drug
def get_R_or_S(isolate,dst_file):
    drug_col=get_drug()
    DRUG=drug_col[0]
    COL=drug_col[1]
    for line in open(GROUPHOME + '/metadata/dst.txt'):
        cols=line.split()
        if isolate == cols[0]:
            return cols[COL]
#get counts for genes and isolates and gene,mutation isolates
def get_count_lists(count_set,count_by):
    gene_mut_list={}
    for mutation in count_set:
        j=0
        if len(mutation)>0:
            if count_by=="gene":
                gene_mut=mutation[0]
            if count_by=="AA":
                gene_mut=mutation
            if gene_mut not in gene_mut_list.keys():
                gene_mut_list[gene_mut]=1
            else:
                gene_mut_list[gene_mut]+=1

    return gene_mut_list

def output_lines(aa_S,aa_R,all_R,all_S,pos_set_S,pos_set_R,aa_isolates):
    gene_mut_count=[("",0)]
    ii_genes=[('')]
    R,S=[],[]
    R_mut,S_mut=[],[]
    output_isolates,output_conseqs=set(),"" # set of isolates/consequences for isolates/consq columns in output
    count_all=[] # gene, mutation, # Rs, # Ss, genomic position, consequence
    count_gene_output=[] # gene, unique # of R per gene, unique # of S per gene, # of isolates per gene that no S has, # of isolates per gene that no R has
                                           # ... , number of unique mutations in the gene,  # of mutations per gene that no S has,  # of mutations per gene that no R has
    gene_mutation_position = {} # this is a list of tuples in the form of (gene,mutation) for all isolates (both R and S)

    for mutation_gene_pos in aa_S:
        gene_mutation_position[mutation_gene_pos] = ''
    for mutation_gene_pos in aa_R:
        gene_mutation_position[mutation_gene_pos] = ''

    mut_w_pos = list(set(pos_set_R) | set(pos_set_S))

    for mutation_gene_pos in gene_mutation_position.keys(): #iterate through all gene and mutation pairs

        mutation = mutation_gene_pos.split(',')[0]
        gene = mutation_gene_pos.split(',')[1]
        position = mutation_gene_pos.split(',')[2]
        if mutation and gene and position:
            z = all_R.copy()
            z.update(all_S)
            output_isolates=set()
            for k,v in z.iteritems():
                for GI in v:
                    if GI[0] == mutation and GI[1]==gene:
                        output_conseqs= GI[2]
                        output_isolates.add(k)

            output_isolates= ','.join(output_isolates)

            if aa_S.has_key(mutation_gene_pos) and aa_R.has_key(mutation_gene_pos):
                count_all.append([gene, mutation, aa_R[mutation_gene_pos], aa_S[mutation_gene_pos], position, output_conseqs,output_isolates])
            elif aa_S.has_key(mutation_gene_pos) and mutation_gene_pos not in aa_R.keys(): # no R's for that gene+mutation
                count_all.append([gene, mutation, 0, aa_S[mutation_gene_pos], position,output_conseqs,output_isolates])
            elif aa_R.has_key(mutation_gene_pos) and mutation_gene_pos not in aa_S.keys(): # no S's for that gene+mutation
                count_all.append([gene, mutation, aa_R[mutation_gene_pos], 0, position, output_conseqs,output_isolates])

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
    f = open('output_4_NN_rif_unexplained.csv', 'w')
    for t in open(args.filename): # need to add a try block, will fail if file doesn't exist
        if "1-0028" in t:
            continue
        isolate=t.strip()
        s_or_r= get_R_or_S(isolate,dst_file)
        output_mut=[0]*len(mut_w_pos)
        output_cols=[0]*(len(mut_w_pos)+1)
        for i,mutt in enumerate(mut_w_pos):
            output_mut[i]=''.join(mutt)

        for gene_col,mut_col,R_count,S_count,pos_col,consq,isolts_col  in count_all:
                for ix,j in enumerate(output_mut):
                    three=str(mut_col)+','+str(gene_col)+','+str(pos_col)
                    if three==j and isolate in isolts_col:
                        output_cols[ix]=1
                    if s_or_r == "R":
                        output_cols[-1]=1
                    elif s_or_r== "S":
                         output_cols[-1]=0
                    else:
                        output_cols[-1]=555
        #print output_cols
        f.write(str(output_cols)+'\n')
    print output_mut

def check_4_dups(a,b,dup_mutation):
    dup=[]
    if len(dup_mutation)>0: # this is for the case where we have a gene inside a gene and we express them as gene1/gene2 but don't want to also print gene1 and gene2
        for d_mut in dup_mutation:
            for j_mut in d_mut:
                if j_mut[0]==a and j_mut[1]==b:
                    dup.append(j_mut)
    return dup

if __name__ == "__main__":

    MAX_DIST=200 # max bp distance from a gene to output
    dst_file=open(GROUPHOME + '/metadata/dst.txt')
    header= dst_file.readline()
    all_islts_dst,dup_mutation,dup_mut=[],[],[] # get the list of all the isolates from the dst file
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Counts mutations/genes for multiple isolates. Input a file with the '
                                                 'list of isolates, a drug name and the what you wish to count by '
                                                 '(mutation/gene/both) and an output file with counts will be generated')
    parser.add_argument("-i", "--input", dest="filename", required=False,action= 'store' ,
                        help="input file with isolate names, one isolate per line", metavar="FILE")
    parser.add_argument("-o", "--output", help="Directs the output to a name of your choice. "
                                               "Will generate an OUTPUT.genes_count or OUTPUT.mutation_count (or both) "
                                               "based on the type of analysis inside a data/GWAS/ dir. Default is to write to "
                                               "standatd output", action='store', dest='output')
    parser.add_argument("-t", "--analysis-type", help="Either 'gene', 'mutation' or 'both'. Default is both ",
                        action='store', dest='type_analysis',required=False, default='both')
    parser.add_argument("-d", "--drug", help="Choose one of the drug names as it appears on the dst table's header ",
                        action='store', dest='drug',required=True,)
    parser.add_argument("-b", "--bed-file", dest="bed", required=False,action= 'store' ,
                        help="input a bed file with 3 columns, tab seperated, CHR START STOP", metavar="FILE")
    parser.add_argument('-c', '--clonality', required=False, action='store_true',
                        help='This prints MIRU and spoligo patterns for the mutation analysis.')

     # i.e the vep file name is 1-0002.vep
    args = parser.parse_args()

    # if the output flag is used, create a data and a GWAS directory to put the data inside
    if args.output:
        try:
            os.makedirs('../data/GWAS/')
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise exc
            pass

    if args.filename:
        isolate_file=open(args.filename)
    else:
        for line in dst_file:
            t=line.split('\t')
            if "1-0028" not in t:
                all_islts_dst.append(t[0])
            isolate_file=all_islts_dst

    bed_range=[]
    # if there's a bed file, get the range of all the positions to check, otherwise use all the positions, need to think of a more efficient way of doing this
    if args.bed:
        bed_file=open(args.bed)
        for POS in bed_file:
            spos = POS.split()
            start=int(spos[1])
            stop=int(spos[2])
            for i in range(start,stop,1):
                bed_range.append(i)
    else:
        for i in range(1,5000000,1):
            bed_range.append(i)

    bed=rm_low_dp(dst_file,bed_range, isolate_file)
    if args.filename:
        isolate_file=open(args.filename)
    else:
        for line in dst_file:
            t=line.split('\t')
            if "1-0028" not in t:
                all_islts_dst.append(t[0])
            isolate_file=all_islts_dst
    count_variation(dst_file,bed,isolate_file)



    """
    # for the mutation argument (or both mutation and gene)
    if (args.type_analysis=="mutation" or args.type_analysis=="both"):
        if not args.clonality:
            if args.output:
                with open('../data/GWAS/'+args.output+".mutation_count", 'w') as output_file:
                    output_file.write("gene" + "\t" + "mutation" + "\t" + "consequence"+ "\t" +"position" +"\t" + "R" + "\t" + "S" + '\n')
                    for gene_col,mut_col,R_count,S_count,pos_col,consq,isolts_col  in count_all:
                        if gene_col!="" and mut_col!="" and len(check_4_dups(gene_col,mut_col,dup_mutation))<1:
                            output_file.write("\t".join((str(gene_col),str(mut_col),str(consq),str(pos_col),str(R_count),str(S_count),str(isolts_col)))+"\n")

            else:  # if output file is not provided, print to standard out
                print "gene" + "\t" + "mutation" + "\t" + "consequence"+ "\t" +"position" +"\t"+ "R" + "\t" + "S"
                for gene_col,mut_col,R_count,S_count,pos_col,consq,isolts_col  in count_all:
                    if gene_col!="" and mut_col!="" and len(check_4_dups(gene_col,mut_col,dup_mutation))<1:  # if dup is not empty, that means this is a gene within a gene and was printed already
                        print "\t".join((str(gene_col),str(mut_col),str(consq),str(pos_col),str(R_count),str(S_count),str(isolts_col)))

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #

    # for the gene argument (or both gene and mutation)
    if (args.type_analysis=="gene" or args.type_analysis=="both"):
            #creating a set of all the genes, to iterate through
            for k,v in all_R.iteritems():
                for i in v:
                    if "/" not in  i[1]:
                        ii_genes.append(i[1])
                    else:
                        [ii_genes.append(x) for x in i[1].split("/")]
            for k,v in all_S.iteritems():
                for i in v:
                    if "/" not in  i[1]:
                        ii_genes.append(i[1])
                    else:
                        [ii_genes.append(x) for x in i[1].split("/")]
            all_uniq_genes=set(ii_genes) # a set of all the unique genes for all the isolates (R and S)

            for i in all_uniq_genes: #iterate through all the genes and get count for each gene
                if i!="":
                    num_of_isolates_R,num_of_isolates_S,num_of_mutations,num_of_mutations_R,num_of_mutations_S,isolt_R_all_count,isolt_R_all_count=0,0,0,0,0,0,0
                    del R[:]
                    del S[:]
                    num_of_mutations=[x[0] for x in count_all].count(i) #count how many unique mutations total
                    num_of_mutations_R= ([x[0] for x in count_all if x[3]==0]).count(i) #if number of S's is 0, count how many mutations in that gene
                    num_of_mutations_S= ([x[0] for x in count_all if x[2]==0]).count(i) # ^ same but for S.

                    isolates_R_in_gene = set([pairs[0] for k,v in all_R.iteritems() for pairs in v if pairs[1]==str(i)]) # all S mutations for that gene
                    isolates_S_in_gene = set([pairs[0] for k,v in all_S.iteritems() for pairs in v if pairs[1]==str(i)])

                    R,S=[],[]
                    for k,v in all_R.iteritems(): # k is the isolate and v is a pair (mutation, gene)
                        for pairs in v:
                            if pairs[1]==i: # for the gene of interest
                                R.append((pairs[0],k)) # append mutation and isolate
                            elif "/" in pairs[1] and i in pairs[1]:
                                R.append((pairs[0],k))
                    for k,v in all_S.iteritems():
                        for pairs in v:
                            if pairs[1]==i:
                                S.append((pairs[0],k))
                            elif "/" in pairs[1] and i in pairs[1]:
                                R.append((pairs[0],k))

                    isolt_R_all_count=len(set(([r[1] for r in R]))) #count how many unique isolates for that gene
                    isolt_S_all_count=len(set(([s[1] for s in S])))

                    isolate_list_col=set(([r[1] for r in R])) | set(([s[1] for s in S])) # set of all the isolates that has (gene,mutation)
                    isolate_list_col= ','.join(isolate_list_col)

                    bloo=set([k for k,v in all_R.iteritems() for pairs in v if pairs[1]==str(i)])
                    Sloo=set([k for k,v in all_S.iteritems() for pairs in v if pairs[1]==str(i)])

                    num_of_isolates_R=len(set([vloo for vloo in bloo for r in R if (str(r[1])== str(vloo) and r[0] not in isolates_S_in_gene)]))
                    num_of_isolates_S=len(set([svloo for svloo in Sloo for s in S if (str(s[1])== str(svloo) and s[0] not in isolates_R_in_gene)]))

                    count_gene_output.append((i,isolt_R_all_count,isolt_S_all_count,num_of_isolates_R,num_of_isolates_S,num_of_mutations,num_of_mutations_R,num_of_mutations_S,isolate_list_col))
            if args.output:
                with open('../data/GWAS/'+args.output+".genes_count", 'w') as output_file:
                    output_file.write("Gene" + "\t" + " Total Number R (unique) [isolate count]" + "\t" + " Total Number S (unique) [isolate count]"+ "\t" + "# R with mutations that No S has [isolate count, subset of Col2]" + "\t" + "# S with mutations that No R has [isolate count, subset of Col3]" + "\t" + "# of unique mutations" + "\t" + " # Mutations that only R has [mutation count]"+ "\t" + " # Mutations that only S has [mutation count]" + '\n')
                    for a,b,c,d,e,f,g,h,isolts_col in count_gene_output: #unpack and print all the counts
                        if a!="" and b!="" and c!="" and e!="" and f!="" and g!="" and h!="":
                            output_file.write("\t".join((str(a),str(b),str(c),str(d),str(e),str(f),str(g),str(h),str(isolts_col)))+"\n")
            else: # if output file is not provided, print to standard out
                print "Gene" + "\t" + " Total Number R (unique) [isolate count]" + "\t" + " Total Number S (unique) [isolate count]"+ "\t" + "# R with mutations that No S has [isolate count, subset of Col2]" + "\t" + "# S with mutations that No R has [isolate count, subset of Col3]" + "\t" + "# of unique mutations" + "\t" + " # Mutations that only R has [mutation count]"+ "\t" + " # Mutations that only S has [mutation count]" + '\n'
                for a,b,c,d,e,f,g,h,isolts_col in count_gene_output: #unpack and print all the counts
                        if a!="" and b!="" and c!="" and e!="" and f!="" and g!="" and h!="":
                            print "\t".join((str(a),str(b),str(c),str(d),str(e),str(f),str(g),str(h),str(isolts_col)))