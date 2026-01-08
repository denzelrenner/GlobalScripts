#!/usr/bin/env python3
# after diamond runs are done, we want to verify that each of the large scaffolds got a hit
# to do that first we need to find which scaffolds are large scaffolds
# then we look in the diamond outdir to see if the file is empty
# if it is, we resample from the scaffold again
import os,argparse,time,json
from PowerhouseOfTheCell import *

vars = argparse.ArgumentParser(description='Filter out contaminated scaffolds')

vars.add_argument('-f','--fa',dest='fasta', type=str, required=True, help='assembly fasat')

vars.add_argument('-b','--blobdir',dest='blobdir', type=str, required=True, help='path to blobdir')

vars.add_argument('-q','--query_taxa',dest='query_taxa', type=str,required=False,default='phylum', help='which taxa rank(s) json we want to use for filtering the fasta file')

vars.add_argument('-v','--query_value',dest='query_value', type=str,nargs='+',required=False,default=['no-hit','Streptophyta'], help='given the taxonomic rank(s), what specifc i.e family are we interested in')

vars.add_argument('-od','--outdir',dest='outdir', type=str, required=False,default=None, help='name of directory for output. will be stored in current working directory if none given')

args = vars.parse_args()

# by default directory to write output to is current wd
outdir = os.getcwd()

# if an actual output directory is given
if args.outdir:

    outdir = args.outdir # set output dir

    # check if dir given exists, if not create one
    if not os.path.isdir(f"{outdir}"):
        os.makedirs(f"{outdir}")

# store all files from blobdir
files = [f"{args.blobdir}/{f}" for f in os.listdir(args.blobdir)]

# set vars for files we will use
identifiers = files[files.index(f"{args.blobdir}/identifiers.json")]
taxa = files[files.index(f"{args.blobdir}/bestsumorder_{args.query_taxa}.json")]


# load files
with open(identifiers,'r') as identifier_file, open(taxa,'r') as taxa_file, open(args.fasta,'r') as FastaFile:

    # store the fasta data for the assembly in a list but do not remove new lines
    FastaData = FastaFile.readlines()

    # chromosomes and their sequence will be stored here
    FastaRefDict = fasta_to_dict(FastaData)

    # load jsons 
    identifier_data = json.load(identifier_file)
    taxa_data = json.load(taxa_file)

    # get list of scaffolds as they appear in blobdir json
    identifier_scaffolds = identifier_data['values']

    # index the taxonomic assignments for each scaffold
    indexed_taxonomic_assignments = index_blobdir_json(taxa_data)

    # index the ordered list of possible taxonomic assignments, different from above which is the actual assignment to scaffolds, this is just they key
    indexed_taxa_keys = index_blobdir_json(taxa_data,key='keys')

    # each scaffold is assigned the index value of one of the keys, but now we want to convert it to the actual taxa name and store it in a dict with its scaffold
    scaffold_taxa = {}

    for i in range(len(identifier_scaffolds)):

        scaf = identifier_scaffolds[i]

        # this looks in the dict of taxa keys, and given the scaffold we are currently on, takes the value in its indexed assignment
        taxonomic_assignment = indexed_taxa_keys[indexed_taxonomic_assignments[i]]

        scaffold_taxa[scaf] = taxonomic_assignment

    # now create a list of scaffolds with and without the query taxa
    with_query_taxa = [k for k,v in scaffold_taxa.items() if v in args.query_value]
    without_query_taxa = [k for k,v in scaffold_taxa.items() if v not in args.query_value]

    # subset fasta file to only contain or not contain scaffolds with query taxa
    fasta_query_taxa = subset_whole_fasta(FastaRefDict=FastaRefDict,records=with_query_taxa)
    fasta_wo_query_taxa = subset_whole_fasta(FastaRefDict=FastaRefDict,records=without_query_taxa)

    # output filtered fasta files with only query taxa
    output_fasta_file(FastaRefDict=fasta_query_taxa ,output_directory=outdir,file_name=f'{args.fasta.split("/")[-1].rstrip(".fa")}.wquery.fa')
    output_fasta_file(FastaRefDict=fasta_wo_query_taxa ,output_directory=outdir,file_name=f'{args.fasta.split("/")[-1].rstrip(".fa")}.woquery.fa')

