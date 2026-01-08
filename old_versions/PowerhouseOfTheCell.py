#!/usr/bin/env python3
class FastaRecord():

    def __init__(self,header:str,array_of_sequences:str):

        self.header = header.rstrip('\n')

        # the sequence is currently a list of sequences so we need to concatenate them together
        self.sequence = ''.join([Line.rstrip() for Line in array_of_sequences])

        # get the sequence length
        self.sequence_length = len(self.sequence)

        # when outputting the sequence we want to preserve the original structure of the input file including new lines
        self.output_sequence = ''.join(array_of_sequences)
        self.output_header = header.rstrip('\n') + '\n'

        # get the number of each base
        self.A_count,self.C_count,self.T_count,self.G_count = find_ACTG_content(sequence=self.sequence)
        self.A_percentage = (self.A_count/self.sequence_length) * 100
        self.C_percentage = (self.C_count/self.sequence_length) * 100
        self.T_percentage = (self.T_count/self.sequence_length) * 100
        self.G_percentage = (self.G_count/self.sequence_length) * 100


def create_multiline_fasta(sequence:str,characters_per_line:int=66) -> list:

    output = []

    for i in range(0,len(sequence),characters_per_line):

        output.append(f"{sequence[i:i+characters_per_line]}\n")
    
    return output

def subset_fasta(sequence:str,start:int,end:int) -> str:

    return sequence[start:end]
    
# useful functions
def find_ACTG_content(sequence:str) -> tuple:

    return sequence.count('A'),sequence.count('C'),sequence.count('T'),sequence.count('G')

# create dict with header name as key and Fasta Record object as value
def fasta_to_dict(FastaData:list) -> dict:

    # scaffold and their sequence will be stored here
    RefDict = {}

    # current seq is a list because we will be storing sequences without removing the new line, so we want each header to have a 'group' of sequences
    CurrentSeq = []
    CurrentHeader = FastaData[0]

    # generate dictionary with gene names and fasta files
    for i in range(1,len(FastaData)):

        # find what the next line is
        Line = FastaData[i]

        # find if the line begins with a gene
        if Line.startswith('>'): 

            # if line begins with a gene and there is some sequence for the previous gene that was found, add it to the dictionary reset the sequence builder
            if CurrentSeq:
                RefDict[CurrentHeader.replace('>','').rstrip()]=FastaRecord(header=CurrentHeader.replace('>','').rstrip(),array_of_sequences=CurrentSeq)
                CurrentHeader=Line

            CurrentSeq = []

        # if a line starts with the last line/sequence in the file, add it to the current sequence builder and then add the current gene to the ref dictionary
        elif i == len(FastaData) - 1:
            CurrentSeq.append(Line)
            RefDict[CurrentHeader.replace('>','').rstrip()]=FastaRecord(header=CurrentHeader.replace('>','').rstrip(),array_of_sequences=CurrentSeq)

        # if none of the above conditions are met then keep building the sequence
        else:
            CurrentSeq.append(Line)

    return RefDict

# create dict with scaffold name as key
def fasta_to_dict_origina(FastaData:list):

    # scaffold and their sequence will be stored here
    RefDict = {}

    # current seq is a list because we will be storing sequences without removing the new line, so we want each header to have a 'group' of sequences
    CurrentSeq = []
    CurrentHeader = FastaData[0]

    # generate dictionary with gene names and fasta files
    for i in range(1,len(FastaData)):

        # find what the next line is
        Line = FastaData[i]

        # find if the line begins with a gene
        if Line.startswith('>'): 

            # if line begins with a gene and there is some sequence for the previous gene that was found, add it to the dictionary reset the sequence builder
            if CurrentSeq:
                RefDict[CurrentHeader.replace('>','')]=CurrentSeq
                CurrentHeader=Line

            CurrentSeq = []

        # if a line starts with the last line/sequence in the file, add it to the current sequence builder and then add the current gene to the ref dictionary
        elif i == len(FastaData) - 1:
            CurrentSeq.append(Line)
            RefDict[CurrentHeader.replace('>','')]=CurrentSeq

        # if none of the above conditions are met then keep building the sequence
        else:
            CurrentSeq.append(Line)

    return RefDict

class GFFfile():

    def __init__(self):
        pass


