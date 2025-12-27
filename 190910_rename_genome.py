# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 17:41:21 2019

@author: dajun
"""
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(prog=None, description='change genome scaffold name and sorted')
parser.add_argument('-i', metavar='if', type=str, help='input genome')
parser.add_argument('-m', help='input genome name')
parser.add_argument('-o', metavar='of', type=str, help='output_file')
args = parser.parse_args()

input_file = args.i
genome_name = args.m
output_file = args.o


#Get the lengths and ids, and sort on length
len_and_ids = sorted((len(rec), rec.id) for rec in SeqIO.parse(input_file,"fasta"))
#sorted and ouput
ids = reversed([id_ for (length, id_) in len_and_ids]) 
record_index = SeqIO.index(input_file, "fasta")
records = (record_index[id_] for id_ in ids)
my_seq=[]
a= genome_name
sca=1
for sor in records:
    sor.id = a + "_scaffold" + str(sca) + "|length_" + str(len(sor.seq))
    sca=sca+1
    sor.description = ""
    my_seq.append(sor)
SeqIO.write(my_seq, output_file, "fasta")