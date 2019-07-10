#!/usr/bin/env python
#import libraries
import argparse
import re

#Set input variables
def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--file", help="File name", type=str)
    return parser.parse_args()

args = get_args()

#Set global variables
FILE = args.file

#Array to store values
header_array = [] #information of the header array
k_cov_mean = [] #K-mer coverage mean array per contig
N50_array = [] #stores every length of a contig to find the N50
contig_dict = dict() # dictionary used to calculate distrubution of contigs KEY = contig length (veries by 100) VALUE = number of contigs
total_length = 0
running_length = 0
num_contigs = 0 #used to keep track of the total number of records
max_contig_len = 0 #maximum length of a contig
record_len = 0 #variable used to store each record
mean_contig = 0 #mean contig value
Minimum_contig = 1000000 #Minimum contig to check distribution

with open(FILE, "r") as fh:
    LN = 0 #Line number
    for line in fh:
        LN+=1
        line = line.strip('\n') #Manipulate the sequence line
        if line.startswith(">"):
            x = re.search("^>[A-Za-z]+_[0-9]+_[A-Za-z]+_([0-9]+)_[A-Za-z]+_(.+)", line)
            kcount = float(x.group(1))
            seq_length = kcount - 1 + 49 #calculation used to get each sequence's length
            running_length += seq_length

            coverage = (float(x.group(2)) * seq_length) / (seq_length - 49 + 1) #calculation of coverage per contigs
            k_cov_mean.append(coverage)

            if record_len != 0: #stores each record length into the array
                N50_array.append(record_len)

            if record_len < Minimum_contig and record_len != 0: #stores the shortest contig length
                Minimum_contig = record_len

            num_contigs += 1
            record_len = 0 #reset the contig length to 0 to keep track of a different contig

        else: #If the line that was read in isnt a header
            record_len += len(line)
            total_length += len(line)
            if record_len > max_contig_len:
                max_contig_len = record_len

#N50 Calculation
N50_array.sort(reverse = True)
N50_val = 0 # The N50 value
length_sum = 0  #stores the length when adding the record's lengths
index = 0 #keeps track of the N50_array's index
N50_sorting_val = total_length // 2 #The value to compare the length_sum to
boolN50 = False
while boolN50 == False:
    length_sum += N50_array[index]
    if length_sum > N50_sorting_val:
        N50_val = N50_array[index]
        boolN50 = True
    index+=1

index = 0
# Find the distribution of the contigs KEY = contig length (veries by 100) VALUE = number of contigs
for element in N50_array:
    Key = N50_array[index] // 100 #sort into a dictionary based on the 100s digit
    Key *= 100  #multiply by 100 to make the key be in the 100s digit
    if Key in contig_dict.keys():
        contig_dict[Key] += 1
    else:
        contig_dict.setdefault(Key,1)
    index += 1

print("# Contig length\t\tNumber of contigs in this category")
for i in sorted (contig_dict.keys()) :
     print("%s\t%s" % (i, contig_dict[i]))

#Calculation of mean contig length
mean_contig = total_length/num_contigs

#Calculation of K-mer coverage mean
mean_depth = sum(k_cov_mean) / len(k_cov_mean)

#print(header_array)
print("running_length", running_length)
print("total_length", total_length)
print ("number of contigs", num_contigs)
print ("maximum contig length", max_contig_len)
print ("mean contig length",mean_contig)
print ("N50", N50_val)
print("minimum contig length", Minimum_contig)
print ("mean depth of coverage", mean_depth)
