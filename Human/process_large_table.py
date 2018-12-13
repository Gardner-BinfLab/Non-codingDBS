#! /usr/bin/env python3.7
import sys
import pandas as pd
import math

csv_file_name=sys.argv[1] # a file prepared by Mini and converted to csv
pathogenic_bed_file_name=sys.argv[2]
benign_bed_file_name=sys.argv[3]
pathogenic_out_file_name=sys.argv[4] # a file to contain chr, region, change, strand info for
benign_out_file_name=sys.argv[5]

def merge_intervals(intervals):
    new_intervals=[]
    min=math.inf
    max=0
    for interval in intervals:
        if interval[0] < min:
            min=interval[0]
        if interval[1] > max:
            max=interval[1]
    bits=[ 0 for x in range(max-min+2) ] # the length of the bits is one bit longer to mark the end of the last interval when converting back
    for interval in intervals: # convert intervals represented by ranges to intervals represented by bits (shifted by min) and simultaneously intersect them
        for i in range(interval[0], interval[1]+1):
            bits[i-min]=1
    reading=False
    for i in range(len(bits)): # find the resulting overlapping intervals and convert to ranges
        if bits[i]==1 and not reading:
            reading=True
            new_interval=[i+min]
        elif bits[i]==0 and reading:
            reading=False
            new_interval.append(i-1+min)
            new_intervals.append(new_interval)
    return new_intervals

def merge_all_intervals(all_intervals):
    new_all_intervals={}
    for key in all_intervals:
        new_all_intervals[key]=merge_intervals(all_intervals[key])
    return new_all_intervals

def get_left_right(locations, change_type):
    L1=int(locations[0])
    if len(locations)>1:
        L2=int(locations[1])
    else:
        L2=L1
    if change_type=="SNP" or change_type=="deletion-insertion" or change_type=="deletion":
        left=L1-101
        right=L2+100
    elif change_type=="insertion":
        left=L1-100
        right=L2+99
    else:
        print("Invalid change type: " + change_type )
    return [left, right]


bases=["A", "G", "C", "T"]

df=pd.read_csv(csv_file_name, sep=",")

pathogenic_df=pd.DataFrame([], columns=['Chr', 'Location', 'Reference', 'Change', 'Change type'])
benign_df=pd.DataFrame([], columns=['Chr', 'Location', 'Reference', 'Change', 'Change type'])

for index, row in df.iterrows():
    if str(row["GRCh38Location"])!="nan":
        chr=str(row["GRCh38Chromosome"])
        location=str(row["GRCh38Location"])
        reference=str(row["Ref allele"])
        change=str(row["Alt allele"])
        type=str(row["Gardner Group Classification"])
        if chr!="nan" and location!="nan" and change!="nan" and reference!="nan":
            if reference in bases and change in bases:
                change_type="SNP"
            elif reference=="-":
                change_type="insertion"
            elif change=="-":
                change_type="deletion"
            else:
                change_type="deletion-insertion" #accounts for duplications
            new_row=pd.DataFrame([[chr, location, reference, change, change_type]], columns=['Chr', 'Location', 'Reference', 'Change', 'Change type'])
            if type.lower()=="pathogenic":
                pathogenic_df=pathogenic_df.append(new_row)
            elif type.lower()=="benign":
                benign_df=benign_df.append(new_row)

pathogenic_intervals_by_chr={}

pathogenic_out_file = open(pathogenic_out_file_name, 'w')

for index, row in pathogenic_df.iterrows():
    change_type=str(row["Change type"])
    locations=row["Location"].split("-")
    left_right=get_left_right(locations, change_type)
    left=left_right[0]
    right=left_right[1]
    chr_name="Chr"+str(row["Chr"])
    pathogenic_out_file.write(chr_name + " " + str(left) + " " + str(right) + " " + str(row["Reference"]) + " " + str(row["Change"]) + " " + change_type + "\n")
    if chr_name in pathogenic_intervals_by_chr:
        pathogenic_intervals_by_chr[chr_name].append([left, right])
    else:
        pathogenic_intervals_by_chr[chr_name]=[[left, right]]

pathogenic_out_file.close()

pathogenic_intervals_by_chr_merged=merge_all_intervals(pathogenic_intervals_by_chr)

pathogenic_bed_file = open(pathogenic_bed_file_name, 'w')

for key, intervals in pathogenic_intervals_by_chr_merged.items():
    for interval in intervals:
        pathogenic_bed_file.write(key + " " + str(interval[0]) + " " + str(interval[1]) + "\n")

pathogenic_bed_file.close()

benign_intervals_by_chr={}

benign_out_file = open(benign_out_file_name, 'w')

for index, row in benign_df.iterrows():
    change_type=str(row["Change type"])
    locations=row["Location"].split("-")
    left_right=get_left_right(locations, change_type)
    left=left_right[0]
    right=left_right[1]
    chr_name="Chr"+str(row["Chr"])
    benign_out_file.write(chr_name + " " + str(left) + " " + str(right) + " " + str(row["Reference"]) + " " + str(row["Change"]) + " " + change_type + "\n")
    if chr_name in benign_intervals_by_chr:
        benign_intervals_by_chr[chr_name].append([left, right])
    else:
        benign_intervals_by_chr[chr_name]=[[left, right]]

benign_out_file.close()

benign_intervals_by_chr_merged=merge_all_intervals(benign_intervals_by_chr)

benign_bed_file = open(benign_bed_file_name, 'w')

for key, intervals in benign_intervals_by_chr_merged.items():
    for interval in intervals:
        benign_bed_file.write(key + " " + str(interval[0]) + " " + str(interval[1]) + "\n")

benign_bed_file.close()

