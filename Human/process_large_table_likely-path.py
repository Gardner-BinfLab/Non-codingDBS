#! /usr/bin/env python3.7
import sys
import pandas as pd
import math

csv_file_name=sys.argv[1] # a file prepared by Mini and converted to csv
likely_path_bed_file_name=sys.argv[2]
likely_path_out_file_name=sys.argv[3] # a file to contain chr, region, change, strand info for

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

likely_path_df=pd.DataFrame([], columns=['Chr', 'Location', 'Reference', 'Change', 'Change type'])

for index, row in df.iterrows():
    if str(row["GRCh38Location"])!="nan":
        chr=str(row["GRCh38Chromosome"])
        location=str(row["GRCh38Location"])
        reference=str(row["Ref allele"])
        change=str(row["Alt allele"])
        type=str(row["Gardner Group Classification"])
        if type.lower()=="likely pathogenic" and chr!="nan" and location!="nan" and change!="nan" and reference!="nan":
            if reference in bases and change in bases:
                change_type="SNP"
            elif reference=="-":
                change_type="insertion"
            elif change=="-":
                change_type="deletion"
            else:
                change_type="deletion-insertion" #accounts for duplications
            new_row=pd.DataFrame([[chr, location, reference, change, change_type]], columns=['Chr', 'Location', 'Reference', 'Change', 'Change type'])
            likely_path_df=likely_path_df.append(new_row)

likely_path_intervals_by_chr={}

likely_path_out_file = open(likely_path_out_file_name, 'w')

for index, row in likely_path_df.iterrows():
    change_type=str(row["Change type"])
    locations=row["Location"].split("-")
    left_right=get_left_right(locations, change_type)
    left=left_right[0]
    right=left_right[1]
    chr_name="Chr"+str(row["Chr"])
    likely_path_out_file.write(chr_name + " " + str(left) + " " + str(right) + " " + str(row["Reference"]) + " " + str(row["Change"]) + " " + change_type + "\n")
    if chr_name in likely_path_intervals_by_chr:
        likely_path_intervals_by_chr[chr_name].append([left, right])
    else:
        likely_path_intervals_by_chr[chr_name]=[[left, right]]

likely_path_out_file.close()

likely_path_intervals_by_chr_merged=merge_all_intervals(likely_path_intervals_by_chr)

likely_path_bed_file = open(likely_path_bed_file_name, 'w')

for key, intervals in likely_path_intervals_by_chr_merged.items():
    for interval in intervals:
        likely_path_bed_file.write(key + " " + str(interval[0]) + " " + str(interval[1]) + "\n")

likely_path_bed_file.close()
