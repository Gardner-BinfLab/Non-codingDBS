#! /usr/bin/env python3.7
import sys
import pandas as pd
import math

csv_file_name=sys.argv[1] # a file with pathogenic or benign variants containing fields: "ID", "chr", "pos", "ref", "alt"
bed_folder_name=sys.argv[2]
out_file_name=sys.argv[3] # a file to contain chr, region, change, strand info

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
    if change_type=="SNP" or change_type=="deletion":
        left=L1-101
        right=L2+100
    elif change_type=="insertion":
        left=L1-100
        right=L2+100
    else:
        print("Invalid change type: " + change_type )
    return [left, right]


bases=["A", "G", "C", "T"]

df=pd.read_csv(csv_file_name, sep=",", encoding = "ISO-8859-1")

variants_df=pd.DataFrame([], columns=['Chr', 'Location', 'Reference', 'Change', 'Change type'])

for index, row in df.iterrows():
    if str(row["ID"])!="nan":
        chr=str(row["chr"])
        location=str(row["pos"])
        reference=str(row["ref"])
        if len(reference)>1:
            end=int(location)+len(reference)-1
            location=location+"-"+str(end)
        change=str(row["alt"])
        if chr!="nan" and location!="nan" and change!="nan" and reference!="nan":
            if reference!="-":
                for nucleotide in reference:
                    if nucleotide.upper() not in bases:
                        print("ERROR: reference sequence at " + location + " is not from the correct alphabet: " + reference)
            if change != "-":
                for nucleotide in change:
                    if nucleotide.upper() not in bases:
                        print("ERROR: alternative sequence at " + location + " is not from the correct alphabet: " + change)
            if reference[0] in bases and change[0] in bases:
                change_type="SNP"
            elif reference=="-":
                change_type="insertion"
            elif change=="-":
                change_type="deletion"
            new_row=pd.DataFrame([[chr, location, reference, change, change_type]], columns=['Chr', 'Location', 'Reference', 'Change', 'Change type'])
            variants_df=variants_df.append(new_row)

intervals_by_chr={}

out_file = open(out_file_name, 'w')

for index, row in variants_df.iterrows():
    change_type=str(row["Change type"])
    locations=row["Location"].split("-")
    left_right=get_left_right(locations, change_type)
    left=left_right[0]
    right=left_right[1]
    chr_name="Chr"+str(row["Chr"])
    out_file.write(chr_name + " " + str(left) + " " + str(right) + " " + str(row["Reference"]) + " " + str(row["Change"]) + " " + change_type + "\n")
    if chr_name in intervals_by_chr:
        intervals_by_chr[chr_name].append([left, right])
    else:
        intervals_by_chr[chr_name]=[[left, right]]

out_file.close()

intervals_by_chr_merged=merge_all_intervals(intervals_by_chr)

i=1
bed_file_name = bed_folder_name+"/bed"+str(i)
bed_file = open(bed_file_name, 'w')
count=0
for key, intervals in intervals_by_chr_merged.items():
    for interval in intervals:
        if count==1000:
            bed_file.close()
            i=i+1
            bed_file_name = bed_folder_name+"/bed"+str(i)
            bed_file = open(bed_file_name, 'w')
            count=0
        bed_file.write(key + " " + str(interval[0]) + " " + str(interval[1]) + "\n")
        count=count+1
bed_file.close()

