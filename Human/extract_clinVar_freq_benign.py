#! /usr/bin/env python3.7
import sys
import pandas as pd
import math

clinVar_file_name=sys.argv[1]
UCSC_file_name=sys.argv[2]
benign_out_file_name=sys.argv[3]

def get_left_right(locations, change_type, ref_len):
    L1=int(locations[0])
    if len(locations)>1:
        L2=int(locations[1])
    elif len(locations)==1 and ref_len>1:
        L2=L1
        L1=L2-ref_len+1
    else:
        L2=L1
    if change_type=="SNP" or change_type=="deletion-insertion" or change_type=="deletion":
        left=L1-101
        right=L2+100
    elif change_type=="insertion":
        left=L1-100
        right=L2+99
    else:
        print("Invalid change type: " + change_type)
    return [left, right]

bases=["A", "G", "C", "T"]

df_clinVar=pd.read_csv(clinVar_file_name, sep=",")

benign_clinVar_df=pd.DataFrame([], columns=['Chr', 'Location', 'Reference', 'Change', 'Change type'])

for index, row in df_clinVar.iterrows():
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
            if type.lower()=="benign":
                benign_clinVar_df=benign_clinVar_df.append(new_row)

df_UCSC=pd.read_csv(UCSC_file_name, sep=",")

benign_UCSC_df=pd.DataFrame([], columns=['Chr', 'Location', 'Reference', 'Change', 'Change type', 'Frequencies'])

for index, row in df_UCSC.iterrows():
    if str(row["GRCh38Location"])!="nan":
        chr=str(row["GRCh38Chromosome"])
        location=str(row["GRCh38Location"])
        reference=str(row["Ref allele"])
        change=str(row["Alt allele"])
        freqs=str(row["allele freqs"])
        if chr!="nan" and location!="nan" and change!="nan" and reference!="nan":
            if reference in bases and change in bases:
                change_type="SNP"
            elif reference=="-":
                change_type="insertion"
            elif change=="-":
                change_type="deletion"
            else:
                change_type="deletion-insertion" #accounts for duplications
            new_row=pd.DataFrame([[chr, location, reference, change, change_type, freqs]], columns=['Chr', 'Location', 'Reference', 'Change', 'Change type', 'Frequencies'])
            benign_UCSC_df=benign_UCSC_df.append(new_row)

benign_out_file = open(benign_out_file_name, 'w')

for index, row in benign_clinVar_df.iterrows():
    change_type=str(row["Change type"])
    locations=row["Location"].split("-")
    left_right=get_left_right(locations, change_type, len(str(row["Reference"])))
    left=left_right[0]
    right=left_right[1]
    chr_name="Chr"+str(row["Chr"])
    benign_out_file.write(chr_name + "_" + str(left) + "-" + str(right) + "\tclinVar\n")

for index, row in benign_UCSC_df.iterrows():
    change_type=str(row["Change type"])
    locations=row["Location"].split("-")
    left_right=get_left_right(locations, change_type, len(str(row["Reference"])))
    left=left_right[0]
    right=left_right[1]
    chr_name="Chr"+str(row["Chr"]).replace("chr", "")
    freqs=str(row["Frequencies"]).split(",")
    freqs_set=set([x.strip() for x in freqs])
    minor_freq=""
    if freqs_set=={"10","90"}:
        minor_freq="10"
    elif freqs_set=={"20", "80"}:
        minor_freq="20"
    elif freqs_set=={"50"}:
        minor_freq="50"
    else:
        print("Frequencies format is non-standard: " + row["Frequencies"])
    if minor_freq!="":
        benign_out_file.write(chr_name + "_" + str(left) + "-" + str(right) + "\t" + minor_freq+ "\n")

benign_out_file.close()
