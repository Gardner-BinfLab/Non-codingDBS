#! /usr/bin/env python3.5
import sys

def extract_score(file_name):
    file=open(file_name)
    ind=0
    for i, line in enumerate(file):
        if ind > 0 and i==ind:
            x=line.split()
            score=x[1]
            break
        elif "Scores for complete hits:" in line:
            ind=i+3
    file.close()
    return score

input_file=sys.argv[1]

seq_name=input_file.split("/")[-1].split(".")[0]

print(seq_name+"\t"+extract_score(input_file))
