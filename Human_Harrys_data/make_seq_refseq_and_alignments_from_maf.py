#! /usr/bin/env python3.7
import sys
import pandas as pd
from Bio import AlignIO
from Bio.AlignIO import MafIO
from Bio import SeqIO
from Bio.Seq import Seq
import copy

def find_chr_name(id, allnames):
    valid_name=""
    for name in allnames:
        if name in id and valid_name in name:
                valid_name=name
    return valid_name

def get_strand(change_type, ref_sequence, left, right, alignment):
    complements={"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    strand=0
    ref_region=str(ref_sequence).upper()
    column_at_position_of_change=alignment.get_spliced([int(left)+100],[int(right)-100], strand=1)
    for rec in column_at_position_of_change:
        if rec.id==chr_name:
            region_from_maf=str(rec.seq.ungap("-")).upper()
    complement=""
    for nucleotide in region_from_maf:
        complement=complements[nucleotide]+complement
    if ref_region==region_from_maf:
        strand=1
    elif ref_region==complement:
        strand=-1
    else:
        print("Reference sequence in table:\t" + ref_region+ "\nSequence in maf alignment:\t" + region_from_maf + "\nComplemented sequence in maf alignment:\t" + complement)
    return strand

def get_alternative_and_reference_sequences(change, alignment, chr_name):
    for record in alignment:
        if chr_name==record.id:
            ref_sequence=record
    alt_sequence=copy.deepcopy(ref_sequence)
    ungaped_seq=str(ref_sequence.seq.ungap("-"))
    length=len(ungaped_seq)
    if str(change)!="-":
        insert_region=str(change)
    else:
        insert_region=""
    ref_sequence.seq=Seq(ungaped_seq)
    alt_seq=ungaped_seq[0:100]+insert_region+ungaped_seq[length-100:length]
    alt_sequence.seq=Seq(alt_seq)
    return [alt_sequence, ref_sequence]

input_file_name=sys.argv[1] #csv file with regions and ref/alt data
maf_file_name=sys.argv[2]
chr_mafs_dir=sys.argv[3]
seq_folder=sys.argv[4]
alignment_folder=sys.argv[5]
refseq_folder=sys.argv[6]

type=input_file_name.split(".")[0]

df=pd.read_csv(input_file_name, sep="\s+", names=["Chr", "Left", "Right", "Reference", "Change", "Change type"])

ma=AlignIO.parse(maf_file_name, "maf")
list_ma=list(ma)
hgnames=set()

for ind_ma in list_ma:
    for record in ind_ma:
        if "hg38" in record.id:
            hgnames.add(record.id)

ma_by_chr={}
for name in hgnames:
    ma_by_chr[name]=list()

for ind_ma in list_ma:
    for record in ind_ma:
        if "hg38" in record.id:
            name=find_chr_name(record.id, hgnames)
            ma_by_chr[name].append(ind_ma)
            break

index_by_chr={}

base_maf_file_name=maf_file_name.split("/")[-1].split(".")[0]

for name in hgnames:
    chr_name=name.split(".")[1]
    chr_maf_path=chr_mafs_dir.strip("/")+"/"+base_maf_file_name+"_"+chr_name
    AlignIO.write(ma_by_chr[name], chr_maf_path+".maf", "maf")
    index_by_chr[name]=MafIO.MafIndex(chr_maf_path+".mafindex", chr_maf_path+".maf", name)
version_count_by_seq_name={}

df.drop_duplicates(inplace=True)

for index, row in df.iterrows():
    left=row["Left"]
    right=row["Right"]
    change_type=row["Change type"]
    seq_name=str(row["Chr"])+ "_" + str(left) + "-" + str(right)
    chr_name="hg38.chr"+str(row["Chr"]).lstrip("Chr")
    if chr_name not in hgnames:
        print("ERROR: the naming of human genomes in alignment " + maf_file_name + " does not comply with the expected standard: hg38.chrNum: "+chr_name + " found.")
        sys.exit(0)
    if change_type=="insertion":
        strand=1
    else:
        strand=get_strand(change_type, row["Reference"], left, right, index_by_chr[chr_name])
    if strand==0:
        print("ERROR: can not identify strand for "+seq_name+ ". Ignoring this sequence")
        #sys.exit(0)
        continue
    else:
        spliced_ma=index_by_chr[chr_name].get_spliced([int(left)],[int(right)], strand=strand)
        alignment_file_name=alignment_folder.rstrip("/")+"/"+seq_name+".sto"
        AlignIO.write(spliced_ma, alignment_file_name, "stockholm")
        if seq_name in version_count_by_seq_name:
            version=version_count_by_seq_name[seq_name]+1
            seq_file_name=seq_folder.rstrip("/")+"/"+seq_name+"v"+str(version)+".fa"
            version_count_by_seq_name[seq_name]=version
        else:
            version_count_by_seq_name[seq_name]=1
            seq_file_name=seq_folder.rstrip("/")+"/"+seq_name+".fa"
        alt_ref_sequences=get_alternative_and_reference_sequences(row["Change"], spliced_ma, chr_name)
        SeqIO.write(alt_ref_sequences[0], seq_file_name, "fasta")
        refseq_file_name=refseq_folder.rstrip("/")+"/"+seq_name+".fa"
        SeqIO.write(alt_ref_sequences[1], refseq_file_name, "fasta")
        
print(maf_file_name + " processed")
