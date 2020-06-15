#run from a directory which contains:
# 1. ncVariation_dataset-testing_labels.csv --- table sheet (converted to csv) with curated benign and pathogenic data
# adjust the name of the table if needed:
path_csv_file='/Users/labadmin/Git/ncVarDB/data/ncVar_pathogenic.csv'
# 2. ncVariation_dataset-benign_UCSC_list.csv --- table sheet (converted to csv) with bening data based on minor allele frequencies
# adjust the name of the table if needed:
benign_csv_file='/Users/labadmin/Git/ncVarDB/data/ncVar_benign.csv'
# 3. process_large_table.py
# 4. process_benign_sheet.py
# 5. make_seq_refseq_and_alignments_from_maf.py
# 6. extract_bitscore.py
# 7. plot_dbss.R
# requiers: python v3.7, HMMER 3.2.1, R version 3.5.0 with package "ggplot2" installed, mafFetch (http://hgdownload.soe.ucsc.edu/admin/exe/)
# To use mafFetch one needs to add specification to $HOME/.hg.conf file (http://genome.ucsc.edu/goldenPath/help/mysql.html)


#python3 process_variants_in_table.py $path_csv_file overBed_pathogenic pathogenic.csv

#mkdir overBed_benign
#mkdir csvs_benign

#python3 process_variants_in_table_max1KperBedfile.py $benign_csv_file overBed_benign csvs_benign

type='pathogenic'

#touch $type'_error'

#mafFetch hg38 multiz100way overBed_$type $type'.maf'

#mkdir $type

#mkdir $type/sequences
#mkdir $type/alignments
#mkdir $type/models
#mkdir $type/scores
#mkdir $type/ref_sequences
#mkdir $type/ref_scores
#mkdir $type'_chr_mafs'


#python3 make_seq_refseq_and_alignments_from_maf.py $type'.csv' $type'.maf' $type'_chr_mafs' $type/sequences $type/alignments $type/ref_sequences > $type'_error'

#touch $type'_bitscores.csv'

#touch $type'_refbitscores.csv'

#for sequence in $type/sequences/*; do

#bname=$(basename "$sequence" .fa)

#bname_without_v=$(echo $bname | cut -f1 -d "v")

#hmmbuild $type/models/$bname_without_v'.hmm' $type/alignments/$bname_without_v'.sto' > $type/models/$bname_without_v'.hmmout'

#nhmmer $type/models/$bname_without_v'.hmm' $sequence > $type/scores/$bname'.out'

#nhmmer $type/models/$bname_without_v'.hmm' $type/ref_sequences/$bname_without_v.'fa' > $type/ref_scores/$bname'.out'

#python3 extract_bitscore.py $type/scores/$bname'.out' >> $type'_bitscores.csv'

#python3 extract_bitscore.py $type/ref_scores/$bname'.out' >> $type'_refbitscores.csv'

#done

type='benign'

touch $type'_error'

#mkdir $type

#mkdir $type/mafs

#for bed_file in overBed_benign/*; do

#bname=$(basename "$bed_file")

#mafFetch hg38 multiz100way $bed_file $type/mafs/$bname'.maf'

#done

mkdir $type/sequences
mkdir $type/alignments
mkdir $type/models
mkdir $type/scores
mkdir $type/ref_sequences
mkdir $type/ref_scores
mkdir $type'_chr_mafs'

for i in {1..10}; do

python3 make_seq_refseq_and_alignments_from_maf.py 'csvs_benign/'$i'.csv' 'benign/mafs/bed'$i'.maf' $type'_chr_mafs' $type/sequences $type/alignments $type/ref_sequences >> $type'_error'

done

#touch $type'_bitscores.csv'

#touch $type'_refbitscores.csv'

#for sequence in $type/sequences/*; do

#bname=$(basename "$sequence" .fa)

#bname_without_v=$(echo $bname | cut -f1 -d "v")

#hmmbuild $type/models/$bname_without_v'.hmm' $type/alignments/$bname_without_v'.sto' > $type/models/$bname_without_v'.hmmout'

#nhmmer $type/models/$bname_without_v'.hmm' $sequence > $type/scores/$bname'.out'

#nhmmer $type/models/$bname_without_v'.hmm' $type/ref_sequences/$bname_without_v.'fa' > $type/ref_scores/$bname'.out'

#python3 extract_bitscore.py $type/scores/$bname'.out' >> $type'_bitscores.csv'

#python3 extract_bitscore.py $type/ref_scores/$bname'.out' >> $type'_refbitscores.csv'

#done

#mkdir plots

#Rscript plot_dbss.R pathogenic_bitscores.csv benign_bitscores.csv pathogenic_refbitscores.csv benign_refbitscores.csv
