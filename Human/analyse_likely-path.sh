# run this script only after analyse.sh is run
# run from a directory which contains:
#   1. process_large_table_likely-path.py
#   2. make_seq_refseq_and_alignments_from_maf.py
#   3. extract_bitscore.py
#   3. plot_dbss_likely-path.R
#   4. ncVariation_dataset-testing_labels.csv  --- table sheet (converted to csv) with curated benign and pathogenic data.
# adjust the name of the table if needed
csv_file="ncVariation_dataset-testing_labels.csv"
#   5. pathogenic_bitscores.csv, benign_bitscores.csv, pathogenic_refbitscores.csv, benign_refbitscores.csv
# requiers: python v3.7, R version 3.5.0 with packages: "ggplot2", "grid", "gridExtra" installed, mafFetch (http://hgdownload.soe.ucsc.edu/admin/exe/)
# To use mafFetch one needs to add specification to $HOME/.hg.conf file (http://genome.ucsc.edu/goldenPath/help/mysql.html)

python3 process_large_table_likely-path.py $csv_file overBed_likely-path likely-path.csv

type='likely-path'

mafFetch hg38 multiz100way overBed_$type $type'.maf'

mkdir $type

mkdir $type/sequences
mkdir $type/alignments
mkdir $type/models
mkdir $type/scores
mkdir $type/ref_sequences
mkdir $type/ref_scores
mkdir $type'_chr_mafs'

python3 make_seq_refseq_and_alignments_from_maf.py $type'.csv' $type'.maf' $type/sequences $type/alignments $type/ref_sequences

touch $type'_bitscores.csv'

touch $type'_refbitscores.csv'

for sequence in $type/sequences/*; do

bname=$(basename "$sequence" .fa)

bname_without_v=$(echo $bname | cut -f1 -d "v")

hmmbuild $type/models/$bname_without_v'.hmm' $type/alignments/$bname_without_v'.sto' > $type/models/$bname_without_v'.hmmout'

nhmmer $type/models/$bname_without_v'.hmm' $sequence > $type/scores/$bname'.out'

nhmmer $type/models/$bname_without_v'.hmm' $type/ref_sequences/$bname_without_v.'fa' > $type/ref_scores/$bname'.out'

python3 extract_bitscore.py $type/scores/$bname'.out' >> $type'_bitscores.csv'

python3 extract_bitscore.py $type/ref_scores/$bname'.out' >> $type'_refbitscores.csv'

done

mkdir plots_likely-path

Rscript plot_dbss_likely-path.R pathogenic_bitscores.csv likely-path_bitscores.csv benign_bitscores.csv pathogenic_refbitscores.csv likely-path_refbitscores.csv benign_refbitscores.csv
