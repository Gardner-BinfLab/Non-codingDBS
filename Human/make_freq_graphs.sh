# run this script after running analyse.sh and analyse_likely-path.sh
# run from a directory which contains:
# 1. ncVariation_dataset-testing_labels.csv --- table sheet (converted to csv) with curated benign and pathogenic data
# adjust the name of the table if needed:
csv_file="ncVariation_dataset-testing_labels.csv"
# 2. ncVariation_dataset-benign_UCSC_list.csv --- table sheet (converted to csv) with bening data based on minor allele frequencies
# adjust the name of the table if needed:
benign_csv_file="ncVariation_dataset-benign_UCSC_list.csv"
# 3. extract_clinVar_freq_benign.py
# 4. pathogenic_bitscores.csv likely-path_bitscores.csv benign_bitscores.csv pathogenic_refbitscores.csv likely-path_refbitscores.csv benign_refbitscores.csv benign_clinVar_freq.csv
# 5. plot_dbs_densities_dif_freq.R

python3 extract_clinVar_freq_benign.py $csv_file $benign_csv_file benign_clinVar_freq.csv > benign_error

mkdir plots_dif_freq

Rscript plot_dbs_densities_dif_freq.R pathogenic_bitscores.csv likely-path_bitscores.csv benign_bitscores.csv pathogenic_refbitscores.csv likely-path_refbitscores.csv benign_refbitscores.csv benign_clinVar_freq.csv

