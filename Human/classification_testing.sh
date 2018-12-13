# run this script after running analyse.sh and analyse_likely-path.sh
# run from a directory which contains:
# 1. logistic_regression.R
# 2. pathogenic_bitscores.csv likely-path_bitscores.csv benign_bitscores.csv pathogenic_refbitscores.csv likely-path_refbitscores.csv benign_refbitscores.csv
# requiers: R version 3.5.0 with packages: "pROC" and "boot" installed

Rscript logistic_regression.R pathogenic_bitscores.csv likely-path_bitscores.csv benign_bitscores.csv pathogenic_refbitscores.csv likely-path_refbitscores.csv benign_refbitscores.csv
