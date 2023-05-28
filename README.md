# EMEC-DBSCAN

EMEC-DBSCAN: DBSCAN algorithm using Eps-Minpts and Eps-Connectivity query

The Eps-Minpts query is based on MLPACK https://mlpack.org/

So you must have MLPACK installed to compile these files here.

# Input data format

any format that compatible with arma::mat

which is the same format used in MLPACK.

e.g., .csv file (with no header and no id).

e.g., three point in dimension 4:

3.0 4.0 5.0 2.0

1.0 2.0 3.0 2.0

6.0 9.0 8.0 7.0

# command line parameters:

$1 $2 $3

$1: dbscan_minpts, $2: dbscan_eps, $3: output_file_name

Change the command line parameters in DBSCAN_DualTraversals.cpp::main() if your want.
