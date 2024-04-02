# GeDual-DBSCAN

GeDual-DBSCAN: A generalized dual-tree algorithm for DBSCAN

The implementation is based on KDTree of MLPACK https://mlpack.org/

So you must have MLPACK installed to compile these files here.

## Algorithm variants

This repo contains tree variants, postfixed by Basic, InterCases and Transitions, respectively.

The last variant, GeDual-DBSCAN-Transitions, is the most powerful.

## Input data format

any format that compatible with arma::mat

which is the same format used in MLPACK.

e.g., .csv file (with no header and no id).

e.g., three point in dimension 4:

3.0 4.0 5.0 2.0

1.0 2.0 3.0 2.0

6.0 9.0 8.0 7.0

# command line parameters:

$1 $2 $3 $4

$1: dbscan_eps, $2: dbscan_minpts, $3: input_file_neme, $4: output_file_name

Change the command line parameters in GeDual-DBSCAN-postfix.cpp::main() if your want.
