﻿// GeDual-DBSCAN-transite.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "HeadOfHead.h"

int main(int argc, char** argv) {
	//test_exptime();
//4 parameters: dbscan_eps dbscan_minpts input_file_in_arma::fmat output_file
	double dbscan_eps = atof(argv[1]);
	int dbscan_minpts = atoi(argv[2]);
	arma::Mat<double> data;
	data.load(argv[3]);

	std::cout << "num " << data.n_cols << " dim " << data.n_rows << "\n";
	std::ofstream outfile(argv[4], std::ios::app);
	//print_mlpack_numclusters(data, dbscan_eps, dbscan_minpts);
	GeDual_DBSCAN xx(data, dbscan_eps, dbscan_minpts);
	xx.DoDBSCAN(outfile);
	std::cout << "Hello World!\n";
}