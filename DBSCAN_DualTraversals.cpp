// DBSCAN_DualTraversals.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "HeadOfHead.h"
void generate_data(arma::mat& data, int num, int dim) {
	//setsize(n_rows,n_cols)
	data.set_size(dim, num);
	//Generate random data
	std::srand((int)time(0));
	double range = num * 10 > RAND_MAX ? RAND_MAX : num * 10;
	//arma::mat is column-wise. a column is a data.
	for (int i = 0; i < num; i++) {
		for (int j = 0; j < dim; j++) {
			data(j, i) = rand() * 1.0 / RAND_MAX * range;
			data(j, i) += rand() * 1.0 / RAND_MAX;
		}
	}
}
int main(int argc, char**argv )
{
	mlpack::RangeSearch aa;
	//test_EpsMinpts_correct();
	arma::mat data;
	int num, dim, dbscan_minpts;
	double dbscan_eps;
	
	//test_EpsMinpts_query_correct();
	//test_EpsMinpts_dual_correct();
#ifdef TEST_CMD
	//if cmdline, arguments
	//num dim dbscan_minpts dbscan_eps
	num = atoi(argv[1]);
	dim = atoi(argv[2]);
	dbscan_minpts = atoi(argv[3]);
	dbscan_eps = atof(argv[4]);
	generate_data(data,num,dim);
	
#endif

#ifdef TEST_SYNTHETIC
	num = 1000000; dim = 3; dbscan_minpts = 5; dbscan_eps = 2000;
	generate_data(data, num, dim);
#endif

#ifdef TEST_FILE
	dbscan_minpts = atoi(argv[1]);
    dbscan_eps = atof(argv[2]);
    std::string ifname = argv[3];
    std::string outname("EMEC_file_time.txt");
    std::ofstream outfile(outname,std::ios::app);
    //generateData(data, num, dim);
    //data::Save("E://datasets//NN/VLG//village_2_50000.txt", input);

    mlpack::data::Load(ifname, data);
	
	denThres = dbscan_eps;
#endif





outfile << data.n_cols << ", " << data.n_rows << ", "<<dbscan_minpts<<", "<<dbscan_eps<<", ";
	//std::cout << outname << std::endl;
	/**/
	//bool loaded = data::Load("E://datasets//NN/VLG//village_2_50000.txt", input);
	//if (!loaded) return -1;
	//cout << "load success";
	//point set
	//string ifile = "E:/datasets/NN/VLG/village_2_50000";
	//string ifile = "E:/datasets/NN/EQ/BiData0";

	//for(int i=0;i<1;i++)	test_EpsMinpts_query_correct();
	DBSCAN(data, dbscan_eps, dbscan_minpts, outfile);



    //std::cout << "Hello World!\n";
}

