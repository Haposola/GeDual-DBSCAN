// DBSCAN_StackSimulation.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "HeadOfHead.h"

int main(int argc, char** argv)
{
	//test_stack_simulation();

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
	generate_data(data, num, dim);
#endif

#ifdef TEST_SYNTHETIC
	num = 1000000; dim = 9; dbscan_minpts = 5; dbscan_eps = 1000000;
	generate_data(data, num, dim);
#endif

#ifdef TEST_FILE
	//std::string ifname = "E:\\datasets\\gangen\\gangen_convert09.txt";
	//std::string ifname = "E:\\datasets\\unigen\\uniint3.txt"; dbscan_minpts = 5; dbscan_eps = 900;//5-900,clusters 28369
	//std::string ifname = "E:\\datasets\\household\\household_multi.txt"; dbscan_minpts = 5; dbscan_eps = 10;//5-10,clusters 8427, 5-5 11736
	//std::string ifname = "E:\\datasets\\corel\\corel_int.txt"; dbscan_minpts = 5; dbscan_eps = 9000;
	std::string ifname = "E:\\datasets\\gangen\\gangen_convert09.txt"; dbscan_minpts = 5; dbscan_eps = 400;//5-10,clusters 8427, 5-5 11736
	//denThres = dbscan_eps;
	mlpack::data::Load(ifname, data);
	//test_all_range_search_clusters(data, dbscan_eps, dbscan_minpts);

	DBSCAN_StackSimulation xx(data, dbscan_eps, dbscan_minpts);
	xx.DoDBSCAN(std::cout);
	//test_rangesearch_cluster(data, dbscan_eps, dbscan_minpts);
#endif
#ifdef RUNEXP
	
	dbscan_minpts = atoi(argv[1]);
	dbscan_eps = atof(argv[2]);
	std::string ifname(argv[3]);
	std::ofstream ofile("DBSCAN_SS_time.txt",std::ios::app);
	mlpack::data::Load(ifname, data);
	DBSCAN_StackSimulation xx(data, dbscan_eps, dbscan_minpts);
	xx.DoDBSCAN(ofile);
	
#endif


}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
