#pragma once



void DBSCAN(arma::mat& data, double dbscan_eps,int dbscan_minpts, std::ofstream& outfile) {

	int dim = data.n_rows;
	int num = data.n_cols;

	std::vector<std::vector<size_t>> neighbors;
	std::vector<std::vector<double>> distances;
	std::set<int> non_cores;
	std::set<int> noises;
	std::vector<bool> cores;
	
	mlpack::UnionFind ufset(num);
	
	int num_cores = 0;
	int num_pre_cluster;

#if (defined  PRINT_TIME) || (defined FILOUT_TIME)
#ifdef _WIN32
	LARGE_INTEGER t1, t2, t3,t4,t5,tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
#elif defined __linux__
	struct timeval t1, t2,t3,t4,t5;
	gettimeofday(&t1, NULL);
#endif
#endif
	EpsMinpts_Query<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> emq(data, dbscan_eps, dbscan_minpts);

	emq.DoQuery(neighbors,distances,cores);
	#if (defined  PRINT_TIME) || (defined FILOUT_TIME)
#ifdef _WIN32
	QueryPerformanceCounter(&t2);
#elif defined __linux__
	gettimeofday(&t2, NULL);
#endif
#endif
	/*mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::EuclideanDistance, arma::mat> nns(data);
	arma::Mat<size_t> nbrs;arma::mat dists;
	nns.Search(dbscan_minpts, nbrs, dists);
	for (int i = 0; i < num; i++) {
		bool xx = (dists(dbscan_minpts - 1, i) <= dbscan_eps);
		if (cores[i] != xx) {
			std::cout << "W" << std::endl;
		}
	}*/

	for (int i = 0; i < num; i++) {
		if (cores[i]) {
			for (int j = 0; j < neighbors[i].size(); j++) {
				if (cores[neighbors[i][j]]) ufset.Union(i, neighbors[i][j]);
			}
		}
	}


	//Use map<int,vec<int>> to collect the clusters
	std::unordered_map<int, std::vector<int>> mapIndexReps;
	for (int i = 0; i < num; i++) {
		if (!cores[i])continue;
		num_cores++;
		int rep = ufset.Find(i);
		if (mapIndexReps.find(rep) == mapIndexReps.end()) {
			std::vector<int> tmp; tmp.push_back(i);
			mapIndexReps.emplace(rep, tmp);
		} else {
			mapIndexReps.at(rep).push_back(i);
		}
	}
#if (defined  PRINT_TIME) || (defined FILOUT_TIME)
#ifdef _WIN32
	QueryPerformanceCounter(&t3);
#elif defined __linux__
	gettimeofday(&t3, NULL);
#endif
#endif
	num_pre_cluster=mapIndexReps.size();


	//Final Clustering process
	//Options: simple_traverse -/- merge_on_find. large_first_visit -/- no specific order 
#ifdef FINAL_SIMPLE_TRAVERSE
#ifdef AKNN_FINAL_CLUSTER
	final_simple_traverse_allknn(data, mapIndexReps, dbscan_eps, ufset);
#endif
#ifdef EPS_DETERMIN_FINAL_CLUSTER
	final_simple_traverse_eps_connectness(data, mapIndexReps, dbscan_eps, ufset);
	//final_simple_traverse_eps_connectness_RTree(data, mapIndexReps, dbscan_eps, ufset);
	//final_simple_traverse_eps_connectness_CoverTree(data, mapIndexReps, dbscan_eps, ufset);
#endif 
#endif


#ifdef FINAL_MERGE_ON_FIND
#ifdef LARGE_TREE_FIRST_VISIT
#ifdef AKNN_FINAL_CLUSTER
	final_merge_on_find_large_first_allknn(data, mapIndexReps, dbscan_eps, ufset);
#endif
#ifdef EPS_DETERMIN_FINAL_CLUSTER
	//final_merge_on_find_large_first_eps_connectness(data, mapIndexReps, dbscan_eps, ufset);
	final_merge_on_find_large_first_eps_connectness_batched(data, mapIndexReps, dbscan_eps, ufset);
#endif

#else
#ifdef AKNN_FINAL_CLUSTER
	final_merge_on_find_arbitrary_order_allknn(data, mapIndexReps, dbscan_eps, ufset);

#endif
#ifdef  EPS_DETERMIN_FINAL_CLUSTER
	//final_simple_traverse_eps_connectness(data, mapIndexReps, dbscan_eps, ufset);//KDTree impl.
	//final_merge_on_find_arbitrary_order_eps_connectness_CoverTree_ExpandOnQuery(data, mapIndexReps, dbscan_eps, ufset);
	//final_merge_on_find_arbitrary_order_eps_connectness(data, mapIndexReps, dbscan_eps, ufset);
	//final_merge_on_find_arbitrary_order_eps_connectness_batched(data, mapIndexReps, dbscan_eps, ufset);
	//final_merge_on_find_arbitrary_order_eps_connectness_RTree(data, mapIndexReps, dbscan_eps, ufset);
	//final_merge_on_find_arbitrary_order_eps_connectness_CoverTree(data, mapIndexReps, dbscan_eps, ufset);
	final_merge_on_find_arbitrary_order_eps_connectness_CoverTree_Lazy(data, mapIndexReps, dbscan_eps, ufset);
	
#endif
#endif
#endif

#ifdef PRINT_TIME
#ifdef _WIN32
	QueryPerformanceCounter(&t4);
	//std::cout << (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << std::endl;
#elif defined __linux__
	gettimeofday(&t4, NULL);
	//std::cout << (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) * 1.0 / 1000000 << std::endl;
#endif
#endif
#if defined FILOUT_TIME
#ifdef _WIN32
	QueryPerformanceCounter(&t4);
#elif defined __linux__
	gettimeofday(&t4, NULL);
#endif
#endif
	//TODO:finally add the non-cores to the cluster
	//Scan the non-core points, if there exists one core point in its kNN ( the kNN distance of one non-core point is >=dbscan_eps)
	//Then it can be added to the cluster.
	//for (int ncpoint : non_cores) {
	//	for (int i = 0; i < neighbors[ncpoint].size(); i++) {
	//		if (distances[ncpoint][i] <= dbscan_eps && cores[neighbors[ncpoint][i]]) {
	//			//mapIndexReps contains core_index
	//			ufset.Union(ncpoint, neighbors[ncpoint][i]);
	//			//mapIndexReps[ufset.Find(aknnRes[i, ncpoint])].push_back(ncpoint);
	//			break;
	//		}
	//	}
	//}
	for (size_t i = 0; i < num; i++) {
		if (cores[i]) continue;
		for (size_t j = 0; j < neighbors[i].size(); j++) {
			if (cores[neighbors[i][j]] && distances[i][j] <= dbscan_eps) {
				ufset.Union(i, neighbors[i][j]); //break;
			}
		}
	}
#if defined FILOUT_TIME
#ifdef _WIN32
	QueryPerformanceCounter(&t5);
#elif defined __linux__
	gettimeofday(&t5, NULL);
	//outfile<<num<<", "<<dim<<", "<<dbscan_minpts<<", "<<dbscan_eps<<", ";
	outfile<< (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) * 1.0 / 1000000 << ", ";
	outfile<< (t3.tv_sec - t2.tv_sec) + (t3.tv_usec - t2.tv_usec) * 1.0 / 1000000 << ", ";
	outfile<< (t4.tv_sec - t3.tv_sec) + (t4.tv_usec - t3.tv_usec) * 1.0 / 1000000 << ", ";
	outfile<< (t5.tv_sec - t4.tv_sec) + (t5.tv_usec - t4.tv_usec) * 1.0 / 1000000 << ", ";
	outfile<< (t5.tv_sec - t1.tv_sec) + (t5.tv_usec - t1.tv_usec) * 1.0 / 1000000 << ", ";
	outfile<<num_cores<<","<<num_pre_cluster<<", ";
#endif
#endif
#ifdef PRINT_NUM_FINAL_CLUSTERS
	mapIndexReps.clear();
	int final_clusters = 0;
	for (int i = 0; i < num; i++) {
		int rep = ufset.Find(i);
		if (mapIndexReps.find(rep) == mapIndexReps.end()) {
			std::vector<int> tmp; tmp.push_back(i);
			mapIndexReps.emplace(rep, tmp);
		} else {
			mapIndexReps.at(rep).push_back(i);
		}
	}
	for (auto& pair : mapIndexReps) {
		if (pair.second.size() >= dbscan_minpts)
			final_clusters++;
	}
	//std::cout << "num clusters after final_clustering: " << final_clusters<< std::endl;
	outfile << final_clusters<< std::endl;
#endif
}