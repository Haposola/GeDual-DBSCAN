#pragma once

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
void test_EpsMinpts_correct() {

	arma::mat data;
	double dbscan_eps = 100000;
	int dbscan_minpts = 5;
	mlpack::data::Load("E://datasets//linux1.txt", data);
	std::cout << "data num: " << data.n_cols << std::endl;
	std::cout << "data dim: " << data.n_rows << std::endl;
	EpsMinpts_Query<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> aa(data, dbscan_eps, dbscan_minpts);
	std::vector<std::vector<size_t>> neighbors;
	std::vector<std::vector<double>> distances;
	std::vector<bool> isCore;
	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	aa.DoQuery(neighbors, distances, isCore);
	QueryPerformanceCounter(&t2);
	std::cout << "my time: " << (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << std::endl;
	std::cout << data.n_cols;
	//mlpack::data::Load("E://datasets//linux1.txt", data);
	 
	mlpack::RangeSearch<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> rse(data);
	QueryPerformanceCounter(&t1);
	rse.Search(mlpack::Range(0.0, dbscan_eps), neighbors, distances);
	QueryPerformanceCounter(&t2);
	std::cout << "all-range-search time: " << (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << std::endl;


}
void test_EpsConnected_correct() {

	for (int rep = 0; rep < 5; rep++) {
		arma::mat data1;
		arma::mat data2;
		int num = 10000, dim = 3, dbscan_minpts = 5;
		double dbscan_eps = 2000;
		generate_data(data1, num, dim);
		generate_data(data2, num, dim);
		mlpack::NeighborSearch< mlpack::NearestNeighborSort, mlpack::EuclideanDistance, arma::mat > bbb(data1);
		arma::mat distances;
		arma::Mat<size_t> neighbors;
		bbb.Search(data2, 1, neighbors, distances);
		bool aknnres = false;
		for (int i = 0; i < distances.n_elem; i++) {
			if (distances.at(i) <= dbscan_eps) aknnres = true;
		}
		EpsConnected_Query<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree>  aa(data1, data2, dbscan_eps);

		aa.DoQuery();

		if (aa.result != aknnres) std::cout << "Hello World!\n";
	}

}

void test_all_range_search_clusters(arma::mat& data,double dbscan_eps,int dbscan_minpts) {
	size_t num = data.n_cols;
	mlpack::RangeSearch rs(data);
	std::vector<std::vector<size_t>> neighbors;
	std::vector<std::vector<double>> distances;
	rs.Search(mlpack::Range(0.0, dbscan_eps), neighbors, distances);
	std::vector<bool> cores(num);
	mlpack::UnionFind ufset(num);
	for (size_t i = 0; i < num; i++) {
		if (neighbors[i].size() >= dbscan_minpts)
			cores[i] = true;
	}
	for (size_t i = 0; i < num; i++) {
		if (!cores[i]) continue;
		for (size_t j = 0; j < neighbors[i].size(); j++) {
			if (cores[neighbors[i][j]]) {
				ufset.Union(i, neighbors[i][j]);
			}
		}
	}
	for (size_t i = 0; i < num; i++) {
		if (cores[i]) continue;
		for (size_t j = 0; j < neighbors[i].size(); j++) {
			if (cores[neighbors[i][j]] && distances[i][j] <= dbscan_eps ) {
				ufset.Union(i, neighbors[i][j]); //break;
			}
		}
	}
	std::unordered_map<int, std::vector<int>>mapIndexReps;
	
	for (int i = 0; i < num; i++) {
		int rep = ufset.Find(i);
		if (mapIndexReps.find(rep) == mapIndexReps.end()) {
			std::vector<int> tmp; tmp.push_back(i);
			mapIndexReps.emplace(rep, tmp);
		} else {
			mapIndexReps.at(rep).push_back(i);
		}
	}
	int final_clusters = 0;
	for (auto& pair : mapIndexReps) {
		if (pair.second.size() >= dbscan_minpts)
			final_clusters++;
	}
	std::cout << "num clusters after final_clustering: " << final_clusters << std::endl;
	
}