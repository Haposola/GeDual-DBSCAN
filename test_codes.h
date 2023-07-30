#pragma once

void generate_data(arma::mat& data, int num, int dim) {
	//setsize(n_rows,n_cols)
	data.set_size(dim, num);
	//Generate random data
	std::srand((int)time(0));
	double range = num * 10 > INT_MAX ? INT_MAX : num * 10;
	//arma::mat is column-wise. a column is a data.
	for (int i = 0; i < num; i++) {
		for (int j = 0; j < dim; j++) {
			data(j, i) = rand() * 1.0 / RAND_MAX * range;
			data(j, i) += rand() * 1.0 / RAND_MAX;
		}
	}
}

void test_stack_simulation() {
	arma::mat data;
	mlpack::data::Load("E://datasets//unigen//uniint3.txt", data);
	scoreStack < mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>*> myStack;
	std::vector<size_t> oldFromNew;
	mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat> searchTree(data,oldFromNew);

	//arma::vec query = arma::zeros(data.n_rows);
	arma::vec query = data.col(0);
	myStack.push( &searchTree);

	double currentNNDist=INT_MAX; size_t currentNN;
	while (!myStack.isEmpty()) {
		//double topScore = myStack.topScore();
		auto topPtr = myStack.topElement();
		myStack.pop();
		//if (topScore > currentNNDist) continue;
		if (topPtr->IsLeaf()) {
			const size_t refEnd = topPtr->Begin() + topPtr->Count();
			for (size_t i = topPtr->Begin(); i < refEnd; ++i) {
				double dist = mlpack::EuclideanDistance::Evaluate(query, data.col(oldFromNew[i]));
				if ( dist< currentNNDist) {
					currentNNDist = dist;
					currentNN = oldFromNew[i];
				}
			}
			continue;
		}
		double leftScore = topPtr->Left()->MinDistance(query);
		double rightScore = topPtr->Right()->MinDistance(query);

		if (leftScore > rightScore) {//push large score first, then the small score will be popped first.
			if(leftScore)
			myStack.push( topPtr->Left());
			myStack.push( topPtr->Right());
		} else {
			myStack.push( topPtr->Right());
			myStack.push( topPtr->Left());
		}
	}
	std::cout << "NNDist simulated: " << currentNNDist << "\n";
	std::cout << "NNPoint simulated: " << currentNN << "\n";
	mlpack::NeighborSearch aa(data);
	arma::mat distances;
	arma::Mat<size_t> indexes;
	aa.Search(query, 1, indexes, distances);
	std::cout << "NNDist true: " << distances(0) << "\n";
	std::cout << "NNPoint true: " << indexes(0) << "\n";
}


void test_rangesearch_cluster(arma::mat& data,double dbscan_eps, int dbscan_minpts) {
	LARGE_INTEGER t1, t2, t3, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);

	mlpack::RangeSearch rr(data,false,true);//naive=false,single=true
	std::vector<std::vector<size_t>> neighbors;
	std::vector<std::vector<double>> distances;
	rr.Search(mlpack::Range(0,dbscan_eps), neighbors, distances);

	bool* cores = new bool[data.n_cols];
	mlpack::UnionFind ufset(data.n_cols);
	int numcores = 0;
	for (size_t i = 0; i < data.n_cols; i++) {
		cores[i] = (neighbors[i].size() >= dbscan_minpts - 1);
	}
	for (size_t i = 0; i < data.n_cols; i++) {
		if (cores[i]) {
			numcores++;
			for (size_t j = 0; j < neighbors[i].size(); j++) {
				if (cores[neighbors[i][j]]) ufset.Union(i, neighbors[i][j]);
			}
		}
	}
	for (size_t i = 0; i < data.n_cols; i++) {
		if (cores[i])continue;
		for (size_t j = 0; j < neighbors[i].size(); j++) {
			if (cores[neighbors[i][j]]) {
				ufset.Union(i, neighbors[i][j]);
				//break;
			}
		}
	}
	std::unordered_map<int, std::vector<int>> mapIndexReps;
	int final_clusters = 0;
	for (int i = 0; i < data.n_cols; i++) {
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
	QueryPerformanceCounter(&t2);
	std::cout << "num cores by rangesearch: " << numcores << "\n";
	std::cout << "num clusters by rangesearch: " << final_clusters << "\n";
	std::cout << "rangesearch time: " << (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << "\n";
}