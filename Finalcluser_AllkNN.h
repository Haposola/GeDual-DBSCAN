#pragma once

void final_simple_traverse_allknn(arma::mat& data, std::unordered_map<int, std::vector<int>>& mapIndexReps, double denThres, mlpack::UnionFind& ufset) {
	typedef mlpack::KDTree<EucDistHere, mlpack::NeighborSearchStat<mlpack::NearestNeighborSort>, arma::mat> KDTreeHere;

	int num = data.n_cols;
	int dim = data.n_rows;

	std::vector<std::pair<int, KDTreeHere*>> vecRepTree;
	for (auto& pair : mapIndexReps) {
		arma::mat tmpmat(dim, pair.second.size());
		for (int i = 0; i < pair.second.size(); i++) {
			tmpmat.col(i) = data.col(pair.second[i]);
		}
		vecRepTree.push_back(std::pair<int, KDTreeHere*>(pair.first, new  KDTreeHere(tmpmat)));
	}
	for (int i = 0; i < vecRepTree.size(); i++) {
		mlpack::NeighborSearch < mlpack::NearestNeighborSort, EucDistHere, arma::mat> a(*vecRepTree[i].second, mlpack::DUAL_TREE_MODE);
		for (int j = i + 1; j < vecRepTree.size(); j++) {
			arma::Mat<size_t> result;
			arma::mat resultDis;
			a.Search(*vecRepTree[j].second, 1, result, resultDis);
			for (int l = 0; l < resultDis.n_elem; l++) {
				if (resultDis.at(l) <= denThres) {
					ufset.Union(vecRepTree[i].first, vecRepTree[j].first);
					break;
				}
			}
		}
	}

}

void final_merge_on_find_large_first_allknn(arma::mat& data, std::unordered_map<int, std::vector<int>>& mapIndexReps, double denThres, mlpack::UnionFind& ufset) {
	typedef mlpack::KDTree<EucDistHere, mlpack::NeighborSearchStat<mlpack::NearestNeighborSort>, arma::mat> KDTreeHere;

	int num = data.n_cols;
	int dim = data.n_rows;

	//The map should be mat - pair<Tree,rep_point>
	std::map<std::pair<int, int>, KDTreeHere*, std::greater<std::pair<int, int>>> mapSizeRepTree;


	//std::set < arma::mat, arma_mat_size_comparator> clustered;
	for (auto& pairReps : mapIndexReps) {
		arma::mat tmpmat(dim, pairReps.second.size());
		for (int i = 0; i < pairReps.second.size(); i++) {
			tmpmat.col(i) = data.col(pairReps.second[i]);
		}
		mapSizeRepTree.emplace(std::pair< int, int>(pairReps.second.size(), pairReps.first), new KDTreeHere(std::move(tmpmat)));
	}
	//std::vector<KDTreeHere*> examined;

	while (!mapSizeRepTree.size() > 1) {
		int lsize = mapSizeRepTree.begin()->first.first;
		int rep = mapSizeRepTree.begin()->first.second;
		KDTreeHere* currentTree = mapSizeRepTree.begin()->second;

		mapSizeRepTree.erase(mapSizeRepTree.begin());
		bool modified = false;

		mlpack::NeighborSearch < mlpack::NearestNeighborSort, EucDistHere, arma::mat> a(*currentTree, mlpack::DUAL_TREE_MODE);
		std::vector<KDTreeHere*> to_merge;
		int merge_size = lsize;
		bool modefied = false;

		for (auto iter = mapSizeRepTree.begin(); iter != mapSizeRepTree.end(); ) {
			arma::Mat<size_t> result;
			arma::mat resultDis;
			a.Search(*iter->second, 1, result, resultDis);
			bool res = false;
			for (int l = 0; l < resultDis.n_elem; l++) {
				if (resultDis.at(l) <= denThres) {
					res = true; break;
				}
			}
			if (res) {
				//MERGE and reinsert
				modefied = true;
				to_merge.push_back(iter->second);
				merge_size += iter->first.first;
				ufset.Union(rep, iter->first.second);
				iter = mapSizeRepTree.erase(iter);
			} else iter++;
		}
		if (modified) {
			arma::mat newmat(dim, merge_size);
			int col_index = 0;
			for (int i = 0; i < currentTree->Dataset().n_cols; i++) {
				newmat.col(col_index++) = currentTree->Dataset().col(i);
			}
			for (auto merged_tree : to_merge) {
				for (int i = 0; i < merged_tree->Dataset().n_cols; i++) {
					newmat.col(col_index++) = merged_tree->Dataset().col(i);
				}
			}
			mapSizeRepTree.emplace(std::pair<int, int>(col_index, rep), new KDTreeHere(newmat));
		}

	}


}


void final_merge_on_find_arbitrary_order_allknn(arma::mat& data, std::unordered_map<int, std::vector<int>>& mapIndexReps, double denThres, mlpack::UnionFind& ufset) {
	typedef mlpack::KDTree<EucDistHere, mlpack::NeighborSearchStat<mlpack::NearestNeighborSort>, arma::mat> KDTreeHere;

	int num = data.n_cols;
	int dim = data.n_rows;

	//It is feasible since representitive is unique;
	std::unordered_map<int, KDTreeHere*>  mapRepTree;

	//std::set < arma::mat, arma_mat_size_comparator> clustered;
	for (auto& pairReps : mapIndexReps) {
		arma::mat tmpmat(dim, pairReps.second.size());
		for (int i = 0; i < pairReps.second.size(); i++) {
			tmpmat.col(i) = data.col(pairReps.second[i]);
		}
		mapRepTree.emplace(pairReps.first, new KDTreeHere(tmpmat));
	}


	while (mapRepTree.size() > 1) {
		int rep = mapRepTree.begin()->first;
		KDTreeHere* currentTree = mapRepTree.begin()->second;

		mapRepTree.erase(mapRepTree.begin());

		mlpack::NeighborSearch < mlpack::NearestNeighborSort, EucDistHere, arma::mat> a(*currentTree, mlpack::DUAL_TREE_MODE);

		arma::mat newmat;
		std::vector<KDTreeHere*> to_merge;
		int merge_size = currentTree->Dataset().n_cols;
		bool modefied = false; //int founded = 0;
		for (auto iter = mapRepTree.begin(); iter != mapRepTree.end(); ) {
			arma::Mat<size_t> result;
			arma::mat resultDis;
			a.Search(*iter->second, 1, result, resultDis);
			bool res = false;
			for (int l = 0; l < resultDis.n_elem; l++) {
				if (resultDis.at(l) <= denThres) {
					res = true; break;
				}
			}
			if (res) {
				//founded++;
				//MERGE and reinsert
				modefied = true;
				to_merge.push_back(iter->second);
				merge_size += iter->second->Dataset().n_cols;
				ufset.Union(rep, iter->first);
				iter = mapRepTree.erase(iter);
			} else iter++;
		}
		if (modefied) {
			arma::mat newmat(dim, merge_size);
			int col_index = 0;
			for (int i = 0; i < currentTree->Dataset().n_cols; i++) {
				newmat.col(col_index++) = currentTree->Dataset().col(i);
			}
			for (auto merged_tree : to_merge) {
				for (int i = 0; i < merged_tree->Dataset().n_cols; i++) {
					newmat.col(col_index++) = merged_tree->Dataset().col(i);
				}
			}
			mapRepTree.emplace(rep, new KDTreeHere(newmat));
		}
		//std::cout << founded << std::endl;
		//else examined.push_back(currentPair.first);

	}

}

