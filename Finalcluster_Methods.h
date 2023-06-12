#pragma once

// typedef mlpack::KDTree<mlpack::EuclideanDistance, EpsConnectedStat, arma::mat> KDTreeECQ;
typedef mlpack::KDTree<EucDistHere, mlpack::EmptyStatistic, arma::mat> KDTreeHere;

void final_simple_traverse_eps_connectness(arma::mat &data, std::unordered_map<int, std::vector<int>> &mapIndexReps, double denThres, mlpack::UnionFind &ufset)
{
	int num = data.n_cols;
	int dim = data.n_rows;

	std::vector<std::pair<int, KDTreeHere *>> vecRepTree;
	// Pre_build the trees
#ifdef PRINT_BUILD_TREE_TIME
	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
#endif
	for (auto &pair : mapIndexReps)
	{
		arma::mat tmpmat(dim, pair.second.size());
		for (int i = 0; i < pair.second.size(); i++)
		{
			tmpmat.col(i) = data.col(pair.second[i]);
		}
		vecRepTree.push_back(std::pair<int, KDTreeHere *>(pair.first, new KDTreeHere(tmpmat)));
	}
#ifdef PRINT_BUILD_TREE_TIME
	QueryPerformanceCounter(&t2);
	std::cout << "build tree time: " << (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << std::endl;
#endif
	for (int i = 0; i < vecRepTree.size(); i++)
	{
		for (int j = i + 1; j < vecRepTree.size(); j++)
		{
			// EpsConnected_Query<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> ecq(vecRepTree[i].second, vecRepTree[j].second, denThres);
			// ecq.DoQuery();
			// if (ecq.result)
			bool found = false;
			eps_connectness_BSTree_dual_DFS(vecRepTree[i].second, vecRepTree[j].second, denThres, found);
			if (found)
				ufset.Union(vecRepTree[i].first, vecRepTree[j].first);
		}
	}
}

void final_simple_traverse_eps_connectness_RTree(arma::mat &data, std::unordered_map<int, std::vector<int>> &mapIndexReps, double denThres, mlpack::UnionFind &ufset)
{
	int num = data.n_cols;
	int dim = data.n_rows;
	typedef mlpack::RPlusTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat> RTreeHere;
	std::vector<std::pair<int, RTreeHere *>> vecRepTree;
	// Pre_build the trees
	for (auto &pair : mapIndexReps)
	{
		arma::mat tmpmat(dim, pair.second.size());
		for (int i = 0; i < pair.second.size(); i++)
		{
			tmpmat.col(i) = data.col(pair.second[i]);
		}
		vecRepTree.push_back(std::pair<int, RTreeHere *>(pair.first, new RTreeHere(tmpmat)));
	}
	for (int i = 0; i < vecRepTree.size(); i++)
	{
		for (int j = i + 1; j < vecRepTree.size(); j++)
		{
			// EpsConnected_Query<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> ecq(vecRepTree[i].second, vecRepTree[j].second, denThres);
			// ecq.DoQuery();
			// if (ecq.result)
			bool found = false;
			eps_connectness_RTree_dual_DFS(*vecRepTree[i].second, *vecRepTree[j].second, denThres, found);
			if (found)
				ufset.Union(vecRepTree[i].first, vecRepTree[j].first);
		}
	}
}
void final_simple_traverse_eps_connectness_CoverTree(arma::mat &data, std::unordered_map<int, std::vector<int>> &mapIndexReps, double denThres, mlpack::UnionFind &ufset)
{
	int num = data.n_cols;
	int dim = data.n_rows;
	typedef mlpack::CoverTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat> CoverTreeHere;
	std::vector<std::pair<int, CoverTreeHere *>> vecRepTree;
	// Pre_build the trees
	for (auto &pair : mapIndexReps)
	{
		arma::mat tmpmat(dim, pair.second.size());
		for (int i = 0; i < pair.second.size(); i++)
		{
			tmpmat.col(i) = data.col(pair.second[i]);
		}
		vecRepTree.push_back(std::pair<int, CoverTreeHere *>(pair.first, new CoverTreeHere(tmpmat)));
	}
	for (int i = 0; i < vecRepTree.size(); i++)
	{
		for (int j = i + 1; j < vecRepTree.size(); j++)
		{
			// EpsConnected_Query<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> ecq(vecRepTree[i].second, vecRepTree[j].second, denThres);
			// ecq.DoQuery();
			// if (ecq.result)
			bool found = false;

			eps_connectness_dual_CoverTree(*vecRepTree[i].second, *vecRepTree[j].second, found);
			if (found)
				ufset.Union(vecRepTree[i].first, vecRepTree[j].first);
		}
	}
}

void final_merge_on_find_arbitrary_order_eps_connectness(arma::mat &data, std::unordered_map<int, std::vector<int>> &mapIndexReps, double denThres, mlpack::UnionFind &ufset)
{

	int num = data.n_cols;
	int dim = data.n_rows;

	std::unordered_map<int, KDTreeHere *> mapRepTree;
	// std::unordered_map<int, KDTreeECQ*>  mapRepTree;
	// std::set < arma::mat, arma_mat_size_comparator> clustered;
	for (auto &pairReps : mapIndexReps)
	{
		arma::mat tmpmat(dim, pairReps.second.size());
		for (int i = 0; i < pairReps.second.size(); i++)
		{
			tmpmat.col(i) = data.col(pairReps.second[i]);
		}
		mapRepTree.emplace(pairReps.first, new KDTreeHere(tmpmat));
		// mapRepTree.emplace(pairReps.first, new KDTreeECQ(tmpmat));
	}

	while (!mapRepTree.empty())
	{

		int rep = mapRepTree.begin()->first;

		auto currentTree = mapRepTree.begin()->second;
		mapRepTree.erase(mapRepTree.begin());

		arma::mat newmat;
		std::vector<KDTreeHere *> to_merge;
		// std::vector<KDTreeECQ*> to_merge;
		int merge_size = currentTree->Dataset().n_cols;
		bool modefied = false;

		for (auto iter = mapRepTree.begin(); iter != mapRepTree.end();)
		{
			// EpsConnected_Query<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> ecq(currentTree, iter->second, denThres);
			// ecq.DoQuery();
			bool found = false;
			eps_connectness_BSTree_dual_DFS(currentTree, iter->second, denThres, found);
			// if(ecq.result){
			if (found)
			{
				// MERGE and reinsert

				modefied = true;
				to_merge.push_back(iter->second);
				merge_size += iter->second->Dataset().n_cols;
				ufset.Union(rep, iter->first);
				iter = mapRepTree.erase(iter);
			}
			else
				iter++;
		}
		if (modefied)
		{
			arma::mat newmat(dim, merge_size);
			int col_index = 0;
			for (int i = 0; i < currentTree->Dataset().n_cols; i++)
			{
				newmat.col(col_index++) = currentTree->Dataset().col(i);
			}
			for (auto merged_tree : to_merge)
			{
				for (int i = 0; i < merged_tree->Dataset().n_cols; i++)
				{
					newmat.col(col_index++) = merged_tree->Dataset().col(i);
				}
			}
			mapRepTree.emplace(rep, new KDTreeHere(newmat));
			// mapRepTree.emplace(rep, new KDTreeECQ(newmat));
		}

		// else examined.push_back(currentPair.first);
	}
}

void final_merge_on_find_arbitrary_order_eps_connectness_RTree(arma::mat &data, std::unordered_map<int, std::vector<int>> &mapIndexReps, double denThres, mlpack::UnionFind &ufset)
{

	typedef mlpack::RPlusTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat> RTreeHere;
	int num = data.n_cols;
	int dim = data.n_rows;

	std::unordered_map<int, RTreeHere *> mapRepTree;

	for (auto &pairReps : mapIndexReps)
	{
		arma::mat tmpmat(dim, pairReps.second.size());
		for (int i = 0; i < pairReps.second.size(); i++)
		{
			tmpmat.col(i) = data.col(pairReps.second[i]);
		}
		mapRepTree.emplace(pairReps.first, new RTreeHere(tmpmat));
		// mapRepTree.emplace(pairReps.first, new KDTreeECQ(tmpmat));
	}

	while (!mapRepTree.empty())
	{

		int rep = mapRepTree.begin()->first;

		auto currentTree = mapRepTree.begin()->second;
		mapRepTree.erase(mapRepTree.begin());

		arma::mat newmat;
		std::vector<RTreeHere *> to_merge;
		// std::vector<KDTreeECQ*> to_merge;
		int merge_size = currentTree->Dataset().n_cols;
		bool modefied = false;

		for (auto iter = mapRepTree.begin(); iter != mapRepTree.end();)
		{
			// EpsConnected_Query<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> ecq(currentTree, iter->second, denThres);
			// ecq.DoQuery();
			bool found = false;
			eps_connectness_RTree_dual_DFS(*currentTree, *iter->second, denThres, found);
			// if(ecq.result){
			if (found)
			{
				// MERGE and reinsert

				modefied = true;
				to_merge.push_back(iter->second);
				// merge_size += iter->second->Dataset().n_cols;
				ufset.Union(rep, iter->first);
				iter = mapRepTree.erase(iter);
			}
			else
				iter++;
		}
		if (modefied)
		{

			arma::mat &newmat = currentTree->Dataset();
			size_t oldsize = newmat.n_cols;
			size_t newsize = newmat.n_cols;
			for (auto merged : to_merge)
			{
				size_t addsize = merged->Dataset().n_cols;
				newmat.set_size(newmat.n_rows, newsize + addsize);

				for (size_t i = 0; i < addsize; i++)
					newmat.col(newsize + i) = merged->Dataset().col(i);

				newsize += addsize;
			}
			for (size_t i = oldsize; i < newsize; i++)
				currentTree->InsertPoint(i); // Hopefully to be correct;

			mapRepTree.emplace(rep, currentTree);
			// mapRepTree.emplace(rep, new KDTreeECQ(newmat));
		}

		// else examined.push_back(currentPair.first);
	}
}

void final_merge_on_find_arbitrary_order_eps_connectness_CoverTree(arma::mat &data, std::unordered_map<int, std::vector<int>> &mapIndexReps, double denThres, mlpack::UnionFind &ufset)
{

	typedef mlpack::CoverTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat> CoverTreeHere;
	int num = data.n_cols;
	int dim = data.n_rows;

	// std::unordered_map<int, CoverTreeHere*>  mapRepTree;
	std::vector<std::pair<int, CoverTreeHere *>> vecRepTree;
	vecRepTree.reserve(mapIndexReps.size());
	size_t ind = 0;
#ifdef PRINT_BUILD_TREE_TIME
	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
#endif
	for (auto &pairReps : mapIndexReps)
	{
		arma::mat tmpmat(dim, pairReps.second.size());
		for (int i = 0; i < pairReps.second.size(); i++)
		{
			tmpmat.col(i) = data.col(pairReps.second[i]);
		}
		vecRepTree[ind++] = std::make_pair(pairReps.first, new CoverTreeHere(std::move(tmpmat), 1.3));
	}
#ifdef PRINT_BUILD_TREE_TIME
	QueryPerformanceCounter(&t2);
	std::cout << "build tree time: " << (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << std::endl;
#endif
	if (vecRepTree.size() <= 1)
		return;
	std::vector<std::pair<int, CoverTreeHere *>> newvec;
	bool modefied = false;
	int rep = vecRepTree[0].first;
	CoverTreeHere *currentTree = vecRepTree[0].second;

	while (vecRepTree.size() > 1)
	{
		newvec.reserve(vecRepTree.size());
		size_t newindex = 0;
		size_t index = (modefied ? 0 : 1);
		std::vector<CoverTreeHere *> to_merge;

		modefied = false;

		for (; index < vecRepTree.size(); index++)
		{
			// EpsConnected_Query<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> ecq(currentTree, iter->second, denThres);
			// ecq.DoQuery();
			bool found = false;
			eps_connectness_dual_CoverTree(*currentTree, *vecRepTree[index].second, found);
			// if(ecq.result){
			if (found)
			{
				modefied = true;
				to_merge.push_back(vecRepTree[index].second);

				ufset.Union(rep, vecRepTree[index].first);
			}
			else
			{
				newvec[newindex++] = std::move(vecRepTree[index]);
			}
		}
		std::swap(vecRepTree, newvec);
		newvec.clear();
		if (modefied)
		{
			CoverTreeHere *merged = to_merge[0];
			for (size_t i = 1; i < to_merge.size(); i++)
			{
				merged = merge_covertree(merged, to_merge[i]);
			}
			currentTree = merged;

			// mapRepTree.emplace(rep, new KDTreeECQ(newmat));
		}
		else
		{
			rep = vecRepTree[0].first;
			currentTree = vecRepTree[0].second;
		}
	}
}

void final_merge_on_find_arbitrary_order_eps_connectness_CoverTree_Lazy(arma::mat &data, std::unordered_map<int, std::vector<int>> &mapIndexReps, double denThres, mlpack::UnionFind &ufset)
{

	int num = data.n_cols;
	int dim = data.n_rows;

	// std::unordered_map<int, CoverTreeHere*>  mapRepTree;
	std::vector<std::pair<int, CoverTree_Lazy *>> vecRepTree;
	vecRepTree.reserve(mapIndexReps.size());
	size_t ind = 0;
#ifdef PRINT_BUILD_TREE_TIME
	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
#endif
	for (auto &pairReps : mapIndexReps)
	{
		arma::mat tmpmat(dim, pairReps.second.size());
		for (int i = 0; i < pairReps.second.size(); i++)
		{
			tmpmat.col(i) = data.col(pairReps.second[i]);
		}
		vecRepTree.push_back(std::make_pair(pairReps.first, new CoverTree_Lazy(std::move(tmpmat), 1.3)));
	}
#ifdef PRINT_BUILD_TREE_TIME
	QueryPerformanceCounter(&t2);
	std::cout << "build tree time: " << (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << std::endl;
#endif
	if (vecRepTree.size() <= 1)	return;
	std::vector<std::pair<int, CoverTree_Lazy *>> newvec;
	bool modefied = false;
	int rep = vecRepTree[0].first;
	CoverTree_Lazy *currentTree = vecRepTree[0].second;

	while (vecRepTree.size() > 1)
	{
		newvec.reserve(vecRepTree.size());
		size_t newindex = 0;
		size_t index = (modefied ? 0 : 1);
		std::vector<CoverTree_Lazy *> to_merge;

		modefied = false;

		for (; index < vecRepTree.size(); index++)
		{
			// EpsConnected_Query<mlpack::EuclideanDistance, arma::mat, mlpack::KDTree> ecq(currentTree, iter->second, denThres);
			// ecq.DoQuery();
			bool found = false;
			// eps_connectness_dual_CoverTree_Lazy(*currentTree, *vecRepTree[index].second, found);
			eps_connectness_dual_CoverTree_Lazy_cendist(*currentTree, *vecRepTree[index].second, 0, found);
			// if(ecq.result){
			if (found)
			{
				modefied = true;
				to_merge.push_back(vecRepTree[index].second);

				ufset.Union(rep, vecRepTree[index].first);
			}
			else
			{
				newvec.push_back( std::move(vecRepTree[index]));
			}
		}
		std::swap(vecRepTree, newvec);
		newvec.clear();
		if (modefied)
		{
			CoverTree_Lazy *merged = to_merge[0];
			for (size_t i = 1; i < to_merge.size(); i++)
			{
				merged = merge_covertree(merged, to_merge[i]);
			}
			currentTree = merged;

			// mapRepTree.emplace(rep, new KDTreeECQ(newmat));
		}
		else
		{
			rep = vecRepTree[0].first;
			currentTree = vecRepTree[0].second;
		}
	}
}
