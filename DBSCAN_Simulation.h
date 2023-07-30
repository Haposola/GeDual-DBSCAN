#pragma once


class DBSCAN_StackSimulation {
	typedef mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat> KDTreeHere;
public:
	DBSCAN_StackSimulation(
		arma::mat& data,
		//arma::vec& query,
		double dbscan_eps,
		int dbscan_minpts
		//scoreStack<mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>*>& simulationStack,
		//mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>& searchTree,
		//std::vector<size_t>& oldFromNew,
	) :
		data(data),
		//query(query),
		dbscan_eps(dbscan_eps),
		dbscan_minpts(dbscan_minpts),
		ufset(mlpack::UnionFind(data.n_cols)) {
		searchTree = new KDTreeHere(data, oldFromNew);
		simulationStacks = (scoreStack< KDTreeHere*>**)
			malloc(data.n_cols * sizeof(scoreStack< KDTreeHere*>*));
		/*for (int i = 0; i < 100; i++) {
			simulationStacks[i] = new scoreStack< mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>*>(10000);

		}*/
		isCore = new bool[data.n_cols];
		for (int i = 0; i < data.n_cols; i++)isCore[i] = false;
		epsNeighbors.resize(data.n_cols);

		for (size_t i = 0; i < data.n_cols; i++) {
			epsNeighbors[i].reserve(2 * dbscan_minpts);
		}
		//ufset = mlpack::UnionFind(data.n_cols);
		/////////////oldFromNew(oldFromNew),

	}
	void DoDBSCAN(std::ostream& ofile) {
		//Identify core-points
#ifdef _WIN32
		LARGE_INTEGER t1, t2, t3, t4, t5, t6, tc;
		QueryPerformanceFrequency(&tc);
		QueryPerformanceCounter(&t1);
#elif defined __linux__
		struct timeval t1, t2, t3, t4, t5, t6, tc;
		gettimeofday(&t1, NULL);
#endif

		KDTreeGrid();//First use KDTree as grid to find core-points
#ifdef _WIN32
		QueryPerformanceCounter(&t2);
#elif defined __linux__
		gettimeofday(&t2, NULL);
#endif

		//Then conduct stack-simulated tree-traversal 
		for (size_t i = 0; i < data.n_cols; i++) {
			if (isCore[i]) continue;
			simulationStacks[i] = new scoreStack< mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>*>();//Size of stack
			EpsMinpts_Simulation(i, *simulationStacks[i]);
		}
#ifdef _WIN32
		QueryPerformanceCounter(&t3);
#elif defined __linux__
		gettimeofday(&t3, NULL);
#endif
		//Pre-cluster using union-find sets
		for (size_t i = 0; i < data.n_cols; i++) {
			if (isCore[i]) {
				for (size_t j = 0; j < epsNeighbors[i].size(); j++) {
					if (isCore[epsNeighbors[i][j]]) ufset.Union(i, epsNeighbors[i][j]);
				}
			}
		}
#ifdef _WIN32
		QueryPerformanceCounter(&t4);
#elif defined __linux__
		gettimeofday(&t4, NULL);
#endif
		//Final clustering. Restore the interrupted range searches.
		for (size_t i = 0; i < data.n_cols; i++) {
			if (isCore[i] && simulationStacks[i] != nullptr) {
				FinalClustering_Simulation(i, *simulationStacks[i]);
			}
		}
		//
		if (tofinal.size() > 0) {
			//FinalClustering_mergeonfind();
			FinalClustering_mergeonfind_simpler();
		}
#ifdef _WIN32
		QueryPerformanceCounter(&t5);
#elif defined __linux__
		gettimeofday(&t5, NULL);
#endif
		//Assigning border/noise points
		for (size_t i = 0; i < data.n_cols; i++) {
			if (isCore[i]) continue;
			for (size_t j = 0; j < epsNeighbors[i].size(); j++) {
				if (isCore[epsNeighbors[i][j]]) {
					ufset.Union(i, epsNeighbors[i][j]);
					//break;
				}
			}
		}
#ifdef _WIN32
		QueryPerformanceCounter(&t6);
#elif defined __linux__
		gettimeofday(&t6, NULL);
#endif
		//The clustering result is stored in the union-find set
		//std::cout << "core found in grid: " << numcoreGrid << ". Total cores: " << numCore << "\n";
		/*std::cout << "KDTreeGrid time: " << (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << "\n";
		std::cout << "Core-point identification time: " << (t3.QuadPart - t2.QuadPart) * 1.0 / tc.QuadPart << "\n";
		std::cout << "Pre-clustering time: " << (t4.QuadPart - t3.QuadPart) * 1.0 / tc.QuadPart << "\n";
		std::cout << "Final-clustering time: " << (t5.QuadPart - t4.QuadPart) * 1.0 / tc.QuadPart << "\n";
		std::cout << "Assinging border/noise time: " << (t6.QuadPart - t5.QuadPart) * 1.0 / tc.QuadPart << "\n";*/
#ifdef _WIN32
		ofile << (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << ", ";
		ofile << (t3.QuadPart - t2.QuadPart) * 1.0 / tc.QuadPart << ", ";
		ofile << (t4.QuadPart - t3.QuadPart) * 1.0 / tc.QuadPart << ", ";
		ofile << (t5.QuadPart - t4.QuadPart) * 1.0 / tc.QuadPart << ", ";
		ofile << (t6.QuadPart - t5.QuadPart) * 1.0 / tc.QuadPart << ", ";
		ofile << (t6.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << ", ";
#elif defined __linux__
		ofile << (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) * 1.0 / 1000000 << ", ";
		ofile << (t3.tv_sec - t2.tv_sec) + (t3.tv_usec - t2.tv_usec) * 1.0 / 1000000 << ", ";
		ofile << (t4.tv_sec - t3.tv_sec) + (t4.tv_usec - t3.tv_usec) * 1.0 / 1000000 << ", ";
		ofile << (t5.tv_sec - t4.tv_sec) + (t5.tv_usec - t4.tv_usec) * 1.0 / 1000000 << ", ";
		ofile << (t6.tv_sec - t5.tv_sec) + (t6.tv_usec - t5.tv_usec) * 1.0 / 1000000 << ", ";
		ofile << (t6.tv_sec - t1.tv_sec) + (t6.tv_usec - t1.tv_usec) * 1.0 / 1000000 << ", ";

#endif
		print_num_final_clusters(ofile);
	}

	void DOKDTreeGrid(KDTreeHere* node) {

		if (node->Bound().Diameter() <= dbscan_eps && node->Count() >= dbscan_minpts) {
			const size_t refEnd = node->Begin() + node->Count();
			for (size_t i = node->Begin(); i < refEnd; ++i) {
				ufset.Union(oldFromNew[node->Begin()], oldFromNew[i]);//epsNeighbors[oldFromNew[i]].push_back(oldFromNew[j]);
				isCore[oldFromNew[i]] = true;
			}
			tofinal.push_back(oldFromNew[node->Begin()]);

		} else if (!node->IsLeaf()) {
			DOKDTreeGrid(node->Left());
			DOKDTreeGrid(node->Right());
		}
	}
	void KDTreeGrid() {
		DOKDTreeGrid(searchTree);
	}
	//Use stack to simulate the tree traversal of EpsMinpts Query
	void EpsMinpts_Simulation(size_t query, scoreStack< mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>*>& simulationStack) {
		// Push root first. 
		// NOTE!: the score of the root is not validated. For KDTree, it is not a big problem

		simulationStack.push(searchTree);

		while (!simulationStack.isEmpty()) {
			if (isCore[query]) {//Rescore functionality
				//numCore++;
				break;
			}
			//double topScore = simulationStack.topScore();
			auto topPtr = simulationStack.topElement();
			simulationStack.pop();
			//if (topScore > dbscan_eps) continue;//Score functionality
			if (topPtr->IsLeaf()) {
				BaseCase_EpsMinpts(query, topPtr, true);//Call BaseCase with distance computation
				continue;
			}
			bool leftValid = true; double leftScore;
			mlpack::Range leftrange = topPtr->Left()->RangeDistance(data.unsafe_col(query));

			if (leftrange.Hi() <= dbscan_eps) {
				BaseCase_EpsMinpts(query, topPtr->Left(), false);//Call BaseCase without distance computation
				leftValid = false;
			} else {
				leftScore = leftrange.Lo();
				if (leftScore > dbscan_eps) leftValid = false;
			}

			// if (isCore[query]) break; It is reasonable for only EpsMinpts simulation
			// However, it may lose the information of the right child which can not be restored in the final-clustering simulation
			// So, delete the rescore check here, and only do rescore check in the beginning of the while-loop.
			bool rightValid = true; double rightScore;
			mlpack::Range rightrange = topPtr->Right()->RangeDistance(data.unsafe_col(query));
			if (rightrange.Hi() <= dbscan_eps) {
				BaseCase_EpsMinpts(query, topPtr->Right(), false);//Call BaseCase without distance computation
				rightValid = false;
			} else {
				rightScore = rightrange.Lo();
				if (rightScore > dbscan_eps) rightValid = false;
			}

			//if (isCore[query]) break;

			if (leftValid && !rightValid) simulationStack.push(topPtr->Left());
			else if (!leftValid && rightValid) simulationStack.push(topPtr->Right());
			else if (!leftValid && !rightValid) continue;
			else {
				if (leftScore > rightScore) {//push large score first, then the small score will be popped first.
					simulationStack.push(topPtr->Left());
					simulationStack.push(topPtr->Right());
				} else {
					simulationStack.push(topPtr->Right());
					simulationStack.push(topPtr->Left());
				}
			}

		}
	}

	//Restore the simulation using the recorded stack. The goal is to finish the traversal to find all eps-neighbors and join them
	void FinalClustering_Simulation(size_t query, scoreStack< mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>*>& simulationStack) {
		while (!simulationStack.isEmpty()) {
			//double topScore = simulationStack.topScore();
			auto topPtr = simulationStack.topElement();
			simulationStack.pop();
			//if (topScore > dbscan_eps) continue;
			if (topPtr->IsLeaf()) {
				BaseCase_EpsConnectivity(query, topPtr, true);
				continue;
			}
			bool leftValid = true; double leftScore;
			mlpack::Range leftrange = topPtr->Left()->RangeDistance(data.unsafe_col(query));
			if (leftrange.Hi() <= dbscan_eps) {
				BaseCase_EpsConnectivity(query, topPtr->Left(), false);//Call BaseCase without distance computation
				leftValid = false;
			} else {
				leftScore = leftrange.Lo();
				if (leftScore > dbscan_eps) leftValid = false;
			}

			// if (isCore[query]) break; It is reasonable for only EpsMinpts simulation
			// However, it may lose the information of the right child which can not be restored in the final-clustering simulation
			// So, delete the rescore check here, and only do rescore check in the beginning of the while-loop.
			bool rightValid = true; double rightScore;
			mlpack::Range rightrange = topPtr->Right()->RangeDistance(data.unsafe_col(query));
			if (rightrange.Hi() <= dbscan_eps) {
				BaseCase_EpsConnectivity(query, topPtr->Right(), false);//Call BaseCase without distance computation
				rightValid = false;
			} else {
				rightScore = rightrange.Lo();
				if (rightScore > dbscan_eps) rightValid = false;
			}
			if (leftValid && !rightValid) simulationStack.push(topPtr->Left());
			else if (!leftValid && rightValid) simulationStack.push(topPtr->Right());
			else if (!leftValid && !rightValid) continue;
			else {
				if (leftScore > rightScore) {//push large score first, then the small score will be popped first.
					simulationStack.push(topPtr->Left());
					simulationStack.push(topPtr->Right());
				} else {
					simulationStack.push(topPtr->Right());
					simulationStack.push(topPtr->Left());
				}
			}

		}

	}
	bool core(size_t point) {
		return isCore[point];
	}
private:

	//hold references to these things

	arma::mat& data;
	double dbscan_eps;
	int dbscan_minpts;
	scoreStack< mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>*>** simulationStacks;
	mlpack::EuclideanDistance metric;
	bool* isCore;

	mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>* searchTree;
	std::vector<size_t> oldFromNew;

	//private things
	std::vector<std::vector<size_t>> epsNeighbors;//epsNeighbors stores the original index of neighbor
	std::vector<size_t> tofinal;

	mlpack::UnionFind ufset;

	void BaseCase_EpsMinpts(size_t query, mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>* node, bool computeDistance) {
		const size_t refEnd = node->Begin() + node->Count();
		for (size_t i = node->Begin(); i < refEnd; ++i) {
			if (computeDistance) {
				//double dist = mlpack::EuclideanDistance::Evaluate(data.unsafe_col(query), data.unsafe_col(oldFromNew[i]));
				double dist = metric.Evaluate(data.unsafe_col(query), data.unsafe_col(oldFromNew[i]));
				if (dist <= dbscan_eps) {
					epsNeighbors[query].push_back(oldFromNew[i]);
				}

			} else {
				epsNeighbors[query].push_back(oldFromNew[i]);
			}
		}
		if (!isCore[query]) {
			isCore[query] = (epsNeighbors[query].size() >= dbscan_minpts);
		}
	}
	void BaseCase_EpsConnectivity(size_t query, mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>* node, bool computeDistance) {
		const size_t refEnd = node->Begin() + node->Count();
		for (size_t i = node->Begin(); i < refEnd; ++i) {
			if (!isCore[oldFromNew[i]]) continue;
			if (ufset.Find(query) == ufset.Find(oldFromNew[i])) continue;
			if (computeDistance) {
				double dist = metric.Evaluate(data.unsafe_col(query), data.unsafe_col(oldFromNew[i]));
				if (dist <= dbscan_eps) {
					//epsNeighbors[query].push_back(oldFromNew[i]);No need to update the neighbors vector.
					ufset.Union(query, oldFromNew[i]);
				}
			} else {
				ufset.Union(query, oldFromNew[i]);
			}
		}

	}
	void FinalClustering_mergeonfind_simpler() {
		//std::unordered_set<size_t> undetermined;// The representatives of the clusters identified in KDTreeGrid. If it is not empty, the MergeOnFind final-clustering process will be invoked.
		//for (size_t i = 0; i < tofinal.size(); i++) {
		//	undetermined.insert(ufset.Find(tofinal[i]));
		//}
		
		std::unordered_map<size_t, std::vector<size_t>> mapIndexReps;
		//std::map<size_t, std::vector<size_t>> mapIndexReps;
		for (size_t i = 0; i < tofinal.size();i++) {
			mapIndexReps.emplace(ufset.Find(tofinal[i]), std::vector<size_t>());
		}
		
		for (size_t i = 0; i < data.n_cols; i++) {
			if (!isCore[i])continue;
			size_t rep = ufset.Find(i);
			
			if (mapIndexReps.contains(rep)) {
				mapIndexReps.at(rep).push_back(i);
			}
		}
		std::vector <std::pair<size_t, CoverTree_Lazy*>> vecRepTree;
		for (auto iter = mapIndexReps.begin(); iter != mapIndexReps.end();iter++) {
			arma::mat* tmpmat = new arma::mat(data.n_rows, iter->second.size());
			for (size_t i = 0; i < iter->second.size(); i++) {
				tmpmat->col(i) = data.col(iter->second[i]);
			}
			vecRepTree.push_back(std::make_pair(iter->first, new CoverTree_Lazy(*tmpmat)));
		}
		std::vector<std::pair<size_t, CoverTree_Lazy*>> newvec;
		std::vector<CoverTree_Lazy*> to_merge;
		bool modefied = false; bool found;
		int rep = vecRepTree[0].first;
		CoverTree_Lazy* currentTree = vecRepTree[0].second;

		while (vecRepTree.size() > 1) {
			newvec.reserve(vecRepTree.size());
			size_t index = (modefied ? 0 : 1);
			modefied = false;

			for (; index < vecRepTree.size(); index++) {
				
				 found= false;
			
				eps_connectness(*currentTree, *vecRepTree[index].second, found);
				
				//if(ecq.result){
				if (found) {
					modefied = true;
					to_merge.push_back(vecRepTree[index].second);
					ufset.Union(rep, vecRepTree[index].first);
				} else {
					newvec.push_back(std::move(vecRepTree[index]));
					//newvec.push_back(vecRepTree[index]);
				}
			}
			std::swap(vecRepTree, newvec);
			newvec.clear();
			if (modefied) {
				CoverTree_Lazy* merged = to_merge[0];
				for (size_t i = 1; i < to_merge.size(); i++) {
					merged = merge_covertree(merged, to_merge[i]);
				}
				currentTree = merged;
				to_merge.clear();
			} else {
				rep = vecRepTree[0].first;
				currentTree = vecRepTree[0].second;
			}


		}
	}
	void FinalClustering_mergeonfind() {
		std::set<size_t> undetermined;// The representatives of the clusters identified in KDTreeGrid. If it is not empty, the MergeOnFind final-clustering process will be invoked.
		for (size_t i = 0; i < tofinal.size(); i++) {
			undetermined.insert(ufset.Find(tofinal[i]));
		}
		size_t num = data.n_cols;
		size_t dim = data.n_rows;
		std::unordered_map<size_t, std::vector<size_t>> mapIndexReps;
		for (size_t i = 0; i < num; i++) {
			if (!isCore[i])continue;
			size_t rep = ufset.Find(i);
			if (mapIndexReps.find(rep) == mapIndexReps.end()) {
				std::vector<size_t> tmp; tmp.push_back(i);
				mapIndexReps.emplace(rep, tmp);
			} else {
				mapIndexReps.at(rep).push_back(i);
			}
		}
		if (mapIndexReps.size() <= 1) return;

		std::vector <std::pair<size_t, CoverTree_Lazy*>> vecRepTree;
		std::unordered_map <size_t, CoverTree_Lazy*> undeterminedRepTree;
		vecRepTree.reserve(mapIndexReps.size());


		for (auto& pairReps : mapIndexReps) {
			arma::mat* tmpmat = new arma::mat(dim, pairReps.second.size());
			for (size_t i = 0; i < pairReps.second.size(); i++) {
				tmpmat->col(i) = data.col(pairReps.second[i]);
			}
			CoverTree_Lazy* tmpctptr = new CoverTree_Lazy(*tmpmat);
			vecRepTree.push_back(std::make_pair(pairReps.first, tmpctptr));
			if (undetermined.contains(pairReps.first)) {
				undeterminedRepTree.emplace(pairReps.first, tmpctptr);
			} 
		}


		//Goal
		std::vector <std::pair<size_t, CoverTree_Lazy*>> newvec;
		bool modefied;	size_t rep;
		CoverTree_Lazy* currentTree;

		while (undeterminedRepTree.size() > 0) {
			rep = undeterminedRepTree.begin()->first;
			currentTree = undeterminedRepTree.begin()->second;

			undeterminedRepTree.erase(rep);

			//newvec.reserve(vecRepTree.size());

			do {
				std::vector<CoverTree_Lazy*> to_merge;
				modefied = false;
				for (size_t index = 0; index < vecRepTree.size(); index++) {
					if (rep == vecRepTree[index].first)
						continue;//correct. this rep is not inserted to new vec, namely it is deleted
					bool epsConnected = false;
					eps_connectness(*currentTree, *vecRepTree[index].second, epsConnected);

					if (epsConnected) {
						modefied = true;
						to_merge.push_back(vecRepTree[index].second);

						ufset.Union(rep, vecRepTree[index].first);
						if (undeterminedRepTree.contains(vecRepTree[index].first)) {
							undeterminedRepTree.erase(vecRepTree[index].first);
						} else {
							std::cout << 1;
						}
					} else {
						newvec.push_back(std::move(vecRepTree[index]));
					}
				}//eps_connectivity search loop

				std::swap(vecRepTree, newvec);
				newvec.clear();
				if (modefied) {
					CoverTree_Lazy* merged = to_merge[0];
					for (size_t i = 1; i < to_merge.size(); i++) {
						merged = merge_covertree(merged, to_merge[i]);
					}
					currentTree = merged;
				}
			} while (modefied);
		}
	}
	inline void eps_connectness(CoverTree_Lazy& queryRoot,
								CoverTree_Lazy& refRoot,
								bool& found) {
		if (mlpack::EuclideanDistance::Evaluate(queryRoot.Dataset().col(queryRoot.Point()), refRoot.Dataset().col(refRoot.Point())) <= dbscan_eps) {
			found = true;
			return;
		} else {
			eps_connectness_dual_CoverTree_Lazy(queryRoot, refRoot, found);
		}
	}

	void eps_connectness_dual_CoverTree_Lazy(CoverTree_Lazy& queryNode,
											 CoverTree_Lazy& refNode,
											 bool& found) {
		//No self-children. Avoid related conditions.
		if (queryNode.IsLeaf() && refNode.IsLeaf())
			return;
		if (!queryNode.IsExpanded())queryNode.Expand();
		if (!refNode.IsExpanded()) refNode.Expand();
		//It is assumed that BaseCase(queryNode, refNode) has been calculated.
		//Then we only need to consider the children node combinations
		//And that is when BaseCase is invoked on the children node combination. Invariant holds.
		if (!queryNode.IsLeaf() && !refNode.IsLeaf()) {
			
			std::vector<std::tuple<size_t, size_t, double>> score_list;
			//score_list.reserve((numq) * (numr));
			for (size_t i = 0; i < queryNode.NumChildren(); i++) {
				auto&  queryChild = queryNode.Child(i);
				for (size_t j = 0; j < refNode.NumChildren(); j++) {
					auto& refChild = refNode.Child(j);
					double dist = mlpack::EuclideanDistance::Evaluate(queryChild.Dataset().col(queryChild.Point()), refChild.Dataset().col(refChild.Point()));
					if (dist <= dbscan_eps) {
						found = true;
						return;
					}
					// Seems no need to check the upper bound, since if the upper bound is within denThres, the base_case above should return true;
					// if (dist + queryNode.Child(i).FurthestDescendantDistance() + refNode.Child(j).FurthestDescendantDistance() <= denThres) {
					// found = true; return;
					//}
					double score = dist - queryChild.FurthestDescendantDistance() - refChild.FurthestDescendantDistance();
					//double score = dist - pow(queryNode.Base(), queryNode.Child(i).Scale() + 2) - pow(refNode.Base(), refNode.Child(j).Scale() + 2);
					if (score > dbscan_eps) continue;
					else score_list.push_back(std::make_tuple(i, j, score));
				}
			}
			std::sort(score_list.begin(), score_list.end(),
					  [](const std::tuple<size_t, size_t, double>& p1, const std::tuple<size_t, size_t, double>& p2) {
						  return std::get<2>(p1) == std::get<2>(p2) ? std::get<0>(p1) < std::get<0>(p2) : std::get<2>(p1) < std::get<2>(p2);
					  });
			for (size_t i = 0; i < score_list.size(); i++) {
				eps_connectness_dual_CoverTree_Lazy(queryNode.Child(std::get<0>(score_list[i])), refNode.Child(std::get<1>(score_list[i])), found);
				if (found)
					return;
			}
			//return;
		} else { // one leaf, another not
			auto& leafOne = queryNode.IsLeaf() ? queryNode : refNode;
			auto& nonLeafOne = queryNode.IsLeaf() ? refNode : queryNode;
			
			std::vector<std::pair<size_t, double>> score_list;
			//score_list.reserve(numc);

			for (size_t i = 0; i < nonLeafOne.NumChildren(); i++) {
				auto& nlChild = nonLeafOne.Child(i);
				double dist = mlpack::EuclideanDistance::Evaluate(nlChild.Dataset().col(nlChild.Point()), leafOne.Dataset().col(leafOne.Point()));
				if (dist <= dbscan_eps) {
					found = true;
					return;
				}
				// No need to check the upper bound
				// if (dist + nonLeafOne.Child(i).FurthestDescendantDistance() <= denThres) {//leafOne.FurthestDescendantDistance()=0
				// found = true; return;
				//}
				
				double score = dist - nlChild.FurthestDescendantDistance(); // leafOne.FurthestDescendantDistance()=0
				//double score = dist - pow(nonLeafOne.Base(), nonLeafOne.Child(i).Scale() + 2);
				if (score > dbscan_eps)	continue;
				else score_list.push_back(std::make_pair(i, score));
			}
			std::sort(score_list.begin(), score_list.end(),
					  [](const std::pair<size_t, double>& p1, const std::pair<size_t, double>& p2) {
						  return p2.second == p1.second ? p1.first < p2.first : p1.second < p2.second;
					  });
			for (size_t i = 0; i < score_list.size(); i++) {
				eps_connectness_dual_CoverTree_Lazy(leafOne, nonLeafOne.Child(score_list[i].first), found);
				if (found)
					return;
			}
		}
	}

	
	void print_num_final_clusters(std::ostream& ofile) {
		int numcore = 0;
		for (size_t i = 0; i < data.n_cols; i++) {
			if (isCore[i]) numcore++;
		}
		ofile << numcore << ", ";
		std::unordered_map<int, std::vector<int>> mapIndexReps;
		int final_clusters = 0;
		for (int i = 0; i < data.n_cols; i++) {
			int rep = ufset.Find(i);
			if (!mapIndexReps.contains(rep)) {
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
		ofile << final_clusters << "\n";
		ofile.flush();
	}
};