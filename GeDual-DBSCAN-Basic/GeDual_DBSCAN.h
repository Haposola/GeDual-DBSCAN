#pragma once


class DBSCAN_GeDual {

public:
	typedef mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::Mat<double>> KDTreeHere;
	DBSCAN_GeDual(
		arma::Mat<double>& data,
		//arma::vec& query,
		double dbscan_eps,
		int dbscan_minpts
	) :
		dataset(data),
		//query(query),
		dbscan_eps(dbscan_eps),
		dbscan_minpts(dbscan_minpts),
		ufset(mlpack::UnionFind(dataset.n_cols)) {

		searchTree = new KDTreeHere(dataset, oldFromNew);

		isCore = new bool[data.n_cols];
		for (int i = 0; i < data.n_cols; i++)isCore[i] = false;

		numNeighbors = new int[data.n_cols];
		for (int i = 0; i < data.n_cols; i++) numNeighbors[i] = 1;

		epsNeighbors.resize(data.n_cols);
		for (int i = 0; i < data.n_cols; i++) epsNeighbors[i].reserve(dbscan_minpts);

	}
	void DoDBSCAN(std::ostream& ofile) {
		auto t1 = ExperimentRun_TimeNow();

		GeDual_Traversal();

		//Assigning border/noise points
		for (size_t i = 0; i < dataset.n_cols; i++) {
			if (isCore[i]) continue;
			for (size_t neighbor : epsNeighbors[i]) {
				if (isCore[neighbor]) {
					ufset.Union(i, neighbor);
					break;
				}
			}
		}
		auto t5 = ExperimentRun_TimeNow();
		ofile << dataset.n_cols << ", " << dataset.n_rows << ", " << dbscan_eps << ", " << dbscan_minpts << ", ";

		ofile << ExperimentRun_TimeCount_s_Seconds(t1, t5) << ", ";// total time

		print_statistics(ofile);
	}
#if defined TEST_STATISTICS
private:
	size_t num_dist_compute = 0;			// #DistCompute= #Pointwise + 2 * #Score
	size_t num_pointwise_dist_compute = 0;	// 
	size_t num_leaf_nodes = 0;				// #LeafNodes of KDTree
	size_t num_BaseCase = 0;				// #BaseCase
	size_t num_Score = 0;					// #Score= #NearCase + #FarCase + #NormalCase
	size_t num_near_case = 0;				// cases where Dmax(N1,N2)\le epsilon
	size_t num_far_case = 0;				// cases where Dmin(N1,N2)\ge epsilon
	size_t num_early_core = 0;				// # of core points identified before it is visited as a query point
	size_t num_samecluster_skip = 0;		// # of skipped distance computation due to ufset.Find(p1)== ufset.Find(p2)
#endif

private:
	arma::Mat<double>& dataset;
	double dbscan_eps;
	int dbscan_minpts;

	struct node_combination_cmp {
		bool operator()(const std::pair<KDTreeHere*, KDTreeHere*>& p1, const std::pair<KDTreeHere*, KDTreeHere*>& p2) const {

			//if same query.Begin(), larger query.Count() value will be larger priority
			if (p1.first->Begin() == p2.first->Begin()) {
				//if same count(the query node is the same for the two pairs), 
				if (p1.first->Count() == p2.first->Count()) {
					//refNode with smaller begin will be large priority
					return p1.second->Begin() > p2.second->Begin();
				} else return p1.first->Count() < p2.first->Count();// < leads to large first
			}
			//smaller query.Begin() value will be with larger priority. Note: > in priority query will result in small first.
			else return p1.first->Begin() > p2.first->Begin();
		}
	};
	std::priority_queue <std::pair<KDTreeHere*, KDTreeHere*>, std::vector<std::pair<KDTreeHere*, KDTreeHere*>>, node_combination_cmp> dual_priqueue;

	bool* isCore;
	KDTreeHere* searchTree;
	std::vector<size_t> oldFromNew;
	int* numNeighbors;
	std::vector<std::vector<size_t>> epsNeighbors;//epsNeighbors stores the original index of neighbor

	mlpack::UnionFind ufset;


	inline void pointwise_basecase_epsminpts(size_t query, size_t refPoint, bool computeDistance) {
#if defined TEST_STATISTICS
		if (computeDistance) {
			num_pointwise_dist_compute++;
			num_dist_compute++;
		}
#endif

		if (isCore[query]) {
			//skip intra-cluster distance computation
			if (ufset.Find(query) == ufset.Find(refPoint)) {
#if defined TEST_STATISTICS
				num_samecluster_skip++;
#endif
				return;
			}
			if (!computeDistance || (computeDistance && mlpack::EuclideanDistance::Evaluate(dataset.unsafe_col(query), dataset.unsafe_col(refPoint)) <= dbscan_eps)) {
				if (isCore[refPoint]) ufset.Union(query, refPoint);
				else {
					//Invariant: After a point is changed into core-point, the neighbors will not matter
					//If a non-core neighbor point is visited, the neighbor of non-core will be updated, 
					//and the merge operation will be executed on the time the non-core point is change into core-point.
					//If a core-point neighbor is visited, the merge operation will be conducted immediately. 
					epsNeighbors[refPoint].push_back(query);
					//if refPoints is changed into a core-point, conduct not-yet-done merge.
					if (++numNeighbors[refPoint] >= dbscan_minpts) {
						isCore[refPoint] = true;
#if defined TEST_STATISTICS
						num_early_core++;
#endif
						for (size_t refnbr : epsNeighbors[refPoint]) {
							if (isCore[refnbr]) ufset.Union(refPoint, refnbr);
						}
					}
				}
			}
		} else {
			if (!computeDistance || (computeDistance && mlpack::EuclideanDistance::Evaluate(dataset.unsafe_col(query), dataset.unsafe_col(refPoint)) <= dbscan_eps)) {
				epsNeighbors[query].push_back(refPoint);
				if (++numNeighbors[query] >= dbscan_minpts) {
					isCore[query] = true;
					for (size_t queryNbr : epsNeighbors[query]) {
						if (isCore[queryNbr]) ufset.Union(query, queryNbr);
					}
				}
				if (!isCore[refPoint]) {
					epsNeighbors[refPoint].push_back(query);
					if (++numNeighbors[refPoint] >= dbscan_minpts) {
						isCore[refPoint] = true;
#if defined TEST_STATISTICS
						num_early_core++;
#endif
						for (size_t refNbr : epsNeighbors[refPoint]) {
							if (isCore[refNbr]) ufset.Union(refPoint, refNbr);
						}
					}
				}
			}
		}
	}
	void BaseCase_EpsMinpts(KDTreeHere* queryNode, KDTreeHere* refNode, bool computeDistance) {
#if defined TEST_STATISTICS
		num_BaseCase++;
#endif
		size_t queryEnd = queryNode->Begin() + queryNode->Count();
		size_t refEnd = refNode->Begin() + refNode->Count();

		for (size_t i = queryNode->Begin(); i < queryEnd; i++) {
			for (size_t j = std::max(i + 1, refNode->Begin()); j < refEnd; j++) {
				pointwise_basecase_epsminpts(oldFromNew[i], oldFromNew[j], computeDistance);
			}
		}

	}
	void GeDual_Traversal();
	void count_leaf_nodes(KDTreeHere* node) {
		if (node->IsLeaf())num_leaf_nodes++;
		else {
			count_leaf_nodes(node->Left());
			count_leaf_nodes(node->Right());

		}
	}
	void print_statistics(std::ostream& ofile) {
		int numcore = 0;
		for (size_t i = 0; i < dataset.n_cols; i++) {
			if (isCore[i]) numcore++;
		}
		
		ofile << numcore << ", ";
		
#if defined  TEST_STATISTICS
		count_leaf_nodes(searchTree);
		ofile << num_dist_compute << ", "
			<< num_pointwise_dist_compute << ", "
			<< num_leaf_nodes<<","
			<< num_BaseCase << ", "
			<< num_Score << ", "
			<< num_near_case << ", "
			<< num_far_case << ", "
			<< num_early_core << ", "
			<< num_samecluster_skip << ", ";
#endif 
		std::unordered_map<int, std::vector<int>> mapIndexReps;
		int final_clusters = 0;
		for (int i = 0; i < dataset.n_cols; i++) {
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
#include "GeDual_Traversal_impl.h"
