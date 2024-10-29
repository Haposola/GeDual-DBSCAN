#pragma once

class KDTreeDecidedState {
public:
	KDTreeDecidedState() :dense_flag(false), onecluster_flag(false), allcore_flag(false), num_core(0) {};
	~KDTreeDecidedState() {};
	template<typename TreeType>KDTreeDecidedState(TreeType& node) :dense_flag(false), onecluster_flag(false), allcore_flag(false), num_core(0) {};

	const bool isDense() const {
		return dense_flag;
	}
	void setDense() {
		//dense is equivalent to: allcore and onecluster
		dense_flag = true;
		//make dense_flag and desc_dense_flag be oppsite
		onecluster_flag = true;
		allcore_flag = true;
	}

	const bool isAllcore() const {
		return allcore_flag;
	}
	void setAllcore() {
		allcore_flag = true;
	}
	const bool isOnecluster() const {
		return onecluster_flag;
	}

	void setOnecluster() {
		onecluster_flag = true;
	}
	void incNumCore() {
		++num_core;
	}
	void setNumCore(size_t num) {
		num_core = num;
	}
	const size_t getNumCore() const {
		return num_core;
	}
	template<typename Archive>
	void serialize(Archive& ar, const uint32_t /* version */) {
		ar(CEREAL_NVP(dense_flag));
	}
private:
	bool dense_flag;
	bool onecluster_flag;
	bool allcore_flag;
	size_t search_start;
	size_t num_core;
};

class GeDual_DBSCAN {

public:

	//typedef mlpack::BinarySpaceTree<mlpack::EuclideanDistance,mlpack::EmptyStatistic,arma::Mat<double>,mlpack::HRectBound,mlpack::MidpointSplit> KDTreeHere;
	typedef mlpack::KDTree<mlpack::EuclideanDistance, KDTreeDecidedState, arma::Mat<double>> KDTreeHere;
	GeDual_DBSCAN(
		arma::Mat<double>& data,
		//arma::vec& query,
		double dbscan_eps,
		int dbscan_minpts
		//ScoreStack<mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>*>& simulationStack,
		//mlpack::KDTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>& searchTree,
		//std::vector<size_t>& oldFromNew,
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

		rangesearch_rightstart.reserve(data.n_cols);
		for (size_t i = 0; i < data.n_cols; i++) rangesearch_rightstart.push_back(i + 1);

	}
	void DoDBSCAN(std::ostream& ofile) {
		auto t1 = ExperimentRun_TimeNow();

		//First use KDTree as grid to find core-points
		FindDenseNode();

		//Conduct stack-simulated tree-traversal 
		auto t2 = ExperimentRun_TimeNow();
		MergeSiblingDense();

		auto t3 = ExperimentRun_TimeNow();
		GeDual_Traversal();
		auto t4 = ExperimentRun_TimeNow();
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
		ofile << ExperimentRun_TimeCount_s_Seconds(t1, t2) << ", ";// find initial dense nodes time
		ofile << ExperimentRun_TimeCount_s_Seconds(t2, t3) << ", ";// merge sibling dense nodes time (using eps-connect dual)
		ofile << ExperimentRun_TimeCount_s_Seconds(t3, t4) << ", ";// Do DBSCAN traversal time
		ofile << ExperimentRun_TimeCount_s_Seconds(t4, t5) << ", ";// Assign border/noise points
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

	size_t num_initial_dense_nodes = 0;
	size_t num_initial_onecluster_nodes = 0;
	size_t num_initial_cores = 0;

	size_t num_intercase_leafdense = 0;
	size_t num_intercase_denseleaf = 0;
	size_t num_intercase_bothdense = 0;

	size_t num_epsminpts_indense = 0;
	size_t num_epsconnect_single = 0;
	size_t num_epsconnect_dual = 0;
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
				} else return p1.first->Count() > p2.first->Count();// < leads to large first
			}
			//smaller query.Begin() value will be with larger priority. Note: > in priority query will result in small first.
			else return p1.first->Begin() > p2.first->Begin();
		}
	};
	std::priority_queue <std::pair<KDTreeHere*, KDTreeHere*>, std::vector<std::pair<KDTreeHere*, KDTreeHere*>>, node_combination_cmp> dual_priqueue;

	bool* isCore;
	KDTreeHere* searchTree;
	std::vector<size_t> oldFromNew;
	std::vector<size_t> rangesearch_rightstart;
	int* numNeighbors;
	std::vector<std::vector<size_t>> epsNeighbors;//epsNeighbors stores the original index of neighbor

	mlpack::UnionFind ufset;

	void DoSetDescOnecluster(KDTreeHere* node, bool isAllcore) {
#if defined TEST_STATISTICS
		if (isAllcore) {
			node->Stat().setDense();
			num_initial_dense_nodes++;
		} else {
			node->Stat().setOnecluster();
			num_initial_onecluster_nodes++;
		}

#else 
		if (isAllcore) node->Stat().setDense();
		else node->Stat().setOnecluster();
#endif 

		if (!node->IsLeaf()) {
			DoSetDescOnecluster(node->Left(), isAllcore);
			DoSetDescOnecluster(node->Right(), isAllcore);
		}
	}
	void DoFindDense(KDTreeHere* node) {
		if (node->Bound().Diameter() <= dbscan_eps) {

			size_t nodeEnd = node->Begin() + node->Count();
			if (node->Count() >= dbscan_minpts) {
				DoSetDescOnecluster(node, true);

				for (size_t i = node->Begin(); i < nodeEnd; i++) {
					//Merge all points in node as one cluster; merge with node->Begin() is enough
					ufset.Union(oldFromNew[node->Begin()], oldFromNew[i]);
					isCore[oldFromNew[i]] = true;
#if defined TEST_STATISTICS
					num_initial_cores++;
#endif
					rangesearch_rightstart[i] = nodeEnd;
				}

			} else {
				DoSetDescOnecluster(node, false);
				for (size_t i = node->Begin(); i < nodeEnd; i++) {
					//Increase numNeighbors, minus self-neighobr
					numNeighbors[oldFromNew[i]] += (node->Count() - 1);
					for (size_t j = node->Begin(); j < nodeEnd; j++) {
						//self-neighbor here is not a big deal
						epsNeighbors[oldFromNew[i]].push_back(oldFromNew[j]);
					}
				}
			}
		} else if (!node->IsLeaf()) {
			DoFindDense(node->Left());
			DoFindDense(node->Right());
		}
	}

	void FindDenseNode() {
		DoFindDense(searchTree);
	}
	void DoMergeSiblingDense(KDTreeHere* node) {
		if (node->IsLeaf())return;
		DoMergeSiblingDense(node->Left());
		DoMergeSiblingDense(node->Right());
		if (node->Left()->Stat().isDense() && node->Right()->Stat().isDense()) {
			mlpack::RangeType<double> lr_range = node->Left()->RangeDistance(*node->Right());
#if defined  TEST_STATISTICS
			num_Score++; num_dist_compute += 2;
#endif 
			if (lr_range.Hi() <= dbscan_eps) {
#if defined  TEST_STATISTICS
				num_near_case++; 
#endif 
				ufset.Union(oldFromNew[node->Left()->Begin()], oldFromNew[node->Right()->Begin()]);
				node->Stat().setDense();
#if defined  TEST_STATISTICS
				num_initial_dense_nodes++;
#endif 
			} else if (lr_range.Lo() <= dbscan_eps) {
				bool found = false;
				BaseCase_EpsConnect_dual(node->Left(), node->Right(), found);

				if (found) {
					ufset.Union(oldFromNew[node->Left()->Begin()], oldFromNew[node->Right()->Begin()]);
					node->Stat().setDense();
#if defined  TEST_STATISTICS
					num_initial_dense_nodes++;
#endif 
				}
			}
		}
	}
	void MergeSiblingDense() {
		DoMergeSiblingDense(searchTree);
	}
	void GeDual_Traversal();
	//intermediate cases, categorized by query leaf or dense, ref leaf or dense, resulting in four InterCases
	void InterCase_BothLeaves(KDTreeHere* queryNode, KDTreeHere* refNode);
	void InterCase_QueryLeafRefDense(KDTreeHere* queryNode, KDTreeHere* refNode);
	void InterCase_QueryDenseRefLeaf(KDTreeHere* queryNode, KDTreeHere* refNode);
	void InterCase_BothDense(KDTreeHere* queryNode, KDTreeHere* refNode);
	//base cases
	inline void BaseCase_PointwiseEpsminpts(size_t query, size_t refPoint);
	void BaseCase_NoDistCompute(KDTreeHere* queryNode, KDTreeHere* refNode);
	void BaseCase_EpsMinptsInDense_single(size_t point, KDTreeHere* node, bool& finished);
	void BaseCase_EpsConnect_single(size_t point, KDTreeHere* node, bool& found);
	void BaseCase_EpsConnect_dual(KDTreeHere* queryNode, KDTreeHere* refNode, bool& found);
	//state transitions
	void PointTransite_ToCore(size_t point);

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
			<< num_leaf_nodes << ","
			<< num_BaseCase << ", "
			<< num_Score << ", "
			<< num_near_case << ", "
			<< num_far_case << ", "
			<< num_early_core << ", "
			<< num_samecluster_skip << ", "
			<< num_initial_dense_nodes << ", "
			<< num_initial_onecluster_nodes << ", "
			<< num_initial_cores << ", "
			<< num_intercase_leafdense << ", "
			<< num_intercase_denseleaf << ", "
			<< num_intercase_bothdense << ", "
			<< num_epsminpts_indense << ", "
			<< num_epsconnect_single << ", "
			<< num_epsconnect_dual << ", ";

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
#include "IntermediateCases_impl.h"
#include "BaseCases_impl.h"
#include "StateTransitions.h"