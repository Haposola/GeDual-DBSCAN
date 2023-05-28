#pragma once
template<typename MetricType = mlpack::EuclideanDistance,
	typename MatType = arma::mat,
	template<typename TreeMetricType = MetricType,
	typename TreeStatType = EpsMinptsStat,
	typename TreeMatType = MatType> class TreeType = mlpack::KDTree>
class EpsMinpts_Query {//Similar with classes such as RangeSearch
	typedef TreeType<MetricType, EpsMinptsStat, MatType> Tree;

public:
	EpsMinpts_Query(MatType referenceSet, double dbscan_eps, int dbscan_minpts,
					const bool naive = false,
					const bool singleMode = false,
					const MetricType metric = MetricType());

	EpsMinpts_Query(Tree* referenceTree, double dbscan_eps, int dbscan_minpts,
					const bool singleMode = false,
					const MetricType metric = MetricType());
	//The other constructors such as copy, move, and operator= constructors are not needed here.
	~EpsMinpts_Query();
private:
	double dbscan_eps;
	int dbscan_minpts;


public:
	//Conduct dual tree search on a different query set
	//This function will be used in the final-cluster stage of DBSCAN.
	void DoQuery(const MatType& querySet,
				std::vector<std::vector<size_t>>& neighbors,
				std::vector<std::vector<double>>& distances,
				 std::vector<bool>& isCore);

	void DoQuery(Tree* queryTree,
				 std::vector<std::vector<size_t>>& neighbors,
				 std::vector<std::vector<double>>& distances,
				 std::vector<bool>& isCore);

	//Conduct self dual tree search on the reference data.
	//This function will be used in the pre-cluster stage of DBSCAN.
	void DoQuery(std::vector<std::vector<size_t>>& neighbors,
				 std::vector<std::vector<double>>& distances,
				 std::vector<bool>& isCore);

private:
	//! Mappings to old reference indices (used when this object builds trees).
	std::vector<size_t> oldFromNewReferences;
	//! Reference tree.
	Tree* referenceTree;
	//! Reference set (data should be accessed using this).  In some situations we
	//! may be the owner of this.
	const MatType* referenceSet;

	//! If true, this object is responsible for deleting the trees.
	bool treeOwner;

	//! If true, O(n^2) naive computation is used.
	bool naive;
	//! If true, single-tree computation is used.
	bool singleMode;

	//! Instantiated distance metric.
	MetricType metric;

	//! The total number of base cases during the last search.
	//size_t baseCases;
	//! The total number of scores during the last search.
	//size_t scores;

};


#include "EpsMinpts_Query_impl.h"