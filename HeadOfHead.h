#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#ifdef _WIN32
#include <Windows.h>
#include <mlpack/mlpack.hpp>
#else
#include <mlpack.hpp>
#endif


double denThres;

//#define PRINT_TIME
#define FILOUT_TIME

//#define PRINT_NUM_PRE_CLUSTERS

//#define TEST_SYNTHETIC
#define TEST_FILE
//#define TEST_CMD

#define PRINT_NUM_FINAL_CLUSTERS

//#define FILEOUT_NUM_DIST_COMPUTE
//#define PRINT_NUM_DIST_COMPUTE
//#define GREEDY_PRE_CLUSTER

//Final cluster methods. 
//choose one from simpletraverse : mergeonfind
//choose to use largefirst or not

//#define AKNN_FINAL_CLUSTER
#define EPS_DETERMIN_FINAL_CLUSTER
//#define FINAL_SIMPLE_TRAVERSE
#define FINAL_MERGE_ON_FIND
//#define LARGE_TREE_FIRST_VISIT

#if (defined FILEOUT_NUM_DIST_COMPUTE)||(defined PRINT_NUM_DIST_COMPUTE)
typedef Counted_Eucdist EucDistHere;
#else
typedef mlpack::EuclideanDistance EucDistHere;
#endif
template<typename T>
T& createNullRef() {return *static_cast<T*>(nullptr);}
struct pair2cmp {
	bool operator()(const std::pair<size_t, double>& p1, const std::pair<size_t, double>& p2) {
		return p1.second == p2.second ? p1.first < p2.first : p1.second < p2.second;
	}
};
#include "CoverTree_Lazy.h"
#include "Merge_CoverTree.h"

#include "EpsMinpts_Statistics.h"
#include "EpsMinpts_Rules.h"
#include "EpsMinpts_Query.h"



#include "EpsConnected_DFS.h"
#include "Finalcluser_AllkNN.h"
#include "Finalcluster_Methods.h"


#include "DBSCAN_DualTraversals.h"
//#include "test_codes.h"

