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



#define PRINT_NUM_PRE_CLUSTERS

//#define TEST_SYNTHETIC
#define TEST_FILE
//#define TEST_CMD
//#define RUNEXP

template<typename T>
T& createNullRef() {return *static_cast<T*>(nullptr);}

#include "scoreStack.h"
#include "CoverTree_Lazy.h"
#include "Merge_CoverTree.h"


#include "DBSCAN_Simulation.h"
#include "test_codes.h"

