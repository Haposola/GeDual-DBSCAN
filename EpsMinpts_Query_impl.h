
template<typename MetricType, typename MatType,
	template<typename TreeMetricType,
	typename TreeStatType,
	typename TreeMatType> class TreeType >
EpsMinpts_Query<MetricType, MatType, TreeType>::EpsMinpts_Query(
	MatType referenceSet,
	double dbscan_eps, int dbscan_minpts,
	const bool naive, const bool singleMode,
	const MetricType metric) :
	dbscan_eps(dbscan_eps), dbscan_minpts(dbscan_minpts),
	referenceTree(naive ? NULL : mlpack::BuildTree<Tree>(std::move(referenceSet),oldFromNewReferences)),
	referenceSet(naive ? new MatType(std::move(referenceSet)) : &referenceTree->Dataset()),
	treeOwner(!naive),
	naive(naive),
	singleMode(!naive && singleMode),
	metric(metric) {
	//Nothing to do
}
template<typename MetricType, typename MatType,
	template<typename TreeMetricType,
	typename TreeStatType,
	typename TreeMatType> class TreeType >
EpsMinpts_Query<MetricType, MatType, TreeType>::EpsMinpts_Query(
	Tree* referenceTree, double dbscan_eps, int dbscan_minpts,
	const bool singleMode,
	const MetricType metric) :
	referenceTree(referenceTree),
	referenceSet(&referenceTree->Dataset()),
	dbscan_eps(dbscan_eps), dbscan_minpts(dbscan_minpts),
	treeOwner(false),
	naive(false),
	singleMode(singleMode),
	metric(metric)
{
	// Nothing else to initialize.
}
template<typename MetricType, typename MatType,
	template<typename TreeMetricType,
	typename TreeStatType,
	typename TreeMatType> class TreeType >
EpsMinpts_Query<MetricType, MatType, TreeType>::~EpsMinpts_Query() {
	if (treeOwner && referenceTree)
		delete referenceTree;
	if (naive && referenceSet)
		delete referenceSet;
}

template<typename MetricType, typename MatType,
	template<typename TreeMetricType,
	typename TreeStatType,
	typename TreeMatType> class TreeType >
void EpsMinpts_Query<MetricType, MatType, TreeType>::DoQuery(
	const MatType& querySet,
	std::vector<std::vector<size_t>>& neighbors,
	std::vector<std::vector<double>>& distances,
	std::vector<bool>& isCore) {
	mlpack::util::CheckSameDimensionality(querySet, *referenceSet,
								  "EpsMinpts_Query::DoQuery()", "query set");
	if (referenceSet->n_cols == 0)
		return;

	// This will hold mappings for query points, if necessary.
	std::vector<size_t> oldFromNewQueries;
	// If we have built the trees ourselves, then we will have to map all the
 // indices back to their original indices when this computation is finished.
 // To avoid extra copies, we will store the unmapped neighbors and distances
 // in a separate object.
	std::vector<std::vector<size_t>>* neighborPtr = &neighbors;
	std::vector<std::vector<double>>* distancePtr = &distances;
	std::vector<bool>* corePtr = &isCore;
	if (mlpack::TreeTraits<Tree>::RearrangesDataset) {
		// Query indices only need to be mapped if we are building the query tree
		// ourselves.
		if (!singleMode && !naive) {
			distancePtr = new std::vector<std::vector<double>>;
			neighborPtr = new std::vector<std::vector<size_t>>;
			corePtr = new std::vector<bool>;
		}

		// Reference indices only need to be mapped if we built the reference tree
		// ourselves.
		else if (treeOwner)
			neighborPtr = new std::vector<std::vector<size_t>>;
	}
	// Resize each vector.
	neighborPtr->clear(); // Just in case there was anything in it.
	neighborPtr->resize(querySet.n_cols);
	distancePtr->clear();
	distancePtr->resize(querySet.n_cols);
	corePtr->clear();
	corePtr->resize(querySet.n_cols);

	typedef EpsMinpts_Rules<MetricType, Tree> RuleType;
	 // Dual-tree recursion.

		// Build the query tree.
		Tree* queryTree = mlpack::BuildTree<Tree>(querySet, oldFromNewQueries);

		// Create the traverser.
		RuleType rules(*referenceSet, queryTree->Dataset(), dbscan_eps, dbscan_minpts,
					   *neighborPtr, *distancePtr, *corePtr, metric);
		typename Tree::template DualTreeTraverser<RuleType> traverser(rules);

		traverser.Traverse(*queryTree, *referenceTree);

		//baseCases += rules.BaseCases();
		//scores += rules.Scores();

		// Clean up tree memory.
		delete queryTree;
	


	// Map points back to original indices, if necessary.
	// The isCore vector is only related to the query set.
	// Thus, it is only mapped when the query set should be mapped.
	if (mlpack::TreeTraits<Tree>::RearrangesDataset) {
		if (!singleMode && !naive && treeOwner) {
			// We must map both query and reference indices.
			neighbors.clear();
			neighbors.resize(querySet.n_cols);
			distances.clear();
			distances.resize(querySet.n_cols);
			isCore.clear();
			isCore.resize(querySet.n_cols);
			for (size_t i = 0; i < distances.size(); ++i) {
				// Map distances (copy a column).
				const size_t queryMapping = oldFromNewQueries[i];
				distances[queryMapping] = (*distancePtr)[i];
				isCore[queryMapping] = (*corePtr)[i];
				// Copy each neighbor individually, because we need to map it.
				neighbors[queryMapping].resize(distances[queryMapping].size());
				for (size_t j = 0; j < distances[queryMapping].size(); ++j)
					neighbors[queryMapping][j] =
					oldFromNewReferences[(*neighborPtr)[i][j]];

			}

			// Finished with temporary objects.
			delete neighborPtr;
			delete distancePtr;
			delete corePtr;
		} else if (!singleMode && !naive) {
			// We must map query indices only.
			neighbors.clear();
			neighbors.resize(querySet.n_cols);
			distances.clear();
			distances.resize(querySet.n_cols);
			isCore.clear();
			isCore.resize(querySet.n_cols);
			for (size_t i = 0; i < distances.size(); ++i) {
				// Map distances and neighbors (copy a column).
				const size_t queryMapping = oldFromNewQueries[i];
				distances[queryMapping] = (*distancePtr)[i];
				neighbors[queryMapping] = (*neighborPtr)[i];
				isCore[queryMapping] = (*corePtr)[i];
			}

			// Finished with temporary objects.
			delete neighborPtr;
			delete distancePtr;
			delete corePtr;
		} else if (treeOwner) {
			// We must map reference indices only.
			neighbors.clear();
			neighbors.resize(querySet.n_cols);

			for (size_t i = 0; i < neighbors.size(); ++i) {
				neighbors[i].resize((*neighborPtr)[i].size());
				for (size_t j = 0; j < neighbors[i].size(); ++j)
					neighbors[i][j] = oldFromNewReferences[(*neighborPtr)[i][j]];
			}

			// Finished with temporary object.
			delete neighborPtr;
		}
	}
}

template<typename MetricType, typename MatType,
	template<typename TreeMetricType,
	typename TreeStatType,
	typename TreeMatType> class TreeType >
void EpsMinpts_Query<MetricType, MatType, TreeType>::DoQuery(
	Tree* queryTree,
	std::vector<std::vector<size_t>>& neighbors,
	std::vector<std::vector<double>>& distances,
	std::vector<bool>& isCore) {

	// If there are no points, there is no search to be done.
	if (referenceSet->n_cols == 0)
		return;

	// Get a reference to the query set.
	const MatType& querySet = queryTree->Dataset();

	// Make sure we are in dual-tree mode.
	if (singleMode || naive)
		throw std::invalid_argument("cannot call RangeSearch::Search() with a "
									"query tree when naive or singleMode are set to true");

	// We won't need to map query indices, but will we need to map distances?
	std::vector<std::vector<size_t>>* neighborPtr = &neighbors;
	if (treeOwner && mlpack::TreeTraits<Tree>::RearrangesDataset)
		neighborPtr = new std::vector<std::vector<size_t>>;
	// Resize each vector.
	neighborPtr->clear(); // Just in case there was anything in it.
	neighborPtr->resize(querySet.n_cols);
	distances.clear();
	distances.resize(querySet.n_cols);
	isCore.clear();
	isCore.clear(querySet.n_cols);

	typedef EpsMinpts_Rules<MetricType, Tree> RuleType;
	//typedef EpsMinpts_Rules RuleType;
	RuleType rules(*referenceSet, queryTree->Dataset(), dbscan_eps, dbscan_minpts, *neighborPtr,
				   distances, isCore, metric /* don't return the query in the results */);
	typename Tree::template DualTreeTraverser<RuleType> traverser(rules);

	traverser.Traverse(*queryTree, *referenceTree);

	// Do we need to map indices?
	if (treeOwner && mlpack::TreeTraits<Tree>::RearrangesDataset) {
		// We must map reference indices only.
		neighbors.clear();
		neighbors.resize(querySet.n_cols);

		for (size_t i = 0; i < neighbors.size(); ++i) {
			neighbors[i].resize((*neighborPtr)[i].size());
			for (size_t j = 0; j < neighbors[i].size(); ++j)
				neighbors[i][j] = oldFromNewReferences[(*neighborPtr)[i][j]];
		}

		// Finished with temporary object.
		delete neighborPtr;
	}

}


template<typename MetricType, typename MatType,
	template<typename TreeMetricType,
	typename TreeStatType,
	typename TreeMatType> class TreeType >
void EpsMinpts_Query<MetricType, MatType, TreeType>::DoQuery(std::vector<std::vector<size_t>>& neighbors,
															 std::vector<std::vector<double>>& distances,
															 std::vector<bool>& isCore) {
	// If there are no points, there is no search to be done.
	if (referenceSet->n_cols == 0)	return;

	// Here, we will use the query set as the reference set.
	std::vector<std::vector<size_t>>* neighborPtr = &neighbors;
	std::vector<std::vector<double>>* distancePtr = &distances;
	std::vector<bool>* corePtr = &isCore;

	if (mlpack::TreeTraits<Tree>::RearrangesDataset && treeOwner) {
		// We will always need to rearrange in this case.
		distancePtr = new std::vector<std::vector<double>>;
		neighborPtr = new std::vector<std::vector<size_t>>;
		corePtr = new std::vector<bool>;
	}

	// Resize each vector.
	neighborPtr->clear(); // Just in case there was anything in it.
	neighborPtr->resize(referenceSet->n_cols);
	distancePtr->clear();
	distancePtr->resize(referenceSet->n_cols);
	corePtr->clear();
	corePtr->resize(referenceSet->n_cols);

	for (size_t nbrIndex = 0; nbrIndex < referenceSet->n_cols; nbrIndex++) {
		neighborPtr->at(nbrIndex).reserve(dbscan_minpts);
		distancePtr->at(nbrIndex).reserve(dbscan_minpts);
	}

	// Create the helper object for the traversal.
	typedef EpsMinpts_Rules<MetricType, Tree> RuleType;

	RuleType rules(*referenceSet, *referenceSet, dbscan_eps, dbscan_minpts, *neighborPtr,
				   *distancePtr, *corePtr, metric, true /* don't return the query in the results */);

	if (naive) {
		// The naive brute-force solution.
		for (size_t i = 0; i < referenceSet->n_cols; ++i)
			for (size_t j = 0; j < referenceSet->n_cols; ++j)
				rules.BaseCase(i, j);

		//baseCases = (referenceSet->n_cols * referenceSet->n_cols);
		//scores = 0;
	} else if (singleMode) {
		// Create the traverser.
		typename Tree::template SingleTreeTraverser<RuleType> traverser(rules);

		// Now have it traverse for each point.
		for (size_t i = 0; i < referenceSet->n_cols; ++i)
			traverser.Traverse(i, *referenceTree);

		//baseCases = rules.BaseCases();
		//scores = rules.Scores();
	} else // Dual-tree recursion.
	{
		// Create the traverser.
		typename Tree::template DualTreeTraverser<RuleType> traverser(rules);

		traverser.Traverse(*referenceTree, *referenceTree);

		//baseCases = rules.BaseCases();
		//scores = rules.Scores();
	}

	// Do we need to map the reference indices?
	if (treeOwner && mlpack::TreeTraits<Tree>::RearrangesDataset) {
		neighbors.clear();
		neighbors.resize(referenceSet->n_cols);
		distances.clear();
		distances.resize(referenceSet->n_cols);
		isCore.clear();
		isCore.resize(referenceSet->n_cols);

		for (size_t i = 0; i < distances.size(); ++i) {
			// Map distances (copy a column).
			const size_t refMapping = oldFromNewReferences[i];
			distances[refMapping] = (*distancePtr)[i];
			//isCore[refMapping] = ((*distancePtr)[i].size() >= dbscan_minpts);
			isCore[refMapping] = (*corePtr)[i];
			// Copy each neighbor individually, because we need to map it.
			neighbors[refMapping].resize(distances[refMapping].size());
			for (size_t j = 0; j < distances[refMapping].size(); ++j) {
				neighbors[refMapping][j] = oldFromNewReferences[(*neighborPtr)[i][j]];
			}
		}

		// Finished with temporary objects.
		delete neighborPtr;
		delete distancePtr;
		delete corePtr;
	}
}