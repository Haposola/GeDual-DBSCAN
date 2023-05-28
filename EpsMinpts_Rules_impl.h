#pragma once

template <typename MetricType, typename TreeType>
EpsMinpts_Rules<MetricType, TreeType>::EpsMinpts_Rules(
	const arma::mat &referenceSet,
	const arma::mat &querySet,
	const double dbscan_eps,
	const int dbscan_minpts,
	std::vector<std::vector<size_t>> &neighbors,
	std::vector<std::vector<double>> &distances,
	std::vector<bool> &isCore,
	MetricType &metric,
	const bool sameSet) : referenceSet(referenceSet),
						  querySet(querySet),
						  dbscan_eps(dbscan_eps),
						  dbscan_minpts(dbscan_minpts),
						  neighbors(neighbors),
						  distances(distances),
						  isCore(isCore),
						  metric(metric),
						  sameSet(sameSet),
						  lastQueryIndex(querySet.n_cols),
						  lastReferenceIndex(referenceSet.n_cols){
							  // Nothing to do.
						  };

//! The base case.  Evaluate the distance between the two points and add to the
//! results if necessary.
template <typename MetricType, typename TreeType>
// double EpsMinpts_Rules::BaseCase(
inline mlpack_force_inline double EpsMinpts_Rules<MetricType, TreeType>::BaseCase(
	const size_t queryIndex,
	const size_t referenceIndex)
{
	// If the datasets are the same, don't return the point as in its own range.
	if (sameSet && (queryIndex == referenceIndex))
		return 0.0;

	// If we have just performed this base case, don't do it again.
	if ((lastQueryIndex == queryIndex) && (lastReferenceIndex == referenceIndex))
		return 0.0; // No value to return... this shouldn't do anything bad.

	// If the query point has enough neighbors, stop base_casing
	// if (neighbors[queryIndex].size() >= dbscan_minpts) return 0.0;//Frequently invoking neighbors[i].size() may be slower
	//if (isCore[queryIndex])		return 0.0;

	const double distance = metric.Evaluate(querySet.unsafe_col(queryIndex),
											referenceSet.unsafe_col(referenceIndex));
	//++baseCases;

	// Update last indices, so we don't accidentally perform a base case twice.
	lastQueryIndex = queryIndex;
	lastReferenceIndex = referenceIndex;

	if (distance <= dbscan_eps)
	{
		neighbors[queryIndex].push_back(referenceIndex);
		distances[queryIndex].push_back(distance);
		if (neighbors[queryIndex].size() >= dbscan_minpts)
			isCore[queryIndex] = true;
		// The .size() function is invoked for at most dbscan_minpts times here.
	}

	return distance;
};

//! Single-tree scoring function.
template <typename MetricType, typename TreeType>
// double EpsMinpts_Rules::Score(const size_t queryIndex,
double EpsMinpts_Rules<MetricType, TreeType>::Score(const size_t queryIndex,
													TreeType &referenceNode)
{
	// Prune if this query point has found enough neighbors.
	// if (neighbors[queryIndex].size() >= dbscan_minpts) return DBL_MAX;
	if (isCore[queryIndex])
		return DBL_MAX;

	// We must get the minimum and maximum distances and store them in this object.
	mlpack::Range distances;
	// mlpack::Range range(0.0, dbscan_eps);
	if (mlpack::TreeTraits<TreeType>::FirstPointIsCentroid)
	{
		// In this situation, we calculate the base case.  So we should check to be
		// sure we haven't already done that.
		double baseCase;
		if (mlpack::TreeTraits<TreeType>::HasSelfChildren &&
			(referenceNode.Parent() != NULL) &&
			(referenceNode.Point(0) == referenceNode.Parent()->Point(0)))
		{
			// If the tree has self-children and this is a self-child, the base case
			// was already calculated.
			baseCase = referenceNode.Parent()->Stat().LastDistance();
			lastQueryIndex = queryIndex;
			lastReferenceIndex = referenceNode.Point(0);
		}
		else
		{
			// We must calculate the base case by hand.
			baseCase = BaseCase(queryIndex, referenceNode.Point(0));
		}

		// This may be possibly loose for non-ball bound trees.
		distances.Lo() = baseCase - referenceNode.FurthestDescendantDistance();
		distances.Hi() = baseCase + referenceNode.FurthestDescendantDistance();

		// Update last distance calculation.
		referenceNode.Stat().LastDistance() = baseCase;
	}
	else
	{
		distances = referenceNode.RangeDistance(querySet.unsafe_col(queryIndex));
		//++scores;
	}

	// Prune if the range does not contain the search range(0,eps)
	if (distances.Lo() > dbscan_eps)
		return DBL_MAX;

	// In this case, all of the points in the reference node will be part of the
	// results.
	if (distances.Hi() <= dbscan_eps)
	{
		AddResult(queryIndex, referenceNode);
		return DBL_MAX; // We don't need to go any deeper.
	}

	// Otherwise the score doesn't matter.  Recursion order is irrelevant in
	// range search.
	return 0.0;
}

//! Single-tree rescoring function.
template <typename MetricType, typename TreeType>
double EpsMinpts_Rules<MetricType, TreeType>::Rescore(
	// double EpsMinpts_Rules::Rescore(
	const size_t queryIndex,
	TreeType & /* referenceNode */,
	const double oldScore) const
{
	// Prune if the query point has found enough neighbors.
	// if (neighbors[queryIndex].size() >= dbscan_minpts) return DBL_MAX;
	if (isCore[queryIndex])
		return DBL_MAX;
	// Otherwise it cant't be pruned by distance score.
	return oldScore;
}

//! Dual-tree scoring function.
template <typename MetricType, typename TreeType>
// double EpsMinpts_Rules::Score(TreeType& queryNode,
double EpsMinpts_Rules<MetricType, TreeType>::Score(TreeType &queryNode,
													TreeType &referenceNode)
{
	// Prune if all the points contained in queryNode have found enough neighbors
	if (queryNode.Stat().IsAllDecided())
		return DBL_MAX;

	mlpack::Range distances ;

	if (mlpack::TreeTraits<TreeType>::FirstPointIsCentroid)
	{
		// It is possible that the base case has already been calculated.
		double baseCase = 0.0;
		if ((traversalInfo.LastQueryNode() != NULL) &&
			(traversalInfo.LastReferenceNode() != NULL) &&
			(traversalInfo.LastQueryNode()->Point(0) == queryNode.Point(0)) &&
			(traversalInfo.LastReferenceNode()->Point(0) == referenceNode.Point(0)))
		{
			baseCase = traversalInfo.LastBaseCase();

			// Make sure that if BaseCase() is called, we don't duplicate results.
			lastQueryIndex = queryNode.Point(0);
			lastReferenceIndex = referenceNode.Point(0);
		}
		else
		{
			// We must calculate the base case.
			baseCase = BaseCase(queryNode.Point(0), referenceNode.Point(0));
		}

		distances.Lo() = baseCase - queryNode.FurthestDescendantDistance() - referenceNode.FurthestDescendantDistance();
		distances.Hi() = baseCase + queryNode.FurthestDescendantDistance() + referenceNode.FurthestDescendantDistance();

		// Update the last distances performed for the query and reference node.
		traversalInfo.LastBaseCase() = baseCase;
	}
	else
	{
		// Just perform the calculation.
		distances = referenceNode.RangeDistance(queryNode);
		//++scores;
	}

	// Prune if the range does not contain the search range(0,eps)
	if (distances.Lo() > dbscan_eps)
		return DBL_MAX;

	// In this case, all of the points in the reference node will be part of all
	// the results for each point in the query node.
	if (distances.Hi() <= dbscan_eps)
	{
		for (size_t i = 0; i < queryNode.NumDescendants(); ++i)
			AddResult(queryNode.Descendant(i), referenceNode);
		return DBL_MAX; // We don't need to go any deeper.
	}

	// Otherwise the score doesn't matter.  Recursion order is irrelevant in range search.

	traversalInfo.LastQueryNode() = &queryNode;
	traversalInfo.LastReferenceNode() = &referenceNode;
	return 0.0;
}

//! Dual-tree rescoring function.
template <typename MetricType, typename TreeType>
double EpsMinpts_Rules<MetricType, TreeType>::Rescore(
	// double EpsMinpts_Rules::Rescore(
	TreeType &queryNode,
	TreeType & /* referenceNode */,
	const double oldScore) const
{
	bool all_decided = true;
	const size_t queryEnd = queryNode.Begin() + queryNode.Count();
	for (size_t query = queryNode.Begin(); query < queryEnd; ++query)
	{
		// if (neighbors[query].size() < dbscan_minpts) {
		all_decided &= isCore[query];
		if (!all_decided)
			break;
	}
	queryNode.Stat().IsAllDecided() = all_decided;
	if (all_decided)
		return DBL_MAX;

	// If it wasn't pruned before, it isn't pruned now.
	return oldScore;
}

//! Add all the points in the given node to the results for the given query point.
//! Stop if dbscan_minpts neighbors are found.
template <typename MetricType, typename TreeType>
void EpsMinpts_Rules<MetricType, TreeType>::AddResult(const size_t queryIndex,
													  // void EpsMinpts_Rules::AddResult(const size_t queryIndex,
													  TreeType &referenceNode)
{
	// Stop if the query point has found enough neighbors.
	// if (neighbors[queryIndex].size() >= dbscan_minpts) return;
	if (isCore[queryIndex])
		return;
	// Some types of trees calculate the base case evaluation before Score() is
	// called, so if the base case has already been calculated, then we must avoid
	// adding that point to the results again.

	size_t baseCaseMod = 0;
	if (mlpack::TreeTraits<TreeType>::FirstPointIsCentroid &&
		(queryIndex == lastQueryIndex) &&
		(referenceNode.Point(0) == lastReferenceIndex))
	{
		baseCaseMod = 1;
	}

	for (size_t i = baseCaseMod; i < referenceNode.NumDescendants(); ++i)
	{
		if ((&referenceSet == &querySet) &&
			(queryIndex == referenceNode.Descendant(i)))
			continue;

		const double distance = metric.Evaluate(querySet.unsafe_col(queryIndex),
												referenceNode.Dataset().unsafe_col(referenceNode.Descendant(i)));

		neighbors[queryIndex].push_back(referenceNode.Descendant(i));
		distances[queryIndex].push_back(distance);
		if (neighbors[queryIndex].size() >= dbscan_minpts)
		{
			isCore[queryIndex] = true;
			return;
		}
	}
}
