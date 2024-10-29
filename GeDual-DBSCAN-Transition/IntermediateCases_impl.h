#pragma once
void GeDual_DBSCAN::InterCase_BothLeaves(KDTreeHere* queryNode, KDTreeHere* refNode) {
#if defined TEST_STATISTICS
	num_BaseCase++;
#endif
	size_t queryEnd = queryNode->Begin() + queryNode->Count();
	size_t refEnd = refNode->Begin() + refNode->Count();

	for (size_t i = queryNode->Begin(); i < queryEnd; i++) {
		for (size_t j = std::max(i+1, refNode->Begin()); j < refEnd; j++) {
			BaseCase_PointwiseEpsminpts(oldFromNew[i], oldFromNew[j]);
		}//TODO:: add ref_turned codes
	}
	if (checkAllcore(queryNode, true)) NodeTransite_toAllcore(queryNode);
	if (checkAllcore(refNode, true)) {
		NodeTransite_toAllcore(refNode);
		// Only do forced check_onecluster for the refNode
		if (checkOnecluster(refNode, true)) 
			NodeTransite_toDense(queryNode, refNode);
	}	
}
void GeDual_DBSCAN::InterCase_QueryLeafRefDense(KDTreeHere* queryNode, KDTreeHere* refNode) {
#if defined TEST_STATISTICS
	num_intercase_leafdense++;
#endif
	size_t queryEnd = queryNode->Begin() + queryNode->Count();
	for (size_t i = queryNode->Begin(); i < queryEnd; i++) {
		size_t queryPoint = oldFromNew[i], refRep = oldFromNew[refNode->Begin()];
		//Find at least one eps-neighbor in refNode to merge the query-core-point with ref-dense-node;
		if (isCore[queryPoint]) {
			if (ufset.Find(queryPoint) == ufset.Find(refRep)) continue;
			bool found = false;
			BaseCase_EpsConnect_single(queryPoint, refNode, found);
			if (found) ufset.Union(queryPoint, refRep);
		}
		// Find enough eps-neighbors to change queryPoint to core-point
		// then interrupt and merge query-core-point with ref-dense-node;
		else {
			bool cored = false;
			BaseCase_EpsMinptsInDense_single(queryPoint, refNode, cored);
			if (cored) {
				ufset.Union(queryPoint, refRep);
				NodeTransite_incNumCore(queryNode);
			}
		}
	}
}
void GeDual_DBSCAN::InterCase_QueryDenseRefLeaf(KDTreeHere* queryNode, KDTreeHere* refNode) {
#if defined TEST_STATISTICS
	num_intercase_denseleaf++;
#endif
	size_t refEnd = refNode->Begin() + refNode->Count();
	for (size_t j = refNode->Begin(); j < refEnd; j++) {
		size_t refPoint = oldFromNew[j], queryRep = oldFromNew[queryNode->Begin()];
		if (isCore[refPoint]) {
			if (ufset.Find(refPoint) == ufset.Find(queryRep)) continue;
			bool found = false;
			BaseCase_EpsConnect_single(refPoint, queryNode, found);
			if (found) ufset.Union(refPoint, queryRep);
		} else {
			bool cored = false;
			BaseCase_EpsMinptsInDense_single(refPoint, queryNode, cored);
			if (cored) {
#if defined TEST_STATISTICS
				num_early_core++;
#endif
				ufset.Union(refPoint, queryRep);
				NodeTransite_incNumCore(refNode);
			}
		}
	}
	if (refNode->Stat().isAllcore() && checkOnecluster(refNode, true)) {
		NodeTransite_toDense(queryNode, refNode);
	}
}
void GeDual_DBSCAN::InterCase_BothDense(KDTreeHere* queryNode, KDTreeHere* refNode) {
#if defined TEST_STATISTICS
	num_intercase_bothdense++;
#endif
	if (queryNode->Parent() == refNode->Parent())return;
	size_t queryRep = oldFromNew[queryNode->Begin()], refRep = oldFromNew[refNode->Begin()];
	if (ufset.Find(queryRep) != ufset.Find(refRep)) {
		mlpack::RangeType<double> lr_range = queryNode->RangeDistance(*refNode);
#if defined TEST_STATISTICS
		num_Score++; num_dist_compute += 2;
#endif
		if (lr_range.Hi() <= dbscan_eps) {
			ufset.Union(queryRep, refRep);
		} else if (lr_range.Lo() <= dbscan_eps) {
			bool found = false;
			BaseCase_EpsConnect_dual(queryNode, refNode, found);
			if (found) ufset.Union(queryRep, refRep);
		}
#if defined TEST_STATISTICS
		else num_far_case++;
#endif
	}
}
