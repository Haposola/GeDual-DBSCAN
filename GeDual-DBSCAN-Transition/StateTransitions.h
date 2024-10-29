#pragma once
inline void GeDual_DBSCAN::PointTransite_toCore(size_t point) {
	isCore[point] = true;
	for (size_t nbr : epsNeighbors[point]) {
		if (isCore[nbr]) ufset.Union(point, nbr);
	}
}

inline bool  GeDual_DBSCAN::checkAllcore(KDTreeHere* node) {
	return (node->Stat().isAllcore() || node->Stat().getNumCore() >= node->Count());
}
inline bool GeDual_DBSCAN::checkAllcore(KDTreeHere* node, bool force) {
	if (node->Stat().isAllcore()) return true;
	else if (node->Stat().getNumCore() >= node->Count()) return true;
	else {
		size_t num_core = 0; bool allcore = true;
		size_t nodeEnd = node->Begin() + node->Count();
		for (size_t i = node->Begin(); i < nodeEnd; i++) {
			if (!isCore[oldFromNew[i]]) allcore = false;
			else num_core++;
		}
		node->Stat().setNumCore(num_core);
		return allcore;
	}
}

inline void GeDual_DBSCAN::NodeTransite_incNumCore(KDTreeHere* node) {
	node->Stat().incNumCore();
	if (checkAllcore(node)) {
		NodeTransite_toAllcore(node);
	}
}
// Forcing to check onecluster by computing the pairwise distances
inline bool GeDual_DBSCAN::checkOnecluster(KDTreeHere* node) {
	if (node->Stat().isOnecluster()) return true;
	size_t rep = ufset.Find(oldFromNew[node->Begin()]);
	size_t nodeEnd = node->Begin() + node->Count();
	bool onecluster = true;
	for (size_t i = node->Begin() + 1; i < nodeEnd; i++) {
		if (ufset.Find(oldFromNew[i]) != rep) {
			onecluster = false; break;
		}
	}
	if (onecluster) {
		node->Stat().setOnecluster();
		return true;
	} else return false;
}
inline bool GeDual_DBSCAN::checkOnecluster(KDTreeHere* node, bool force) {
	//do base case here
	size_t nodeEnd = node->Begin() + node->Count();
	for (size_t i = node->Begin(); i < nodeEnd; i++) {
		size_t p1 = oldFromNew[i];
		for (size_t j = i + 1; j < nodeEnd; j++) {
			size_t p2 = oldFromNew[j];
			if (ufset.Find(p1) == ufset.Find(p2)) continue;
			double dist = mlpack::EuclideanDistance::Evaluate(dataset.unsafe_col(p1), dataset.unsafe_col(p2));
#if defined TEST_STATISTICS
			num_dist_compute++;
#endif
			if (dist <= dbscan_eps) ufset.Union(p1, p2);
		}
	}
	bool onecluster = true; size_t rep = ufset.Find(oldFromNew[node->Begin()]);
	for (size_t i = node->Begin() + 1; i < nodeEnd; i++) {
		if (ufset.Find(oldFromNew[i]) != rep) {
			onecluster = false; break;
		}
	}
	if (onecluster) {
		node->Stat().setOnecluster();
		return true;
	}
	return node->Stat().isOnecluster();
}


inline void GeDual_DBSCAN::NodeTransite_toAllcore(KDTreeHere* node) {
#if defined TEST_STATISTICS
	num_node_toAllCore++;
#endif
	node->Stat().setAllcore();
	node->Stat().setNumCore(node->Count());
	if (checkOnecluster(node)) {
		NodeTransite_toDense(node);
	}
}
inline void GeDual_DBSCAN::NodeTransite_toDense(KDTreeHere* node) {
#if defined TEST_STATISTICS
	num_node_toDense++;
#endif
	node->Stat().setDense();
	
}
inline void GeDual_DBSCAN::NodeTransite_toDense(KDTreeHere* queryNode, KDTreeHere* node) {
	//cascading transition to dense
	KDTreeHere* parent = node->Parent();
	while (parent != nullptr) {
		if (parent->Begin() > queryNode->Begin()) {
			if (parent->Left()->Stat().isDense() && parent->Right()->Stat().isDense()) {
				size_t leftRep = ufset.Find(oldFromNew[parent->Left()->Begin()]);
				size_t rightRep = ufset.Find(oldFromNew[parent->Right()->Begin()]);
				if (leftRep == rightRep) {
					parent->Stat().setDense();
					parent = parent->Parent();
					continue;
				} else {
					// The way is to check whether one of the two dense siblings is the ancestor of the current queryNode
					// If so, we do not use eps_connect to merge them, since duplicate distance computations may be incurred.
					// To check ancestorship, simply check node->Begin()<=queryNode->Begin()
					// For KDTreeHere ancestorship, node->Begin>=ancestor.Begin, and node->Begin+Count<= ancestor->Begin+Count
					// The other side, node->Begin+Count> not_ancestor->Begin+count need not to be checked, 
					// The case node->Begin+Count> not_ancestor->Begin+count is impossible since queryNode->Begin<refNode.Begin
					mlpack::RangeType<double> lr_range = parent->Left()->RangeDistance(*parent->Right());
#if defined TEST_STATISTICS
					num_Score++; num_dist_compute += 2;
#endif
					if (lr_range.Hi() <= dbscan_eps) {
						ufset.Union(leftRep, rightRep);
						parent->Stat().setDense();
#if defined TEST_STATISTICS
						num_node_toDense++;
#endif
						parent = parent->Parent();
					} else if (lr_range.Lo() <= dbscan_eps) {
						bool found = false;
						BaseCase_EpsConnect_dual(parent->Left(), parent->Right(), found);
						if (found) {
							ufset.Union(leftRep, rightRep);
							parent->Stat().setDense();
#if defined TEST_STATISTICS
							num_node_toDense++;
#endif
							parent = parent->Parent();
						} else break;
					} else break;
				}
			} else break;
		} else break;

		
	}
}

