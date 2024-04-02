#pragma once
inline void GeDual_DBSCAN::BaseCase_PointwiseEpsminpts(size_t queryPoint, size_t refPoint) {
	if (isCore[queryPoint]) {
		//skip intra-cluster distance computation
		if (ufset.Find(queryPoint) == ufset.Find(refPoint)) return;
		if (mlpack::EuclideanDistance::Evaluate(dataset.unsafe_col(queryPoint), dataset.unsafe_col(refPoint)) <= dbscan_eps) {
			if (isCore[refPoint]) ufset.Union(queryPoint, refPoint);
			else {
				//Invariant: After a point is changed into core-point, the neighbors will not matter
				//If a non-core neighbor point is visited, the neighbor of non-core will be updated, 
				//and the merge operation will be executed on the time the non-core point is change into core-point.
				//If a core-point neighbor is visited, the merge operation will be conducted immediately. 
				epsNeighbors[refPoint].push_back(queryPoint);
				//if refPoints is changed into a core-point, conduct not-yet-done merge.
				if (++numNeighbors[refPoint] >= dbscan_minpts) {
					PointTransite_ToCore(refPoint);
				}
			}
		}
	} else {
		if (mlpack::EuclideanDistance::Evaluate(dataset.unsafe_col(queryPoint), dataset.unsafe_col(refPoint)) <= dbscan_eps) {
			epsNeighbors[queryPoint].push_back(refPoint);
			if (++numNeighbors[queryPoint] >= dbscan_minpts) {
				PointTransite_ToCore(queryPoint);
			}
			if (!isCore[refPoint]) {
				epsNeighbors[refPoint].push_back(queryPoint);
				if (++numNeighbors[refPoint] >= dbscan_minpts) {
					PointTransite_ToCore(refPoint);
				}
			}
		}
	}
}


void GeDual_DBSCAN::BaseCase_NoDistCompute(KDTreeHere* queryNode, KDTreeHere* refNode) {
	size_t queryEnd = queryNode->Begin() + queryNode->Count();
	size_t refEnd = refNode->Begin() + refNode->Count();
	//first change non-core to cores if possible
	if (!queryNode->Stat().isDense()) {
		for (size_t i = queryNode->Begin(); i < queryEnd; i++) {
			size_t queryPoint = oldFromNew[i];
			if (!isCore[queryPoint]) {
				numNeighbors[queryPoint] += refNode->Count();
				if (numNeighbors[queryPoint] >= dbscan_minpts) {
					//Only merge existing neighbors;
					//The newly added neighbors in refNode will be processed outside this loop;
					PointTransite_ToCore(queryPoint);
				} else {//Not changed to core-point; record the neighobrs
					for (size_t j = refNode->Begin(); j < refEnd; j++) 
						epsNeighbors[queryPoint].push_back(oldFromNew[j]);
				}
			}
		}
	}
	//process refNode similarly;
	if (!refNode->Stat().isDense()) {
		for (size_t i = refNode->Begin(); i < refEnd; i++) {
			size_t refPoint = oldFromNew[i];
			if (!isCore[refPoint]) {
				numNeighbors[refPoint] += queryNode->Count();
				if (numNeighbors[refPoint] >= dbscan_minpts) {
					PointTransite_ToCore(refPoint);
				} else {
					for (size_t j = queryNode->Begin(); j < queryEnd; j++) 
						epsNeighbors[refPoint].push_back(oldFromNew[j]);
				}
			}
		}
	}
	//Merge all queryCore and refCore pairs, unless query/ref is knows to be of onecluster
	for (size_t i = queryNode->Begin(); i < queryEnd; i++) {
		size_t queryPoint = oldFromNew[i];
		if (!isCore[queryPoint])continue;
		for (size_t j = refNode->Begin(); j < refEnd; j++) {
			size_t refPoint = oldFromNew[j];
			if (!isCore[refPoint]) continue;
			ufset.Union(queryPoint, refPoint);
			if (refNode->Stat().isOnecluster())break;//the break sentence will be executed for the first core-point visited in a onecluster node
		}
		if (queryNode->Stat().isOnecluster())break;
	}
}

//invoked in QueryDenseRefLeaf or QueryLeafRefDense, when para_point is non-core
// Functionality: 
// find enough eps-neighbors for the point be become a core-point
// then interrupt the search immediately and merge core-point with dense-node
void GeDual_DBSCAN::BaseCase_EpsMinptsInDense_single(size_t point, KDTreeHere* node, bool& become_core) {
	if (node->IsLeaf()) {
		size_t refEnd = node->Begin() + node->Count();
		for (size_t i = node->Begin(); i < refEnd; i++) {
			size_t refPoint = oldFromNew[i];
			double dist = mlpack::EuclideanDistance::Evaluate(dataset.unsafe_col(point), dataset.unsafe_col(refPoint));
			if (dist <= dbscan_eps) {
				if (++numNeighbors[point] >= dbscan_minpts) {
					become_core = true;
					PointTransite_ToCore(point);
					return;
				} else 	epsNeighbors[point].push_back(refPoint);
			}
		}
		return;
	} else {
		bool leftValid = true;
		mlpack::RangeType<double> l_range = node->Left()->RangeDistance(dataset.unsafe_col(point));
		if (l_range.Hi() <= dbscan_eps) {
			leftValid = false;
			numNeighbors[point] += node->Left()->Count();
			if (numNeighbors[point] >= dbscan_minpts) {
				PointTransite_ToCore(point);
				become_core = true; return;
			} else {//Not changed to core-point; record the neighobrs
				size_t leftEnd = node->Left()->Begin() + node->Left()->Count();
				for (size_t j = node->Left()->Begin(); j < leftEnd; j++) {
					epsNeighbors[point].push_back(oldFromNew[j]);
				}
			}
		}
		if (l_range.Lo() > dbscan_eps) leftValid = false;

		bool rightValid = true;
		mlpack::RangeType<double> r_range = node->Right()->RangeDistance(dataset.unsafe_col(point));
		if (r_range.Hi() <= dbscan_eps) {
			rightValid = false;
			numNeighbors[point] += node->Right()->Count();
			if (numNeighbors[point] >= dbscan_minpts) {
				PointTransite_ToCore(point);
				become_core = true; return;
			} else {//Not changed to core-point; record the neighobrs
				size_t rightEnd = node->Right()->Begin() + node->Right()->Count();
				for (size_t j = node->Right()->Begin(); j < rightEnd; j++) {
					epsNeighbors[point].push_back(oldFromNew[j]);
				}
			}
		}
		if (r_range.Lo() > dbscan_eps) rightValid = false;

		if (leftValid && rightValid) {
			if (r_range.Lo() < l_range.Lo()) {
				//visit right first
				BaseCase_EpsMinptsInDense_single(point, node->Right(), become_core);
				if (become_core) return;
				BaseCase_EpsMinptsInDense_single(point, node->Left(), become_core);
				if (become_core) return;
			} else {//visit left first
				BaseCase_EpsMinptsInDense_single(point, node->Left(), become_core);
				if (become_core) return;
				BaseCase_EpsMinptsInDense_single(point, node->Right(), become_core);
				if (become_core) return;
			}
		} else if (leftValid && !rightValid) BaseCase_EpsMinptsInDense_single(point, node->Left(), become_core);
		else if (!leftValid && rightValid) BaseCase_EpsMinptsInDense_single(point, node->Right(), become_core);
		else return;
	}
}
//Invoked in QueryDenseRefLeaf or QueryLeafRefDense, when point is core-point
void GeDual_DBSCAN::BaseCase_EpsConnect_single(size_t point, KDTreeHere* node, bool& found) {
	if (node->IsLeaf()) {
		size_t refEnd = node->Begin() + node->Count();
		for (size_t i = node->Begin(); i < refEnd; i++) {
			double dist = mlpack::EuclideanDistance::Evaluate(dataset.unsafe_col(point), dataset.unsafe_col(oldFromNew[i]));
			if (dist <= dbscan_eps) {
				found = true; return;
			}
		}
	} else {
		mlpack::RangeType<double> l_range = node->Left()->RangeDistance(dataset.unsafe_col(point));
		if (l_range.Hi() <= dbscan_eps) {
			found = true; return;
		}
		mlpack::RangeType<double> r_range = node->Right()->RangeDistance(dataset.unsafe_col(point));
		if (r_range.Hi() <= dbscan_eps) {
			found = true; return;
		}
		if (r_range.Lo() < l_range.Lo()) {
			if (r_range.Lo() <= dbscan_eps) {
				BaseCase_EpsConnect_single(point, node->Right(), found);
				if (found) return;
				if (l_range.Lo() <= dbscan_eps) {
					BaseCase_EpsConnect_single(point, node->Left(), found);
					if (found) return;
				}
			}
		} else {
			if (l_range.Lo() <= dbscan_eps) {
				BaseCase_EpsConnect_single(point, node->Left(), found);
				if (found) return;
				if (r_range.Lo() <= dbscan_eps) {
					BaseCase_EpsConnect_single(point, node->Right(), found);
					if (found) return;
				}
			}
		}
	}
}


void GeDual_DBSCAN::BaseCase_EpsConnect_dual(KDTreeHere* queryNode, KDTreeHere* refNode, bool& found) {
	if (queryNode->IsLeaf() && refNode->IsLeaf()) {
		size_t queryEnd = queryNode->Begin() + queryNode->Count();
		size_t refEnd = refNode->Begin() + refNode->Count();
		for (size_t i = queryNode->Begin(); i < queryEnd; i++) {
			for (size_t j = refNode->Begin(); j < refEnd; j++) {
				double dist = mlpack::EuclideanDistance::Evaluate(dataset.unsafe_col(oldFromNew[i]), dataset.unsafe_col(oldFromNew[j]));
				if (dist <= dbscan_eps) {
					found = true; return;
				}
			}
		}
	} else if (!queryNode->IsLeaf() && !refNode->IsLeaf()) {
		std::vector<std::tuple<KDTreeHere*, KDTreeHere*, double>> score_list;
		//visit LR
		mlpack::RangeType<double> lr_range = queryNode->Left()->RangeDistance(*refNode->Right());
		if (lr_range.Hi() <= dbscan_eps) {
			found = true; return;
		} else if (lr_range.Lo() <= dbscan_eps) score_list.push_back(std::make_tuple(queryNode->Left(), refNode->Right(), lr_range.Lo()));
		//visit LL
		mlpack::RangeType<double> ll_range = queryNode->Left()->RangeDistance(*refNode->Left());
		if (ll_range.Hi() <= dbscan_eps) {
			found = true; return;
		} else if (ll_range.Lo() <= dbscan_eps) score_list.push_back(std::make_tuple(queryNode->Left(), refNode->Left(), ll_range.Lo()));
		//visit RL
		mlpack::RangeType<double> rl_range = queryNode->Right()->RangeDistance(*refNode->Left());
		if (rl_range.Hi() <= dbscan_eps) {
			found = true; return;
		} else if (rl_range.Lo() <= dbscan_eps) score_list.push_back(std::make_tuple(queryNode->Right(), refNode->Left(), rl_range.Lo()));
		//visit RR
		mlpack::RangeType<double> rr_range = queryNode->Right()->RangeDistance(*refNode->Right());
		if (rr_range.Hi() <= dbscan_eps) {
			found = true; return;
		} else if (rr_range.Lo() <= dbscan_eps) score_list.push_back(std::make_tuple(queryNode->Right(), refNode->Right(), rr_range.Lo()));

		std::sort(score_list.begin(), score_list.end(),
				  [](const std::tuple<KDTreeHere*, KDTreeHere*, double>& p1, const std::tuple<KDTreeHere*, KDTreeHere*, double>& p2) {
					  return std::get<2>(p1) < std::get<2>(p2);
				  });

		for (size_t i = 0; i < score_list.size(); i++) {
			BaseCase_EpsConnect_dual(std::get<0>(score_list[i]), std::get<1>(score_list[i]), found);
			if (found) return;
		}
	} else {//one leaf and one not
		auto leafOne = queryNode->IsLeaf() ? queryNode : refNode;
		auto nonLeafOne = queryNode->IsLeaf() ? refNode : queryNode;

		mlpack::RangeType<double> r_range = nonLeafOne->Right()->RangeDistance(*leafOne);
		if (r_range.Hi() <= dbscan_eps) {
			found = true; return;
		}
		mlpack::RangeType<double> l_range = nonLeafOne->Left()->RangeDistance(*leafOne);
		if (l_range.Hi() <= dbscan_eps) {
			found = true; return;
		}
		if (r_range.Lo() < l_range.Lo()) {
			if (r_range.Lo() <= dbscan_eps) {
				BaseCase_EpsConnect_dual(leafOne, nonLeafOne->Right(), found);
				if (found) return;
				if (l_range.Lo() <= dbscan_eps) {
					BaseCase_EpsConnect_dual(leafOne, nonLeafOne->Left(), found);
					if (found) return;
				}
			}
		} else {
			if (l_range.Lo() <= dbscan_eps) {
				BaseCase_EpsConnect_dual(leafOne, nonLeafOne->Left(), found);
				if (found) return;
				if (r_range.Lo() <= dbscan_eps) {
					BaseCase_EpsConnect_dual(leafOne, nonLeafOne->Right(), found);
					if (found) return;
				}
			}
		}
	}
}