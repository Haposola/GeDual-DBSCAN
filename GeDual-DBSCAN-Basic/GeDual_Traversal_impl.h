#pragma once
#if defined TEST_STATISTICS
void DBSCAN_GeDual::GeDual_Traversal() {

	dual_priqueue.push(std::make_pair(searchTree, searchTree));
	while (!dual_priqueue.empty()) {

		auto& top_combination = dual_priqueue.top();
		KDTreeHere* queryNode = top_combination.first;
		KDTreeHere* refNode = top_combination.second;
		dual_priqueue.pop();

		if (queryNode->IsLeaf() && refNode->IsLeaf()) {
			BaseCase_EpsMinpts(queryNode, refNode, true);
		} else if (!queryNode->IsLeaf()) {
			if (queryNode == refNode) {
				dual_priqueue.push(std::make_pair(queryNode->Left(), queryNode->Left()));
				dual_priqueue.push(std::make_pair(queryNode->Right(), queryNode->Right()));
				//visit LR
				mlpack::RangeType<double> lr_range = queryNode->Left()->RangeDistance(*queryNode->Right());
				num_Score++; num_dist_compute += 2;
				if (lr_range.Hi() <= dbscan_eps) {
					num_near_case++;
					BaseCase_EpsMinpts(queryNode->Left(), queryNode->Right(), false);
				} else if (lr_range.Lo() <= dbscan_eps) {
					dual_priqueue.push(std::make_pair(queryNode->Left(), queryNode->Right()));
				} else num_far_case++;

			} else {//descend queryNode
				mlpack::RangeType<double> lref_range = queryNode->Left()->RangeDistance(*refNode);
				num_Score++; num_dist_compute += 2;
				if (lref_range.Hi() <= dbscan_eps) {
					num_near_case++;
					BaseCase_EpsMinpts(queryNode->Left(), refNode, false);
				} else if (lref_range.Lo() <= dbscan_eps) {
					dual_priqueue.push(std::make_pair(queryNode->Left(), refNode));
				} else num_far_case++;

				mlpack::RangeType<double> rref_range = queryNode->Right()->RangeDistance(*refNode);
				num_Score++; num_dist_compute += 2;

				if (rref_range.Hi() <= dbscan_eps) {
					num_near_case++;
					BaseCase_EpsMinpts(queryNode->Right(), refNode, false);
				} else if (rref_range.Lo() <= dbscan_eps) {
					dual_priqueue.push(std::make_pair(queryNode->Right(), refNode));
				} else num_far_case++;

			}
		} else {//descend refNode
			mlpack::RangeType<double> qr_range = queryNode->RangeDistance(*refNode->Right());
			num_Score++;num_dist_compute += 2;
			if (qr_range.Hi() <= dbscan_eps) {
				num_near_case++;
				BaseCase_EpsMinpts(queryNode, refNode->Right(), false);
			} else if (qr_range.Lo() <= dbscan_eps) {
				dual_priqueue.push(std::make_pair(queryNode, refNode->Right()));
			} else num_far_case++;

			mlpack::RangeType<double> ql_range = queryNode->RangeDistance(*refNode->Left());
			num_Score++;num_dist_compute += 2;
			if (ql_range.Hi() <= dbscan_eps) {
				num_near_case++;
				BaseCase_EpsMinpts(queryNode, refNode->Left(), false);
			} else if (ql_range.Lo() <= dbscan_eps) {
				dual_priqueue.push(std::make_pair(queryNode, refNode->Left()));
			} else num_far_case++;
		}
	}
}
#else 
void DBSCAN_GeDual::GeDual_Traversal() {

	dual_priqueue.push(std::make_pair(searchTree, searchTree));
	while (!dual_priqueue.empty()) {

		auto& top_combination = dual_priqueue.top();
		KDTreeHere* queryNode = top_combination.first;
		KDTreeHere* refNode = top_combination.second;
		dual_priqueue.pop();

		if (queryNode->IsLeaf() && refNode->IsLeaf()) {
			BaseCase_EpsMinpts(queryNode, refNode, true);
		} else if (!queryNode->IsLeaf()) {
			if (queryNode == refNode) {
				dual_priqueue.push(std::make_pair(queryNode->Left(), queryNode->Left()));
				dual_priqueue.push(std::make_pair(queryNode->Right(), queryNode->Right()));
				//visit LR
				mlpack::RangeType<double> lr_range = queryNode->Left()->RangeDistance(*queryNode->Right());
				if (lr_range.Hi() <= dbscan_eps) BaseCase_EpsMinpts(queryNode->Left(), queryNode->Right(), false);
				else if (lr_range.Lo() <= dbscan_eps) dual_priqueue.push(std::make_pair(queryNode->Left(), queryNode->Right()));
			} else {//descend queryNode
				mlpack::RangeType<double> lref_range = queryNode->Left()->RangeDistance(*refNode);
				if (lref_range.Hi() <= dbscan_eps) BaseCase_EpsMinpts(queryNode->Left(), refNode, false);
				else if (lref_range.Lo() <= dbscan_eps) dual_priqueue.push(std::make_pair(queryNode->Left(), refNode));
				
				mlpack::RangeType<double> rref_range = queryNode->Right()->RangeDistance(*refNode);
				if (rref_range.Hi() <= dbscan_eps) BaseCase_EpsMinpts(queryNode->Right(), refNode, false);
				else if (rref_range.Lo() <= dbscan_eps) dual_priqueue.push(std::make_pair(queryNode->Right(), refNode));
			}
		} else {//descend refNode
			mlpack::RangeType<double> qr_range = queryNode->RangeDistance(*refNode->Right());
			if (qr_range.Hi() <= dbscan_eps) BaseCase_EpsMinpts(queryNode, refNode->Right(), false);
			else if (qr_range.Lo() <= dbscan_eps) dual_priqueue.push(std::make_pair(queryNode, refNode->Right()));
			
			mlpack::RangeType<double> ql_range = queryNode->RangeDistance(*refNode->Left());
			if (ql_range.Hi() <= dbscan_eps) BaseCase_EpsMinpts(queryNode, refNode->Left(), false);
			else if (ql_range.Lo() <= dbscan_eps) dual_priqueue.push(std::make_pair(queryNode, refNode->Left()));
		}
	}
}
#endif 
