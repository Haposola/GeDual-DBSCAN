#pragma once
void GeDual_DBSCAN::GeDual_Traversal() {
	// Generalized Dual-tree priority traversal
	// First point to do dual-tree
	// When visiting nodes: 
	// if query==ref, push (query.left, query.left), (query.right,query.left) (query.right,query.right) into priqueue
	//  LR check range distance. If highBound< dbscan_eps, mark as core, merge, and move to restore_stack
	// Otherwise, push LL and RR into stack. Check lowBound, if satisfied, push LR into stack
	// else, visit query.both child and ref. both child, incurring four combinations.
	// Push root first. 
	// While traversing: use node->Begin() as the priority. always recurse the one with the same Begin.
	// In this way, while the top node has a different Begin bound with the last visited combination, we know the last node is decided.
	// This applies to those nodes containing non-core points
	// 
	dual_priqueue.push(std::make_pair(searchTree, searchTree));
	//KDTreeHere* last_visit_query = searchTree;
	while (!dual_priqueue.empty()) {

		auto& top_combination = dual_priqueue.top();
		KDTreeHere* queryNode = top_combination.first;
		KDTreeHere* refNode = top_combination.second;
		dual_priqueue.pop();

		if (queryNode->IsLeaf() || queryNode->Stat().isDense()) { //not descendable
			if (refNode->IsLeaf() || refNode->Stat().isDense()) {//both not descendable
				//Both dense has the Highest priority
				if (queryNode == refNode && (queryNode->Stat().isOnecluster() || queryNode->Stat().isAllcore())) continue;
				if (queryNode->Stat().isDense() && refNode->Stat().isDense()) InterCase_BothDense(queryNode, refNode);
				else if (queryNode->Stat().isDense() && refNode->IsLeaf()) InterCase_QueryDenseRefLeaf(queryNode, refNode);
				else if (queryNode->IsLeaf() && refNode->Stat().isDense()) InterCase_QueryLeafRefDense(queryNode, refNode);
				else InterCase_BothLeaves(queryNode, refNode);
			} else {//descend refNode
				mlpack::RangeType<double> qr_range = queryNode->RangeDistance(*refNode->Right());
				if (qr_range.Hi() <= dbscan_eps) BaseCase_NoDistCompute(queryNode, refNode->Right());
				else if (qr_range.Lo() <= dbscan_eps) dual_priqueue.push(std::make_pair(queryNode, refNode->Right()));

				mlpack::RangeType<double> ql_range = queryNode->RangeDistance(*refNode->Left());
				if (ql_range.Hi() <= dbscan_eps) BaseCase_NoDistCompute(queryNode, refNode->Left());
				else if (ql_range.Lo() <= dbscan_eps) dual_priqueue.push(std::make_pair(queryNode, refNode->Left()));
			}
		} else {//query is descendable
			if (queryNode == refNode) {
				if (queryNode->Stat().isOnecluster() || queryNode->Stat().isAllcore()) continue;
				//Push LL and RR to ensure the traversal is correct.
				dual_priqueue.push(std::make_pair(queryNode->Left(), queryNode->Left()));
				dual_priqueue.push(std::make_pair(queryNode->Right(), queryNode->Right()));
				//visit LR
				mlpack::RangeType<double> lr_range = queryNode->Left()->RangeDistance(*queryNode->Right());
				if (lr_range.Hi() <= dbscan_eps) BaseCase_NoDistCompute(queryNode->Left(), queryNode->Right());
				else if (lr_range.Lo() <= dbscan_eps) dual_priqueue.push(std::make_pair(queryNode->Left(), queryNode->Right()));

			} else {//descend queryNode
				mlpack::RangeType<double> lref_range = queryNode->Left()->RangeDistance(*refNode);
				if (lref_range.Hi() <= dbscan_eps) BaseCase_NoDistCompute(queryNode->Left(), refNode);
				else if (lref_range.Lo() <= dbscan_eps) dual_priqueue.push(std::make_pair(queryNode->Left(), refNode));

				mlpack::RangeType<double> rref_range = queryNode->Right()->RangeDistance(*refNode);
				if (rref_range.Hi() <= dbscan_eps) BaseCase_NoDistCompute(queryNode->Right(), refNode);
				else if (rref_range.Lo() <= dbscan_eps) dual_priqueue.push(std::make_pair(queryNode->Right(), refNode));
			}


		}
	}
}