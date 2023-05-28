#pragma once

template <typename MetricType, typename StatisticType, typename MatType>
void eps_connectness_BSTree_dual_DFS(mlpack::BinarySpaceTree<MetricType, StatisticType, MatType> *queryNode,
									 mlpack::BinarySpaceTree<MetricType, StatisticType, MatType> *refNode,
									 const double &denThres, bool &found)
{
	/*double minDis = queryNode->MinDistance(*refNode);
	double maxDis = queryNode->MaxDistance(*refNode);
	double qFur = queryNode->FurthestDescendantDistance();
	double rFur = refNode->FurthestDescendantDistance();
	if (minDis > denThres + qFur + rFur) return false;
	if (maxDis + qFur + rFur <= denThres) return true;*/
	if (found)
		return;
	if (queryNode->IsLeaf() && refNode->IsLeaf())
	{

		double scoreLow = queryNode->MinDistance(*refNode) - queryNode->FurthestDescendantDistance() - refNode->FurthestDescendantDistance();
		if (scoreLow > denThres)
			return;
		double scoreHigh = queryNode->MaxDistance(*refNode) + queryNode->FurthestDescendantDistance() + refNode->FurthestDescendantDistance();
		if (scoreHigh <= denThres)
		{
			found = true;
			return;
		}
		int num1 = queryNode->NumPoints();
		int num2 = refNode->NumPoints();
		for (int i = 0; i < num1; i++)
		{
			for (int j = 0; j < num2; j++)
				if (MetricType::Evaluate(queryNode->Dataset().col(queryNode->Point(i)), refNode->Dataset().col(refNode->Point(j))) <= denThres)
				{
					found = true;
					return;
				}
		}
		return;
	}
	else
	{
		// otherwise traverse children
		if (queryNode->IsLeaf() && !refNode->IsLeaf())
		{
			double scoreLLow = queryNode->MinDistance(*refNode->Left()) - queryNode->FurthestDescendantDistance() - refNode->Left()->FurthestDescendantDistance();
			double scoreRLow = queryNode->MinDistance(*refNode->Right()) - queryNode->FurthestDescendantDistance() - refNode->Right()->FurthestDescendantDistance();
			if (scoreLLow > denThres && scoreRLow > denThres)
				return;

			double scoreLHigh = queryNode->MaxDistance(*refNode->Left()) + queryNode->FurthestDescendantDistance() + refNode->Left()->FurthestDescendantDistance();
			if (scoreLHigh <= denThres)
			{
				found = true;
				return;
			}
			double scoreRHigh = queryNode->MaxDistance(*refNode->Right()) + queryNode->FurthestDescendantDistance() + refNode->Right()->FurthestDescendantDistance();
			if (scoreRHigh <= denThres)
			{
				found = true;
				return;
			}

			if (scoreLLow <= scoreRLow)
			{
				eps_connectness_BSTree_dual_DFS(queryNode, refNode->Left(), denThres, found);
				if (found)
					return;
				if (scoreRLow <= denThres)
				{
					eps_connectness_BSTree_dual_DFS(queryNode, refNode->Right(), denThres, found);
					if (found)
						return;
				}
			}
			else
			{
				eps_connectness_BSTree_dual_DFS(queryNode, refNode->Right(), denThres, found);
				if (found)
					return;
				if (scoreLLow <= denThres)
				{
					eps_connectness_BSTree_dual_DFS(queryNode, refNode->Left(), denThres, found);
					if (found)
						return;
				}
			}
			return;
		}
		if (!queryNode->IsLeaf() && refNode->IsLeaf())
		{
			double scoreLLow = refNode->MinDistance(*queryNode->Left()) - queryNode->Left()->FurthestDescendantDistance() - refNode->FurthestDescendantDistance();
			double scoreRLow = refNode->MinDistance(*queryNode->Right()) - queryNode->Right()->FurthestDescendantDistance() - refNode->FurthestDescendantDistance();
			if (scoreLLow > denThres && scoreRLow > denThres)
				return;

			double scoreLHigh = refNode->MaxDistance(*queryNode->Left()) + queryNode->Left()->FurthestDescendantDistance() + refNode->FurthestDescendantDistance();
			if (scoreLHigh <= denThres)
			{
				found = true;
				return;
			}
			double scoreRHigh = refNode->MaxDistance(*queryNode->Right()) + queryNode->Right()->FurthestDescendantDistance() + refNode->FurthestDescendantDistance();
			if (scoreRHigh <= denThres)
			{
				found = true;
				return;
			}

			if (scoreLLow <= scoreRLow)
			{
				eps_connectness_BSTree_dual_DFS(queryNode->Left(), refNode, denThres, found);
				if (found)
					return;
				if (scoreRLow <= denThres)
				{
					eps_connectness_BSTree_dual_DFS(queryNode->Right(), refNode, denThres, found);
					if (found)
						return;
				}
			}
			else
			{
				eps_connectness_BSTree_dual_DFS(queryNode->Right(), refNode, denThres, found);
				if (found)
					return;
				if (scoreLLow <= denThres)
				{
					eps_connectness_BSTree_dual_DFS(queryNode->Left(), refNode, denThres, found);
					if (found)
						return;
				}
			}
			return;
		}
		if (!queryNode->IsLeaf() && !refNode->IsLeaf())
		{
			double scoreLLLow = queryNode->Left()->MinDistance(*refNode->Left()) - queryNode->Left()->FurthestDescendantDistance() - refNode->Left()->FurthestDescendantDistance();
			double scoreLRLow = queryNode->Left()->MinDistance(*refNode->Right()) - queryNode->Left()->FurthestDescendantDistance() - refNode->Right()->FurthestDescendantDistance();
			double scoreRLLow = queryNode->Right()->MinDistance(*refNode->Left()) - queryNode->Right()->FurthestDescendantDistance() - refNode->Left()->FurthestDescendantDistance();
			double scoreRRLow = queryNode->Right()->MinDistance(*refNode->Right()) - queryNode->Right()->FurthestDescendantDistance() - refNode->Right()->FurthestDescendantDistance();
			if (scoreLLLow > denThres && scoreLRLow > denThres && scoreRLLow > denThres && scoreRRLow > denThres)
				return;

			double scoreLLHigh = queryNode->Left()->MaxDistance(*refNode->Left()) + queryNode->Left()->FurthestDescendantDistance() + refNode->Left()->FurthestDescendantDistance();
			if (scoreLLHigh <= denThres)
			{
				found = true;
				return;
			}
			double scoreLRHigh = queryNode->Left()->MaxDistance(*refNode->Right()) + queryNode->Left()->FurthestDescendantDistance() + refNode->Right()->FurthestDescendantDistance();
			if (scoreLRHigh <= denThres)
			{
				found = true;
				return;
			}
			double scoreRLHigh = queryNode->Right()->MaxDistance(*refNode->Left()) + queryNode->Right()->FurthestDescendantDistance() + refNode->Left()->FurthestDescendantDistance();
			if (scoreRLHigh <= denThres)
			{
				found = true;
				return;
			}
			double scoreRRHigh = queryNode->Right()->MaxDistance(*refNode->Right()) + queryNode->Right()->FurthestDescendantDistance() + refNode->Right()->FurthestDescendantDistance();
			if (scoreRRHigh <= denThres)
			{
				found = true;
				return;
			}

			std::vector<std::pair<int, double>> score_list(4);
			score_list[0] = std::make_pair(0, scoreLLLow);
			score_list[1] = std::make_pair(1, scoreLRLow);
			score_list[2] = std::make_pair(2, scoreRLLow);
			score_list[3] = std::make_pair(3, scoreRRLow);
			std::sort(score_list.begin(), score_list.end(), [](const std::pair<int, double> &p1, const std::pair<int, double> &p2)
					  { return p2.second == p1.second ? p1.first < p2.first : p1.second < p2.second; });

			switch (score_list[0].first)
			{
			case 0:
				eps_connectness_BSTree_dual_DFS(queryNode->Left(), refNode->Left(), denThres, found);
				if (found)
					return;
				break;
			case 1:
				eps_connectness_BSTree_dual_DFS(queryNode->Left(), refNode->Right(), denThres, found);
				if (found)
					return;
				break;
			case 2:
				eps_connectness_BSTree_dual_DFS(queryNode->Right(), refNode->Left(), denThres, found);
				if (found)
					return;
				break;
			case 3:
				eps_connectness_BSTree_dual_DFS(queryNode->Right(), refNode->Right(), denThres, found);
				if (found)
					return;
				break;
			}
			if (score_list[1].second <= denThres)
			{
				switch (score_list[1].first)
				{
				case 0:
					eps_connectness_BSTree_dual_DFS(queryNode->Left(), refNode->Left(), denThres, found);
					if (found)
						return;
					break;
				case 1:
					eps_connectness_BSTree_dual_DFS(queryNode->Left(), refNode->Right(), denThres, found);
					if (found)
						return;
					break;
				case 2:
					eps_connectness_BSTree_dual_DFS(queryNode->Right(), refNode->Left(), denThres, found);
					if (found)
						return;
					break;
				case 3:
					eps_connectness_BSTree_dual_DFS(queryNode->Right(), refNode->Right(), denThres, found);
					if (found)
						return;
					break;
				}
				if (score_list[2].second <= denThres)
				{
					switch (score_list[2].first)
					{
					case 0:
						eps_connectness_BSTree_dual_DFS(queryNode->Left(), refNode->Left(), denThres, found);
						if (found)
							return;
						break;
					case 1:
						eps_connectness_BSTree_dual_DFS(queryNode->Left(), refNode->Right(), denThres, found);
						if (found)
							return;
						break;
					case 2:
						eps_connectness_BSTree_dual_DFS(queryNode->Right(), refNode->Left(), denThres, found);
						if (found)
							return;
						break;
					case 3:
						eps_connectness_BSTree_dual_DFS(queryNode->Right(), refNode->Right(), denThres, found);
						if (found)
							return;
						break;
					}
					if (score_list[3].second <= denThres)
					{
						switch (score_list[3].first)
						{
						case 0:
							eps_connectness_BSTree_dual_DFS(queryNode->Left(), refNode->Left(), denThres, found);
							if (found)
								return;
							break;
						case 1:
							eps_connectness_BSTree_dual_DFS(queryNode->Left(), refNode->Right(), denThres, found);
							if (found)
								return;
							break;
						case 2:
							eps_connectness_BSTree_dual_DFS(queryNode->Right(), refNode->Left(), denThres, found);
							if (found)
								return;
							break;
						case 3:
							eps_connectness_BSTree_dual_DFS(queryNode->Right(), refNode->Right(), denThres, found);
							if (found)
								return;
							break;
						}
					}
				}
			}
			return;
		}
		// If we arrive here unfortunately, return false
		return;
	}
}

template <typename MetricType,
		  typename StatisticType,
		  typename MatType,
		  typename SplitType,
		  typename DescentType,
		  template <typename> class AuxiliaryInformationType>
void eps_connectness_RTree_dual_DFS(mlpack::RectangleTree<MetricType, StatisticType, MatType, SplitType, DescentType, AuxiliaryInformationType> &queryNode,
									mlpack::RectangleTree<MetricType, StatisticType, MatType, SplitType, DescentType, AuxiliaryInformationType> &refNode,
									// void eps_connectness_RTree_dual_DFS(mlpack::RectangleTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>& queryNode,
									//								mlpack::RectangleTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>& refNode,
									const double &denThres, bool &found)
{
	if (found)
		return;
	// If both are leaf, conduct BaseCase
	if (queryNode.IsLeaf() && refNode.IsLeaf())
	{
		double scoreLow = queryNode.MinDistance(refNode) - queryNode.FurthestDescendantDistance() - refNode.FurthestDescendantDistance();
		if (scoreLow > denThres)
			return;
		double scoreHigh = queryNode.MaxDistance(refNode) + queryNode.FurthestDescendantDistance() + refNode.FurthestDescendantDistance();
		if (scoreHigh <= denThres)
		{
			found = true;
			return;
		}
		int num1 = queryNode.NumPoints();
		int num2 = refNode.NumPoints();
		for (int i = 0; i < num1; i++)
		{
			for (int j = 0; j < num2; j++)
				if (MetricType::Evaluate
					// mlpack::EuclideanDistance::Evaluate
					(queryNode.Dataset().col(queryNode.Point(i)), refNode.Dataset().col(refNode.Point(j))) <= denThres)
				{
					found = true;
					return;
				}
		}
		return;
	}
	// if one is leaf, recurse the non-leaf one
	if ((queryNode.IsLeaf() && !refNode.IsLeaf()) || (!queryNode.IsLeaf() && refNode.IsLeaf()))
	{
		auto &leafOne = queryNode.IsLeaf() ? queryNode : refNode;
		auto &nonLeafOne = queryNode.IsLeaf() ? refNode : queryNode;

		for (size_t i = 0; i < nonLeafOne.NumChildren(); i++)
		{
			double scoreHigh = leafOne.MaxDistance(nonLeafOne.Child(i)) + leafOne.FurthestDescendantDistance() + nonLeafOne.Child(i).FurthestDescendantDistance();
			if (scoreHigh <= denThres)
			{
				found = true;
				return;
			}
		}
		std::vector<std::pair<size_t, double>> score_list(nonLeafOne.NumChildren());
		for (size_t i = 0; i < nonLeafOne.NumChildren(); i++)
		{
			score_list[i] = std::make_pair(i, leafOne.MinDistance(nonLeafOne.Child(i)) - leafOne.FurthestDescendantDistance() - nonLeafOne.Child(i).FurthestDescendantDistance());
		}
		std::sort(score_list.begin(), score_list.end(),
				  [](const std::pair<size_t, double> &p1, const std::pair<size_t, double> &p2)
				  {
					  return p2.second == p1.second ? p1.first < p2.first : p1.second < p2.second;
				  });
		for (size_t i = 0; i < score_list.size(); i++)
		{
			if (score_list[i].second > denThres)
				return;
			eps_connectness_RTree_dual_DFS(leafOne, nonLeafOne.Child(score_list[i].first), denThres, found);
			if (found)
				return;
		}
		return;
	}
	if (!queryNode.IsLeaf() && !refNode.IsLeaf())
	{
		for (size_t i = 0; i < queryNode.NumChildren(); i++)
		{
			for (size_t j = 0; j < refNode.NumChildren(); j++)
			{
				double scoreHigh = queryNode.Child(i).MaxDistance(refNode.Child(j)) + queryNode.Child(i).FurthestDescendantDistance() + refNode.Child(j).FurthestDescendantDistance();
				if (scoreHigh <= denThres)
				{
					found = true;
					return;
				}
			}
		}
		std::vector<std::tuple<size_t, size_t, double>> score_list(queryNode.NumChildren() * refNode.NumChildren());
		size_t index = 0;
		for (size_t i = 0; i < queryNode.NumChildren(); i++)
		{
			for (size_t j = 0; j < refNode.NumChildren(); j++)
			{
				score_list[index++] = std::make_tuple(i, j,
													  queryNode.Child(i).MinDistance(refNode.Child(j)) - queryNode.Child(i).FurthestDescendantDistance() - refNode.Child(j).FurthestDescendantDistance());
			}
		}
		std::sort(score_list.begin(), score_list.end(),
				  [](const std::tuple<size_t, size_t, double> &p1, const std::tuple<size_t, size_t, double> &p2)
				  {
					  return std::get<2>(p1) == std::get<2>(p2) ? std::get<0>(p1) < std::get<0>(p2) : std::get<2>(p1) < std::get<2>(p2);
				  });
		for (size_t i = 0; i < score_list.size(); i++)
		{
			if (std::get<2>(score_list[i]) > denThres)
				return;
			eps_connectness_RTree_dual_DFS(queryNode.Child(std::get<0>(score_list[i])), refNode.Child(std::get<1>(score_list[i])), denThres, found);
			if (found)
				return;
		}
		return;
	}
	return;
}
//
void eps_connectness_dual_CoverTree(mlpack::CoverTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat> &queryNode,
									mlpack::CoverTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat> &refNode,
									bool &found)
{
	double centerDist = mlpack::EuclideanDistance::Evaluate(queryNode.Dataset().col(queryNode.Point()), refNode.Dataset().col(refNode.Point()));
	if (centerDist <= denThres)
	{
		// Avoid duplicate base_cases
		found = true;
		return;
	}
	if (queryNode.IsLeaf() && refNode.IsLeaf())
		return;
	if (!queryNode.IsLeaf() && !refNode.IsLeaf())
	{
		size_t numq = queryNode.NumChildren(), numr = refNode.NumChildren();
		std::vector<std::tuple<size_t, size_t, double>> score_list;
		score_list.reserve((numq) * (numr));
		size_t recurse_index = 0;
		for (size_t i = 0; i < numq; i++)
		{
			// if (queryNode.Child(i).Point() == queryNode.Point())continue;
			size_t queryIndex = queryNode.Child(i).Point();
			for (size_t j = 0; j < numr; j++)
			{
				size_t refIndex = refNode.Child(j).Point();
				double dist;
				if (queryIndex == queryNode.Point() && refIndex == refNode.Point())
					dist = centerDist;
				// if (refNode.Child(j).Point() == refNode.Point())continue;
				else
					dist = mlpack::EuclideanDistance::Evaluate(queryNode.Dataset().col(queryIndex), refNode.Dataset().col(refIndex));
				if (dist <= denThres)
				{
					found = true;
					return;
				}
				// Seems no need to check the upper bound, since if the upper bound is within denThres, the base_case above should return true;
				// if (dist + queryNode.Child(i).FurthestDescendantDistance() + refNode.Child(j).FurthestDescendantDistance() <= denThres) {
				// found = true; return;
				//}
				double score = dist - queryNode.Child(i).FurthestDescendantDistance() - refNode.Child(j).FurthestDescendantDistance();
				if (score > denThres)
					continue;
				score_list[recurse_index++] = std::make_tuple(i, j, score);
			}
		}
		std::sort(score_list.begin(), score_list.end(),
				  [](const std::tuple<size_t, size_t, double> &p1, const std::tuple<size_t, size_t, double> &p2)
				  {
					  return std::get<2>(p1) == std::get<2>(p2) ? std::get<0>(p1) < std::get<0>(p2) : std::get<2>(p1) < std::get<2>(p2);
				  });
		for (size_t i = 0; i < score_list.size(); i++)
		{
			eps_connectness_dual_CoverTree(queryNode.Child(std::get<0>(score_list[i])), refNode.Child(std::get<1>(score_list[i])), found);
			if (found)
				return;
		}
		return;
	}
	else
	{ // one leaf, another not
		auto &leafOne = queryNode.IsLeaf() ? queryNode : refNode;
		auto &nonLeafOne = queryNode.IsLeaf() ? refNode : queryNode;
		size_t numc = nonLeafOne.NumChildren();
		std::vector<std::pair<size_t, double>> score_list;
		score_list.reserve(numc);
		size_t recurse_index = 0;
		for (size_t i = 0; i < numc; i++)
		{
			double dist;
			if (nonLeafOne.Child(i).Point() == nonLeafOne.Point())
				dist = centerDist;
			else
				dist = mlpack::EuclideanDistance::Evaluate(nonLeafOne.Dataset().col(nonLeafOne.Child(i).Point()), leafOne.Dataset().col(leafOne.Point()));
			if (dist <= denThres)
			{
				found = true;
				return;
			}
			// No need to check the upper bound
			// if (dist + nonLeafOne.Child(i).FurthestDescendantDistance() <= denThres) {//leafOne.FurthestDescendantDistance()=0
			// found = true; return;
			//}
			double score = dist - nonLeafOne.Child(i).FurthestDescendantDistance(); // leafOne.FurthestDescendantDistance()=0
			if (score > denThres)
				continue;
			score_list[recurse_index++] = std::make_pair(i, score);
		}
		std::sort(score_list.begin(), score_list.end(),
				  [](const std::pair<size_t, double> &p1, const std::pair<size_t, double> &p2)
				  {
					  return p2.second == p1.second ? p1.first < p2.first : p1.second < p2.second;
				  });
		for (size_t i = 0; i < score_list.size(); i++)
		{
			eps_connectness_dual_CoverTree(leafOne, nonLeafOne.Child(score_list[i].first), found);
			if (found)
				return;
		}
	}
}


void eps_connectness_dual_CoverTree_Lazy(CoverTree_Lazy& queryNode,
										 CoverTree_Lazy& refNode,
									bool& found) {
	
	double centerDist = mlpack::EuclideanDistance::Evaluate(queryNode.Dataset().col(queryNode.Point()), refNode.Dataset().col(refNode.Point()));
	if (centerDist <= denThres) {
		// Avoid duplicate base_cases
		found = true;
		return;
	}
	if (queryNode.IsLeaf() && refNode.IsLeaf())
		return;
	if (!queryNode.IsExpanded())queryNode.Expand();
	if (!refNode.IsExpanded()) refNode.Expand();
	if (!queryNode.IsLeaf() && !refNode.IsLeaf()) {
		size_t numq = queryNode.NumChildren(), numr = refNode.NumChildren();
		std::vector<std::tuple<size_t, size_t, double>> score_list;
		score_list.reserve((numq) * (numr));
		size_t recurse_index = 0;
		for (size_t i = 0; i < numq; i++) {
			// if (queryNode.Child(i).Point() == queryNode.Point())continue;
			size_t queryIndex = queryNode.Child(i).Point();
			for (size_t j = 0; j < numr; j++) {
				size_t refIndex = refNode.Child(j).Point();
				double dist;
				if (queryIndex == queryNode.Point() && refIndex == refNode.Point())
					dist = centerDist;
				// if (refNode.Child(j).Point() == refNode.Point())continue;
				else
					dist = mlpack::EuclideanDistance::Evaluate(queryNode.Dataset().col(queryIndex), refNode.Dataset().col(refIndex));
				if (dist <= denThres) {
					found = true;
					return;
				}
				// Seems no need to check the upper bound, since if the upper bound is within denThres, the base_case above should return true;
				// if (dist + queryNode.Child(i).FurthestDescendantDistance() + refNode.Child(j).FurthestDescendantDistance() <= denThres) {
				// found = true; return;
				//}
				double score = dist - queryNode.Child(i).FurthestDescendantDistance() - refNode.Child(j).FurthestDescendantDistance();
				if (score > denThres)
					continue;
				score_list[recurse_index++] = std::make_tuple(i, j, score);
			}
		}
		std::sort(score_list.begin(), score_list.end(),
				  [](const std::tuple<size_t, size_t, double>& p1, const std::tuple<size_t, size_t, double>& p2) {
					  return std::get<2>(p1) == std::get<2>(p2) ? std::get<0>(p1) < std::get<0>(p2) : std::get<2>(p1) < std::get<2>(p2);
				  });
		for (size_t i = 0; i < score_list.size(); i++) {
			eps_connectness_dual_CoverTree_Lazy(queryNode.Child(std::get<0>(score_list[i])), refNode.Child(std::get<1>(score_list[i])), found);
			if (found)
				return;
		}
		return;
	} else { // one leaf, another not
		auto& leafOne = queryNode.IsLeaf() ? queryNode : refNode;
		auto& nonLeafOne = queryNode.IsLeaf() ? refNode : queryNode;
		size_t numc = nonLeafOne.NumChildren();
		std::vector<std::pair<size_t, double>> score_list;
		score_list.reserve(numc);
		size_t recurse_index = 0;
		for (size_t i = 0; i < numc; i++) {
			double dist;
			if (nonLeafOne.Child(i).Point() == nonLeafOne.Point())
				dist = centerDist;
			else
				dist = mlpack::EuclideanDistance::Evaluate(nonLeafOne.Dataset().col(nonLeafOne.Child(i).Point()), leafOne.Dataset().col(leafOne.Point()));
			if (dist <= denThres) {
				found = true;
				return;
			}
			// No need to check the upper bound
			// if (dist + nonLeafOne.Child(i).FurthestDescendantDistance() <= denThres) {//leafOne.FurthestDescendantDistance()=0
			// found = true; return;
			//}
			double score = dist - nonLeafOne.Child(i).FurthestDescendantDistance(); // leafOne.FurthestDescendantDistance()=0
			if (score > denThres)
				continue;
			score_list[recurse_index++] = std::make_pair(i, score);
		}
		std::sort(score_list.begin(), score_list.end(),
				  [](const std::pair<size_t, double>& p1, const std::pair<size_t, double>& p2) {
					  return p2.second == p1.second ? p1.first < p2.first : p1.second < p2.second;
				  });
		for (size_t i = 0; i < score_list.size(); i++) {
			eps_connectness_dual_CoverTree_Lazy(leafOne, nonLeafOne.Child(score_list[i].first), found);
			if (found)
				return;
		}
	}
}


void eps_connectness_dual_CoverTree_Lazy_cendist(CoverTree_Lazy& queryNode,
										 CoverTree_Lazy& refNode,
										 double cenDist,bool& found) {
	if (queryNode.IsLeaf() && refNode.IsLeaf())
		return;
	
	double centerDist;
	if (queryNode.Parent() == nullptr || refNode.Parent() == nullptr)
		centerDist = mlpack::EuclideanDistance::Evaluate(queryNode.Dataset().col(queryNode.Point()), refNode.Dataset().col(refNode.Point()));
	else centerDist = cenDist;

	
	if (!queryNode.IsExpanded())queryNode.Expand();
	if (!refNode.IsExpanded()) refNode.Expand();
	if (!queryNode.IsLeaf() && !refNode.IsLeaf()) {
		size_t numq = queryNode.NumChildren(), numr = refNode.NumChildren();
		std::vector<std::tuple<size_t, size_t, double,double>> score_list;
		score_list.reserve((numq) * (numr));
		size_t recurse_index = 0;
		for (size_t i = 0; i < numq; i++) {
			// if (queryNode.Child(i).Point() == queryNode.Point())continue;
			size_t queryIndex = queryNode.Child(i).Point();
			for (size_t j = 0; j < numr; j++) {
				size_t refIndex = refNode.Child(j).Point();
				double dist;
				if (queryIndex == queryNode.Point() && refIndex == refNode.Point())
					dist = centerDist;
				// if (refNode.Child(j).Point() == refNode.Point())continue;
				else
					dist = mlpack::EuclideanDistance::Evaluate(queryNode.Dataset().col(queryIndex), refNode.Dataset().col(refIndex));
				if (dist <= denThres) {
					found = true;
					return;
				}
				// Seems no need to check the upper bound, since if the upper bound is within denThres, the base_case above should return true;
				// if (dist + queryNode.Child(i).FurthestDescendantDistance() + refNode.Child(j).FurthestDescendantDistance() <= denThres) {
				// found = true; return;
				//}
				double score = dist - queryNode.Child(i).FurthestDescendantDistance() - refNode.Child(j).FurthestDescendantDistance();
				if (score > denThres)
					continue;
				score_list[recurse_index++] = std::make_tuple(i, j, score,dist);
			}
		}
		std::sort(score_list.begin(), score_list.end(),
				  [](const std::tuple<size_t, size_t, double,double>& p1, const std::tuple<size_t, size_t, double,double>& p2) {
					  return std::get<2>(p1) == std::get<2>(p2) ? std::get<0>(p1) < std::get<0>(p2) : std::get<2>(p1) < std::get<2>(p2);
				  });
		for (size_t i = 0; i < score_list.size(); i++) {
			eps_connectness_dual_CoverTree_Lazy_cendist(queryNode.Child(std::get<0>(score_list[i])), refNode.Child(std::get<1>(score_list[i])), std::get<3>(score_list[i]), found);
			if (found)
				return;
		}
		return;
	} else { // one leaf, another not
		auto& leafOne = queryNode.IsLeaf() ? queryNode : refNode;
		auto& nonLeafOne = queryNode.IsLeaf() ? refNode : queryNode;
		size_t numc = nonLeafOne.NumChildren();
		std::vector<std::tuple<size_t, double,double>> score_list;
		score_list.reserve(numc);
		size_t recurse_index = 0;
		for (size_t i = 0; i < numc; i++) {
			double dist;
			if (nonLeafOne.Child(i).Point() == nonLeafOne.Point())
				dist = centerDist;
			else
				dist = mlpack::EuclideanDistance::Evaluate(nonLeafOne.Dataset().col(nonLeafOne.Child(i).Point()), leafOne.Dataset().col(leafOne.Point()));
			if (dist <= denThres) {
				found = true;
				return;
			}
			// No need to check the upper bound
			// if (dist + nonLeafOne.Child(i).FurthestDescendantDistance() <= denThres) {//leafOne.FurthestDescendantDistance()=0
			// found = true; return;
			//}
			double score = dist - nonLeafOne.Child(i).FurthestDescendantDistance(); // leafOne.FurthestDescendantDistance()=0
			if (score > denThres)
				continue;
			score_list[recurse_index++] = std::make_tuple(i, score,dist);
		}
		std::sort(score_list.begin(), score_list.end(),
				  [](const std::tuple<size_t, double,double>& p1, const std::tuple<size_t, double,double>& p2) {
					  return std::get<1>(p1) == std::get<1>(p2) ? std::get<0>(p1) < std::get<0>(p2) : std::get<1>(p1) < std::get<1>(p2);
				  });
		for (size_t i = 0; i < score_list.size(); i++) {
			eps_connectness_dual_CoverTree_Lazy_cendist(leafOne, nonLeafOne.Child(std::get<0>( score_list[i])), std::get<2>(score_list[i]), found);
			if (found)
				return;
		}
	}
}



// void eps_connectness_dual_CoverTree(mlpack::CoverTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>& queryNode,
//									mlpack::CoverTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>& refNode,
//									//void eps_connectness_RTree_dual_DFS(mlpack::RectangleTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>& queryNode,
//										//								mlpack::RectangleTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat>& refNode,
//									double cenDist, bool& found) {//Avoid duplicate base_cases
//	double centerDist;
//	if (queryNode.Parent() != NULL || refNode.Parent() != NULL)
//		centerDist = cenDist;
//	//If has parent, then this base case has bee already evaluated in the parent node.
//	else {
//		centerDist = mlpack::EuclideanDistance::Evaluate(queryNode.Dataset().col(queryNode.Point()), refNode.Dataset().col(refNode.Point()));
//		if (centerDist <= denThres) {
//			found = true; return;
//		}
//	}
//	if (queryNode.IsLeaf() && refNode.IsLeaf()) return;
//	if (!queryNode.IsLeaf() && !refNode.IsLeaf()) {
//		size_t numq = queryNode.NumChildren(), numr = refNode.NumChildren();
//		std::vector<std::tuple<size_t, size_t, double>> score_list;
//		score_list.reserve((numq) * (numr)); size_t recurse_index = 0;
//		double* distmat = new double(numq * numr);
//		for (size_t i = 0; i < numq; i++) {
//			//if (queryNode.Child(i).Point() == queryNode.Point())continue;
//			size_t queryIndex = queryNode.Child(i).Point();
//			for (size_t j = 0; j < numr; j++) {
//				size_t refIndex = refNode.Child(j).Point();
//				double dist;
//				if (queryIndex == queryNode.Point() && refIndex == refNode.Point())
//					dist = centerDist;
//				//if (refNode.Child(j).Point() == refNode.Point())continue;
//				else dist = mlpack::EuclideanDistance::Evaluate(queryNode.Dataset().col(queryIndex), refNode.Dataset().col(refIndex));
//				if (dist <= denThres) {
//					found = true;  return;
//				}
//				//Seems no need to check the upper bound, since if the upper bound is within denThres, the base_case above should return true;
//				//if (dist + queryNode.Child(i).FurthestDescendantDistance() + refNode.Child(j).FurthestDescendantDistance() <= denThres) {
//					//found = true; return;
//				//}
//				double score = dist - queryNode.Child(i).FurthestDescendantDistance() - refNode.Child(j).FurthestDescendantDistance();
//				if (score > denThres)continue;
//				score_list[recurse_index++] = std::make_tuple(i, j, score);
//				distmat[i*numq+j] = dist;
//			}
//		}
//		std::sort(score_list.begin(), score_list.end(),
//				  [](const std::tuple<size_t, size_t, double>& p1, const std::tuple<size_t, size_t, double>& p2) {
//					  return std::get<2>(p1) == std::get<2>(p2) ? std::get<0>(p1) < std::get<0>(p2) : std::get<2>(p1) < std::get<2>(p2);
//				  });
//		for (size_t i = 0; i < score_list.size(); i++) {
//			size_t qind = std::get<0>(score_list[i]), rind = std::get<1>(score_list[i]);
//			eps_connectness_dual_CoverTree(queryNode.Child(qind), refNode.Child(rind),  distmat[qind*numq+rind], found);
//			if (found) {
//				return;
//			}
//		}
//
//		return;
//	} else {//one leaf, another not
//		auto& leafOne = queryNode.IsLeaf() ? queryNode : refNode;
//		auto& nonLeafOne = queryNode.IsLeaf() ? refNode : queryNode;
//		size_t numc = nonLeafOne.NumChildren();
//		std::vector<std::pair<size_t, double>> score_list;
//		score_list.reserve(numc); size_t recurse_index = 0;
//		double* distvec = new double[numc];
//		for (size_t i = 0; i < numc; i++) {
//			double dist;
//			if (nonLeafOne.Child(i).Point() == nonLeafOne.Point())
//				dist = centerDist;
//			else dist = mlpack::EuclideanDistance::Evaluate(nonLeafOne.Dataset().col(nonLeafOne.Child(i).Point()), leafOne.Dataset().col(leafOne.Point()));
//			if (dist <= denThres) {
//				found = true; return;
//			}
//			//No need to check the upper bound
//			//if (dist + nonLeafOne.Child(i).FurthestDescendantDistance() <= denThres) {//leafOne.FurthestDescendantDistance()=0
//				//found = true; return;
//			//}
//			double score = dist - nonLeafOne.Child(i).FurthestDescendantDistance();//leafOne.FurthestDescendantDistance()=0
//			if (score > denThres)continue;
//			score_list[recurse_index++] = std::make_pair(i, score);
//			distvec[i] = dist;
//		}
//		std::sort(score_list.begin(), score_list.end(),
//				  [](const std::pair<size_t, double>& p1, const std::pair<size_t, double>& p2) {
//					  return p1.second == p2.second ? p1.first < p2.first : p1.second < p2.second;
//				  });
//		for (size_t i = 0; i < score_list.size(); i++) {
//			eps_connectness_dual_CoverTree(leafOne, nonLeafOne.Child(score_list[i].first),  distvec[score_list[i].first], found);
//			if (found) return;
//		}
//	}
// }
