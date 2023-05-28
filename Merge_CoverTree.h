#pragma once

//typedef mlpack::CoverTree<mlpack::EuclideanDistance, mlpack::EmptyStatistic, arma::mat> CoverTreeHere;
template<typename MetricType, typename StatType, typename MatType>
mlpack::CoverTree<MetricType, StatType, MatType>* merge_covertree(mlpack::CoverTree<MetricType, StatType, MatType>* root1, mlpack::CoverTree<MetricType, StatType, MatType>* root2) {
	if (root2->Scale() > root1->Scale())std::swap(root1, root2);
	double dis = mlpack::EuclideanDistance::Evaluate(root1->Dataset().col(root1->Point()), root2->Dataset().col(root2->Point()));
	if (dis <= pow(root1->Base(), root1->Scale())) {
		//The cover_dis of root1 is greater than dis, then node1 can cover root2
		//Let root2 be a child of root1
		//TODO! :update the information of the two nodes if needed
		root1->Children().push_back(root2);
		//!! Not correct. furthestDescendant, not furthest Child
		root1->FurthestDescendantDistance() = dis > root1->FurthestDescendantDistance() ? dis : root1->FurthestDescendantDistance();
		//!!-----------------------------
		root2->Parent() = root1;
		root2->ParentDistance() = dis;
		//------------------------
		return root1;
	} else {
		int new_scale = std::ceil(log(dis) / log(root1->Base()));
		//Format of the constructor
		// CoverTree 	( 	const MatType &  	dataset,
		//    const ElemType  	base,
		//    const size_t  	pointIndex,
		//    const int  	scale,
		//    CoverTree< MetricType, StatisticType, MatType, RootPointPolicy >* parent,
		//    const ElemType  	parentDistance,
		//    const ElemType  	furthestDescendantDistance,
		//    MetricType* metric=NULL
		//    )
		//Constuct a new node with scale =new_scale, let the rep_point of the new node be node1.Point
		//let node1 and node2 be the children of the new_node;

		mlpack::CoverTree<MetricType, StatType, MatType>* new_root = new mlpack::CoverTree<MetricType, StatType, MatType>(root1->Dataset(), root1->Base(), root1->Point(), new_scale, NULL, 0, dis);

		new_root->Children().push_back(root1);
		new_root->Children().push_back(root2);
		return new_root;
	}
}

CoverTree_Lazy* merge_covertree(CoverTree_Lazy* root1, CoverTree_Lazy* root2) {
	if (root2->Scale() > root1->Scale())std::swap(root1, root2);
	double dis = mlpack::EuclideanDistance::Evaluate(root1->Dataset().col(root1->Point()), root2->Dataset().col(root2->Point()));
	if (dis <= pow(root1->Base(), root1->Scale())) {
		//The cover_dis of root1 is greater than dis, then node1 can cover root2
		//Let root2 be a child of root1
		//TODO! :update the information of the two nodes if needed
		root1->Children().push_back(root2);
		//!! Not correct. furthestDescendant, not furthest Child
		root1->FurthestDescendantDistance() = dis > root1->FurthestDescendantDistance() ? dis : root1->FurthestDescendantDistance();
		//!!-----------------------------
		root2->Parent() = root1;
		
		//------------------------
		return root1;
	} else {
		int new_scale = std::ceil(log(dis) / log(root1->Base()));
		//Format of the constructor
		// CoverTree 	( 	const MatType &  	dataset,
		//    const ElemType  	base,
		//    const size_t  	pointIndex,
		//    const int  	scale,
		//    CoverTree< MetricType, StatisticType, MatType, RootPointPolicy >* parent,
		//    const ElemType  	parentDistance,
		//    const ElemType  	furthestDescendantDistance,
		//    MetricType* metric=NULL
		//    )
		//Constuct a new node with scale =new_scale, let the rep_point of the new node be node1.Point
		//let node1 and node2 be the children of the new_node;

		CoverTree_Lazy* new_root = new CoverTree_Lazy(root1->Dataset(),  root1->Point(), new_scale, nullptr,root1->Base());

		new_root->Children().push_back(root1);
		new_root->Children().push_back(root2);
		return new_root;
	}
}

//
//template<typename MetricType, typename StatType, typename MatType>
//CoverTree_ExpandOnQuery<MetricType, StatType, MatType>* merge_covertree(CoverTree_ExpandOnQuery<MetricType, StatType, MatType>* root1, CoverTree_ExpandOnQuery<MetricType, StatType, MatType>* root2) {
//	if (root2->Scale() > root1->Scale())std::swap(root1, root2);
//	double dis = mlpack::EuclideanDistance::Evaluate(root1->Dataset().col(root1->Point()), root2->Dataset().col(root2->Point()));
//	if (dis <= pow(root1->Base(), root1->Scale())) {
//		//The cover_dis of root1 is greater than dis, then node1 can cover root2
//		//Let root2 be a child of root1
//		//TODO! :update the information of the two nodes if needed
//		root1->Children().push_back(root2);
//		//!! Not correct. furthestDescendant, not furthest Child
//		root1->FurthestDescendantDistance() = dis > root1->FurthestDescendantDistance() ? dis : root1->FurthestDescendantDistance();
//		//!!-----------------------------
//		root2->Parent() = root1;
//		root2->ParentDistance() = dis;
//		//------------------------
//		return root1;
//	} else {
//		int new_scale = std::ceil(log(dis) / log(root1->Base()));
//		//Format of the constructor
//		// CoverTree 	( 	const MatType &  	dataset,
//		//    const ElemType  	base,
//		//    const size_t  	pointIndex,
//		//    const int  	scale,
//		//    CoverTree< MetricType, StatisticType, MatType, RootPointPolicy >* parent,
//		//    const ElemType  	parentDistance,
//		//    const ElemType  	furthestDescendantDistance,
//		//    MetricType* metric=NULL
//		//    )
//		//Constuct a new node with scale =new_scale, let the rep_point of the new node be node1.Point
//		//let node1 and node2 be the children of the new_node;
//
//		CoverTree_ExpandOnQuery<MetricType, StatType, MatType>* new_root = new CoverTree_ExpandOnQuery<MetricType, StatType, MatType>(root1->Dataset(), root1->Base(), root1->Point(), new_scale, NULL, 0, dis);
//
//		new_root->Children().push_back(root1);
//		new_root->Children().push_back(root2);
//		return new_root;
//	}
//}
