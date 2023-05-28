#pragma once


CoverTree_Lazy::CoverTree_Lazy(
	const arma::mat& dataset,
	const double base) :
	dataset(&dataset),
	point(0),
	base(base),
	scale(INT_MAX),
	numDescendants(dataset.n_cols),
	isExpanded(false),
	parent(NULL)

{//First point is root
	for (size_t i = 1; i < dataset.n_cols; i++) {
		double dist = mlpack::EuclideanDistance::Evaluate(dataset.col(0), dataset.col(i));
		add_descendant(i, dist);
	}
	setup_info();
	//Expand the root on initialization. May be uncessary.
	//Expand();
}
CoverTree_Lazy::CoverTree_Lazy(
	const arma::mat& dataset,
	size_t index,
	int scale,
	CoverTree_Lazy* parent,
	const double base) :
	dataset(&dataset),
	point(index),
	scale(scale),
	parent(parent),
	base(base),
	isExpanded(false) {

}
void CoverTree_Lazy::setup_info() {
	scale = ceil(log(furthestDescendantDistance) / log(base));
	numDescendants = descendants.size();
}
#ifdef CTL_USE_SET
void CoverTree_Lazy::add_descendant(size_t index, double dist) {
	//descendants.push_back(index);
	descendants.insert(index);
	furthestDescendantDistance = dist > furthestDescendantDistance ? dist : furthestDescendantDistance;
}
void CoverTree_Lazy::add_child(size_t index) {
	CoverTree_Lazy* child = new CoverTree_Lazy(*dataset, index, scale - 1, this, base);
	children.push_back(child);
	for (auto iter = descendants.begin(); iter != descendants.end();) {
		double dist = mlpack::EuclideanDistance::Evaluate(dataset->col(index), dataset->col(*iter));
		if (dist <= std::pow(base, scale - 1)) {
			child->add_descendant(*iter, dist);
			//desc_distances.erase(*iter);
			iter = descendants.erase(iter);
		} else iter++;
	}
	child->setup_info();
}
void CoverTree_Lazy::Expand() {
	//Since we perform full expand for each expand, we do not need to delete an arbitrary element from sescendants. 
	//Thus, it need not to be a std::set, std::vector is enough
	isExpanded = true;
	//Add self-child first.
	add_child(point);
	//Loop to add the others.
	while (!descendants.empty()) {
		size_t index = *descendants.begin();//A better Heuristic need
		descendants.erase(descendants.begin());
		add_child(index);
	}
}
#endif
#ifdef CTL_USE_VECTOR
void CoverTree_Lazy::add_descendant(size_t index, double dist) {
	descendants.push_back(index);
	furthestDescendantDistance = dist > furthestDescendantDistance ? dist : furthestDescendantDistance;
}
void CoverTree_Lazy::add_child(size_t index, size_t& pointLeft) {
	CoverTree_Lazy* child = new CoverTree_Lazy(*dataset, index, scale - 1, this, base);
	children.push_back(child);
	for (size_t iter = 0; iter < pointLeft;) {
		double dist = mlpack::EuclideanDistance::Evaluate(dataset->col(index), dataset->col(descendants[iter]));
		if (dist <= std::pow(base, scale - 1)) {
			child->add_descendant(descendants[iter], dist);
			//desc_distances.erase(*iter);
			std::swap(descendants[iter], descendants[pointLeft]);
			pointLeft--;
		} else iter++;
	}
	child->setup_info();
}
void CoverTree_Lazy::Expand() {
	//Since we perform full expand for each expand, we do not need to delete an arbitrary element from sescendants. 
	//Thus, it need not to be a std::set, std::vector is enough
	isExpanded = true;
	size_t pointLeft = numDescendants - 1;
	//Add self-child first.
	add_child(point, pointLeft);
	//Loop to add the others.
	while (pointLeft > 0) {
		size_t index = descendants[0];//A better Heuristic need
		std::swap(descendants[0], descendants[pointLeft]);
		pointLeft--;
		add_child(index, pointLeft);
	}
	
}




#endif