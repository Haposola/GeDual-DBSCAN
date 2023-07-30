#pragma once


CoverTree_Lazy::CoverTree_Lazy(
	const arma::mat& dataset,
	const double base) :
	dataset(&dataset),
	point(0),
	base(base),
	scale(INT_MAX),
	numDescendants(dataset.n_cols - 1),
	isExpanded(false),
	parent(NULL)

{//First point is root. No self-child is added here.
	for (size_t i = 1; i < dataset.n_cols; i++) {
		double dist = mlpack::EuclideanDistance::Evaluate(dataset.col(0), dataset.col(i));
		add_descendant(i, dist);
	}
	setup_info();
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
	furthestDescendantDistance(0),
	numDescendants(0),
	base(base),
	isExpanded(false) {

}
void CoverTree_Lazy::setup_info() {
	scale = ceil(log(furthestDescendantDistance) / log(base));
	numDescendants = descendants.size();
}

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
			
			/*{
				size_t temp = descendants[iter]; descendants[iter] = descendants[pointLeft - 1]; descendants[pointLeft - 1] = temp;
			}*/

			std::swap(descendants[iter], descendants[pointLeft - 1]);
			pointLeft--;
		} else iter++;
	}
	child->setup_info();
}
void CoverTree_Lazy::Expand() {
	//Since we perform full expand for each expand, we do not need to delete an arbitrary element from descendants. 
	//Thus, it need not to be a std::set, std::vector is enough
	isExpanded = true;
	//pointLeft is index of the first invalid element.
	size_t pointLeft = numDescendants;
	//Loop to add the others.
	while (pointLeft > 0) {
		size_t index = descendants[0];
		std::swap(descendants[0], descendants[pointLeft - 1]);//correct, and self-child is avoided too.
		/*{
			size_t temp = descendants[0]; descendants[0] = descendants[pointLeft - 1]; descendants[pointLeft - 1] = temp;
		}*/
		pointLeft--;
		add_child(index, pointLeft);
	}

}



