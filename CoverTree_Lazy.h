#pragma once
//#define CTL_USE_SET
#define CTL_USE_VECTOR
class CoverTree_Lazy {
public:
	CoverTree_Lazy(const arma::mat& dataset,
				   const double base = 2.0);
	CoverTree_Lazy(const arma::mat& dataset,
				   size_t index,
				   int scale,
				   CoverTree_Lazy* parent,
				   const double base = 2.0);
	//! Get a reference to the dataset.
	const arma::mat& Dataset() const {
		return *dataset;
	}

	//! Get the index of the point which this node represents. 
	size_t Point() const {
		return point;
	}
	//! For compatibility with other trees; the argument is ignored.
	size_t Point(const size_t) const {
		return point;
	}

	bool IsLeaf() const {
		return (children.size() == 0);
	}
	size_t NumPoints() const {
		return 1;
	}

	//! Get a particular child node.
	const CoverTree_Lazy& Child(const size_t index) const {
		return *children[index];
	}
	//! Modify a particular child node.
	CoverTree_Lazy& Child(const size_t index) {
		return *children[index];
	}

	CoverTree_Lazy*& ChildPtr(const size_t index) {
		return children[index];
	}

	//! Get the number of children.
	size_t NumChildren() const {
		return children.size();
	}

	//! Get the children.
	const std::vector<CoverTree_Lazy*>& Children() const {
		return children;
	}
	//! Modify the children manually (maybe not a great idea).
	std::vector<CoverTree_Lazy*>& Children() {
		return children;
	}

	//! Get the number of descendant points.
	size_t NumDescendants() const {
		return numDescendants;
	};

	//! Get the index of a particular descendant point.
	size_t Descendant(const size_t index) const;

	//! Get the scale of this node.
	int Scale() const {
		return scale;
	}
	//! Modify the scale of this node.  Be careful...
	int& Scale() {
		return scale;
	}

	//! Get the base.
	double Base() const {
		return base;
	}
	//! Modify the base; don't do this, you'll break everything.
	double& Base() {
		return base;
	}
	//! Get the parent node.
	CoverTree_Lazy* Parent() const {
		return parent;
	}
	//! Modify the parent node.
	CoverTree_Lazy*& Parent() {
		return parent;
	}
	double FurthestDescendantDistance() const {
		return furthestDescendantDistance;
	}
	double& FurthestDescendantDistance() {
		return furthestDescendantDistance;
	}
	bool IsExpanded()const {
		return isExpanded;
	}
	bool& IsExpanded() {
		return isExpanded;
	}
	void Expand();
	void add_descendant(size_t index, double dist);
	void setup_info();
private:


	std::vector<size_t> descendants;


	//std::unordered_map<size_t, double> desc_distances;
	bool isExpanded;

	//! Reference to the matrix which this tree is built on.
	const arma::mat* dataset;
	//! Index of the point in the matrix which this node represents.
	size_t point;
	//! The list of children; the first is the self-child.
	std::vector<CoverTree_Lazy*> children;
	//! Scale level of the node.
	int scale;
	//! The base used to construct the tree.
	double base;

	//! The number of descendant points.
	size_t numDescendants;
	//! The parent node (NULL if this is the root of the tree).
	CoverTree_Lazy* parent;

	//! Distance to the furthest descendant.
	double furthestDescendantDistance;

	//void compute_distances();
	//void delete_from_list(size_t index);
#ifdef CTL_USE_SET
	void add_child(size_t index);
#endif
#ifdef CTL_USE_VECTOR
	void add_child(size_t index, size_t& pointLeft);
#endif


	//CoverTree_Lazy* create_childptr(size_t index);
};
#include "CoverTree_Lazy_impl.h"