#pragma once
/**
 * Statistic class for EpsMinptsQuery. to be set to the StatisticType of the tree type that 
 * eps_minpts query is being performed with.
 * eps_minpts query is a truncted range query, where the range search is truncated 
 * when all query points in the query node are filled with minpts neighbors
 * Thus, this class holds the last visited node and the corresponding base case result, 
 * which is the same with the RangeSearchStat class.
 * Besides, this class holds also a boolean value isAllDecided, 
 * indicating whether the points in the query nodes are filled with minpts neighbors
 */

class EpsMinptsStat {
public:
	/**
	 * Initialize the statistic.
	 */
	EpsMinptsStat() : lastDistance(0.0), isAllDecided(false) {}

	/**
	 * Initialize the statistic given a tree node that this statistic belongs to.
	 * In this case, we ignore the node.
	 */
	template <typename TreeType>
	EpsMinptsStat(TreeType& /* node */) : lastDistance(0.0), isAllDecided(false) {}

	//! Get the last distance evaluation.
	double LastDistance() const {
		return lastDistance;
	}
	//! Modify the last distance evaluation.
	double& LastDistance() {
		return lastDistance;
	}
	bool IsAllDecided() const {
		return isAllDecided;
	}
	bool& IsAllDecided() {
		return isAllDecided;
	}


private:
	//! The last distance evaluation.
	double lastDistance;
	//! The boolean flag indicating whether all the points in this node have already found minpts number of neighbors.
	bool isAllDecided;
};