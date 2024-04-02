#pragma once
void GeDual_DBSCAN::PointTransite_ToCore(size_t point) {
	isCore[point] = true; 
	for (size_t refNbr : epsNeighbors[point]) {
		if (isCore[refNbr]) ufset.Union(point, refNbr);
	}
}