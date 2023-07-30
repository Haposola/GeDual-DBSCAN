#pragma once

CoverTree_Lazy* merge_covertree(CoverTree_Lazy* root1, CoverTree_Lazy* root2) {

	if (root2->Scale() > root1->Scale()) {
		double dis = mlpack::EuclideanDistance::Evaluate(root1->Dataset().unsafe_col(root1->Point()), root2->Dataset().unsafe_col(root2->Point()));
		double f1 = root2->FurthestDescendantDistance();
		double f2 = dis + root1->FurthestDescendantDistance();
		double fm = f1 > f2 ? f1 : f2;
		if (dis > pow(root2->Base(), root2->Scale())) {
			int new_scale = std::ceil(log(fm) / log(root2->Base()));
			root2->Scale() = new_scale;
		}

		root2->Children().push_back(root1);
		root2->FurthestDescendantDistance() = fm;

		root1->Parent() = root2;
		return root2;


	}else{
		double dis = mlpack::EuclideanDistance::Evaluate(root1->Dataset().unsafe_col(root1->Point()), root2->Dataset().unsafe_col(root2->Point()));
		double f1 = root1->FurthestDescendantDistance();
		double f2 = dis + root2->FurthestDescendantDistance();
		double fm = f1 > f2 ? f1 : f2;
		if (dis > pow(root1->Base(), root1->Scale())) {
			int new_scale = std::ceil(log(fm) / log(root1->Base()));
			root1->Scale() = new_scale;
		}

		root1->Children().push_back(root2);
		root1->FurthestDescendantDistance() = fm;

		root2->Parent() = root1;
		return root1;
	
	
	}
	


}