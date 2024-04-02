#pragma once
void test_exptime() {
	auto t1 = ExperimentRun_TimeNow();
	int j = 0;
	for (int i = 0; i < 10000000; i++) j += i*100;
	auto t2 = ExperimentRun_TimeNow();
	std::cout << ExperimentRun_TimeCount_s_Seconds(t1, t2)<<" hhh\n"		;
	std::cout << ExperimentRun_TimeCount_ms_MicroSeconds(t1, t2) << " hhh\n";


}
void generate_data(arma::mat& data, int num, int dim) {
	//setsize(n_rows,n_cols)
	data.set_size(dim, num);
	//Generate random data
	std::srand((int)time(0));
	double range = num * 10 > INT_MAX ? INT_MAX : num * 10;
	//arma::mat is column-wise. a column is a data.
	for (int i = 0; i < num; i++) {
		for (int j = 0; j < dim; j++) {
			data(j, i) = rand() * 1.0 / RAND_MAX * range;
			data(j, i) += rand() * 1.0 / RAND_MAX;
		}
	}
}

