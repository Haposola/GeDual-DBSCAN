#pragma once
#include <chrono>


#define ExperimentRun_TimeNow() std::chrono::steady_clock::now()
#define ExperimentRun_TimeCount_s_Seconds(t1,t2)  (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count())*1.0/1000000
#define ExperimentRun_TimeCount_ms_MicroSeconds(t1,t2)  (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count())

class StopW {
	std::chrono::steady_clock::time_point time_begin;
public:
	StopW() {
		time_begin = std::chrono::steady_clock::now();
	}

	double getElapsedTimeMicro() {
		std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
		return (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_begin).count());
	}
	
	void reset() {
		time_begin = std::chrono::steady_clock::now();
	}

};