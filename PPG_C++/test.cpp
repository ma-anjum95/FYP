#include <iostream>
#include <fstream>

#include "ppg_analysis.h"

using namespace std;

int main(void) 
{
	double *arr = new double[6000];
	int samp_freq = 300;
	int window_size = 6000;

	PPG::PPGAnalysis ppg_analysis;
	ifstream file("data.txt");

	
	if (!file.is_open()) {
		cout << "File not opened" << endl;
		return 1;
	}

	for (int i = 0; i < window_size; i++) {
		file >> arr[i];
	}


	ppg_analysis.run(arr, arr, window_size, samp_freq);
	cout << ppg_analysis.get_hr() << " " << ppg_analysis.get_rr() << " " << ppg_analysis.get_spo2() << endl;


	ppg_analysis.run(arr, arr, window_size, samp_freq);
	cout << ppg_analysis.get_hr() << " " << ppg_analysis.get_rr() << " " << ppg_analysis.get_spo2() << endl;

	int tmp;
	cin >> tmp;
}