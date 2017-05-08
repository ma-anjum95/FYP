#include <iostream>
#include <stdio.h>
#include <wiringPi.h>
#include <bcm2835.h>
#include <time.h>
#include <unistd.h>
#include <fstream>
#include <vector>
#include <thread>

#include <chrono>

using namespace std::chrono;

#include "PPG_C++/ppg_analysis.h"
#include "MAX117/max30102.h"

using namespace std;
using namespace PPG;

double *linear_interp_10(vector<double> to_interp, int start_index);

int main() 
{
	cout << "Please wait 20 seconds before the values start to display." << endl;
	cout << "Keep your finger stable and do not talk during the session." << endl;
	uint32_t tmp1, tmp2;
	signed int i = 0, update=125;
	PPGAnalysis ppg_analysis;
	vector<double> ppg_red, ppg_ir;
	double *interp_red, *interp_ir;
	if (!maxim_max30102_init())
		cout << "failure initializing max30102" << endl;
	int t = 1;

	uint8_t cr1, cr2;
	while(ppg_red.size() < 25*60*2) {
		if (maxim_max30102_read_fifo(&tmp1, &tmp2)) {
		
			ppg_red.push_back((double)tmp1);
			ppg_ir.push_back((double)tmp2);

			while (i < (signed) ppg_red.size() - 500) {
				cout << "############## " << i/25 << " ###############" << endl;
				interp_red = linear_interp_10(ppg_red, i);
				interp_ir = linear_interp_10(ppg_ir, i);

				ppg_analysis.run(interp_red, interp_ir, 5000, 250);
				cout << "hr: " << ppg_analysis.get_hr() << endl;
				cout << "rr: " << ppg_analysis.get_rr() << endl;
				cout << "rr std: " << ppg_analysis.get_rr_std() << endl;
				cout << "spo2: " << ppg_analysis.get_spo2() << endl;


				delete interp_red;
				delete interp_ir;
				i += update;
			}
		}
	}
 
}

double *linear_interp_10(vector<double> to_interp, int start_index)
{
	// The formula used is:
	// y = (x - x1) (y2 - y1)
	//              ========= + y1
	//              (x2 - x1)
	// we can assume x starts at zero every time so
	// x - x1 = j
	// x2 -x1 = 10

	// be sure to delete it somewhere
	// it will always return an array of 5000. 
	// incase the samples are less than 500 we return an error
	double *to_return = new double[5000];

	for (int i = start_index; i < 500 + start_index; i++) {
		for (int j = 0; j < 10; j++) {
			int new_index = (i - start_index) * 10 + j;
			
			to_return[new_index] = ((double)j / 10.0) * (to_interp[i + 1] - to_interp[i]) + to_interp[i];
		}
	}

	return to_return;
}
