#include <iostream>
#include <stdio.h>
#include <wiringPi.h>
#include <bcm2835.h>
#include "max30102.h"
#include <time.h>
#include <unistd.h>
#include <fstream>
#include <vector>

using namespace std;

// -I/usr/local/include -L/usr/local/lib -lwiringPi

vector<double> linear_intep(vector<double> to_interp, int start_index, int end_index);

int main(void) 
{	
	int i = 0; // index for the start of ppg array
	int window = 500;
	int start_index, end_index; // index used in extracting a specific window
	vector<double> interpolated;
	uint32_t tmp1, tmp2;
	vector<double> ppg_data_red;
	vector<double> ppg_data_ir;
	
	ofstream ofs("data.txt", fstream::out);
	ofstream ofs_interp("data_interp.txt", fstream::out);
	
	
	if (!maxim_max30102_init())
		printf("failure initializing max30102. \n");
	
	printf("Initialzing the sensor. Place finger on the sensor head.");
	sleep(2000);

	
	while(ppg_data_red.size() < 25 * 60 * 5) {

		if (maxim_max30102_read_fifo(&tmp1, &tmp2)) {
			ppg_data_red.push_back((double)tmp1);
			ppg_data_ir.push_back((double)tmp2);
		}
		
		if (ppg_data_red.size() >= (i + 1) * window ) {
			start_index = i * window;
			end_index = (i + 1) * window; // this index is like length. its own value will not be included
			
			vector<double> interpolated_red = linear_intep(ppg_data_red, start_index, end_index);
			vector<double> interpolated_ir = linear_intep(ppg_data_ir, start_index, end_index);
			
			for (int j = 0, n = interpolated_red.size(); j < n; j++)
				ofs_interp << interpolated_red[j] << " " << interpolated_ir[j] << endl;
			
			i++;
		}
	}
	
	for (int j = 0, n = ppg_data_red.size(); j < n; j++)
		ofs << ppg_data_red[i] << " " << ppg_data_ir[i] << endl;
	
	return true;
	
}

vector<double> linear_intep(vector<double> to_interp, int start_index, int end_index)
{
	// interpolates the to_interp vector by a factor of 10
	// since our device operates at 25Hz we want our data to be at 250Hz
	vector<double> to_return;
	double y_diff, x_diff;
	// we will be using two point formula
	// y = (y1 - y0)(x - x0) 
	//     ================= + y0
	//         (x1 - x0)
	for (int i = start_index, n = end_index - 1; i < n; i++) {
			//double y0 = to_interp[i];
			//double y1 = to_interp[i + 1];
			y_diff = to_interp[i + 1] - to_interp[i];
			
			//int x0 = 0;
			//int x1 = 10;
			x_diff = 10.0;

		for (int j = 0; j < 10; j++) {
			double y_interp = (y_diff/x_diff) * j + to_interp[i];

			to_return.push_back(y_interp);
		}
	}	

	return to_return;
}
