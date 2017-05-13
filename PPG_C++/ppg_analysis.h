
#ifndef PPG_ANALYSIS_H
#define PPG_ANALYSIS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#include "kiss_fft.h"

#include "ppg_lines.h"
#include "point.h"

using namespace std;

namespace PPG 
{
	class PPGAnalysis 
	{
		int samp_freq;
		int window_size;

		double *ppg_red_hpf, *ppg_ir_hpf;

		double *ppg_red_ac, *ppg_ir_ac;
		double ppg_red_dc, ppg_ir_dc;

		Point ppg_red_intensity_fft_max;
		Point ppg_ir_intensity_fft_max;

		Point ppg_red_amplitude_fft_max;
		Point ppg_ir_amplitude_fft_max;
		
		vector<PPGLines> *ppg_red_lines, *ppg_ir_lines;
		vector<PPGLines> *ppg_red_lines_processed, *ppg_ir_lines_processed;
		vector<Point> *ppg_red_intensity_waveform, *ppg_ir_intensity_waveform;
		vector<Point> *ppg_red_amplitude_waveform, *ppg_ir_amplitude_waveform;
		vector<double> *ppg_red_intensity_waveform_interp, *ppg_ir_intensity_waveform_interp;
		vector<double> *ppg_red_amplitude_waveform_interp, *ppg_ir_amplitude_waveform_interp;
		vector<Point> *ppg_red_intensity_waveform_fft, *ppg_ir_intensity_waveform_fft;
		vector<Point> *ppg_red_amplitude_waveform_fft, *ppg_ir_amplitude_waveform_fft;

		// Results of the algorithm
		double hr;
		double spo2;
		double rr;
		double rr_std;

		// cleans the contents of the vectors
		void clean_vectors(void);

		// high pass filter
		void ppg_hpf(double *ppg_raw, double *&ppg_hpf);

		// dc filter
		void ppg_dcf(double *ppg_hpf, double *&ppg_ac, double &ppg_dc);

		// ppg line segementation algorithm
		void ppg_line_segmentation(vector<PPGLines> *to_return, double *ppg_ac, int m = 10);

		// ppg artifact removal algorithm
		void ppg_artifact_removal(vector<PPGLines> *to_return, vector<PPGLines> *ppg_lines);

		// ppg intensity waveform
		void ppg_intensity_waveform(vector<Point> *to_return, vector<PPGLines> *ppg_lines_processed);

		// ppg amplitude waveform
		void ppg_amplitude_waveform(vector<Point> *to_return, vector<PPGLines> *ppg_lines_processed);

		// linear interpolation function which converts old freq to new
		void linear_interpolation_with_freq(vector<double> *to_return, vector<Point> *data, int old_freq, int new_freq);

		// extract the heart rate from the intensity waveform
		// we will be using the ir waveform since it is much clearer and noise free
		double ppg_heart_rate(vector<PPGLines> *ppg_lines_processed);

		// extracting the spo2 from both the signal
		double ppg_spo2(double *ppg_red_ac, double ppg_red_dc, double *ppg_ir_ac, double ppg_ir_dc);

		// helper function for calculating the ppg
		double calculate_rms(double *ac_signal, int length);

		// caculates the fft of the waveform signals
		void ppg_waveform_fft(vector<Point> *to_return, vector<double> *waveform, int freq);

		// helper function for finding maximum in an array within a given range
		Point find_max_in_range(vector<Point> *arr, double low, double high);

		// calculates respiratory rate from the max freq from fft
		// The x in point is the rr
		// The y in point is the standard deviation
		Point ppg_calculate_rr(Point ppg_intensity_fft_max, Point ppg_amplitude_fft_max);

		// helper function for calculating the magnitude of a complex number
		double cplx_magnitude(kiss_fft_cpx num);
	public:

		PPGAnalysis();
		~PPGAnalysis();

		void run(double *ppg_red_raw, double *ppg_ir_raw, int window_size, int samp_freq);
		
		// the functions which return the results
		double get_hr();
		double get_spo2();
		double get_rr();
		double get_rr_std();
	};
}

#endif