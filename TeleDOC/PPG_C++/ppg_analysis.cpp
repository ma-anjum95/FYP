#include "ppg_analysis.h"

using namespace std;

namespace PPG
{

	PPGAnalysis::PPGAnalysis()
	{
		this->ppg_red_hpf = NULL;
		this->ppg_ir_hpf = NULL;

		this->ppg_red_ac = NULL;
		this->ppg_ir_ac = NULL;
		
		this->ppg_red_dc = 0;
		this->ppg_ir_dc = 0;

		this->ppg_red_lines = new vector<PPGLines>;
		this->ppg_ir_lines = new vector<PPGLines>;
		this->ppg_red_lines_processed = new vector<PPGLines>;
		this->ppg_ir_lines_processed = new vector<PPGLines>;
		this->ppg_red_intensity_waveform = new vector<Point>;
		this->ppg_ir_intensity_waveform = new vector<Point>;
		this->ppg_red_amplitude_waveform = new vector<Point>;
		this->ppg_ir_amplitude_waveform = new vector<Point>;
        this->ppg_red_period_waveform = new vector<Point>;
        this->ppg_ir_period_waveform = new vector<Point>;
		this->ppg_red_intensity_waveform_interp = new vector<double>;
		this->ppg_ir_intensity_waveform_interp = new vector<double>;
		this->ppg_red_amplitude_waveform_interp = new vector<double>;
		this->ppg_ir_amplitude_waveform_interp = new vector<double>;
        this->ppg_red_period_waveform_interp = new vector<double>;
        this->ppg_ir_period_waveform_interp = new vector<double>;
		this->ppg_red_intensity_waveform_fft = new vector<Point>;
		this->ppg_ir_intensity_waveform_fft = new vector<Point>;
		this->ppg_red_amplitude_waveform_fft = new vector<Point>;
		this->ppg_ir_amplitude_waveform_fft = new vector<Point>;
        this->ppg_red_period_waveform_fft = new vector<Point>;
        this->ppg_ir_period_waveform_fft = new vector<Point>;

		this->hr = 0;
		this->spo2 = 0;
		this->rr = 0;
	}

	PPGAnalysis::~PPGAnalysis()
	{
		delete this->ppg_red_lines;
		delete this->ppg_ir_lines;
		delete this->ppg_red_lines_processed;
		delete this->ppg_ir_lines_processed;
		delete this->ppg_red_intensity_waveform;
		delete this->ppg_ir_intensity_waveform;
		delete this->ppg_red_amplitude_waveform;
		delete this->ppg_ir_amplitude_waveform;
        delete this->ppg_red_period_waveform;
        delete this->ppg_ir_period_waveform;
		delete this->ppg_red_intensity_waveform_interp;
		delete this->ppg_ir_intensity_waveform_interp;
		delete this->ppg_red_amplitude_waveform_interp;
		delete this->ppg_ir_amplitude_waveform_interp;
        delete this->ppg_red_period_waveform_interp;
        delete this->ppg_ir_period_waveform_interp;
        delete this->ppg_red_intensity_waveform_fft;
        delete this->ppg_ir_intensity_waveform_fft;
        delete this->ppg_red_amplitude_waveform_fft;
        delete this->ppg_ir_amplitude_waveform_fft;
        delete this->ppg_red_period_waveform_fft;
        delete this->ppg_ir_period_waveform_fft;

	}

	void PPGAnalysis::clean_vectors(void)
	{
		this->ppg_red_lines->clear();
		this->ppg_ir_lines->clear();
		this->ppg_red_lines_processed->clear();
		this->ppg_ir_lines_processed->clear();
		this->ppg_red_intensity_waveform->clear();
		this->ppg_ir_intensity_waveform->clear();
		this->ppg_red_amplitude_waveform->clear();
		this->ppg_ir_amplitude_waveform->clear();
		this->ppg_red_intensity_waveform_interp->clear();
		this->ppg_ir_intensity_waveform_interp->clear();
		this->ppg_red_amplitude_waveform_interp->clear();
		this->ppg_ir_amplitude_waveform_interp->clear();
		this->ppg_red_intensity_waveform_fft->clear();
		this->ppg_ir_intensity_waveform_fft->clear();
		this->ppg_red_amplitude_waveform_fft->clear();
		this->ppg_ir_amplitude_waveform_fft->clear();

	}

	void PPGAnalysis::ppg_hpf(double *ppg_raw, double *&ppg_hpf)
	{
		ppg_hpf = ppg_raw;
	}

	void PPGAnalysis::ppg_dcf(double *ppg_hpf, double *&ppg_ac, double &ppg_dc)
	{
		double sum = 0;
		ppg_ac = new double[this->window_size];

		for (int i = 0; i < this->window_size; i++) {
			sum += ppg_hpf[i];
		}

		ppg_dc = (double)sum / (double)this->window_size;

		for (int i = 0; i < this->window_size; i++) {
			ppg_ac[i] = ppg_hpf[i] - ppg_dc;
		}
	}

	void PPGAnalysis::ppg_line_segmentation(vector<PPGLines> *to_return, double *ppg_ac, int m)
	{
		int seg = 1;
		int z = 0;
		int seg_in_line = 1;

		// adding the first element 
		to_return->push_back(PPGLines(seg, ppg_ac[seg - 1], seg + m, ppg_ac[seg + m - 1]));

		z += 1;
		seg += 1;

		while ((seg + 1) * m < this->window_size) {

			// this insertion is at index z in case the z doesnot exist we push back else we replace
			if (z < to_return->size())
				to_return[0][z] = PPGLines(seg*m, ppg_ac[seg*m - 1], (seg + 1)*m, ppg_ac[(seg + 1)*m - 1]);
			else
				to_return->push_back(PPGLines(seg*m, ppg_ac[seg*m - 1], (seg + 1)*m, ppg_ac[(seg + 1)*m - 1]));


			if (to_return[0][z].get_slope_dir() * to_return[0][z - 1].get_slope_dir() > 0) {
				to_return[0][z - 1] = PPGLines((seg - seg_in_line)*m, ppg_ac[(seg - seg_in_line)*m - 1], (seg + 1)*m, ppg_ac[(seg + 1)*m - 1]);
				seg += 1;
				seg_in_line += 1;
			} else {
				z += 1;
				seg += 1;
				seg_in_line = 1;
			}
		}
	}

	void PPGAnalysis::ppg_artifact_removal(vector<PPGLines> *to_return, vector<PPGLines> *ppg_lines)
	{
		// if the number of segments are less than or equal to 2 just return as there is no use in running artifact removal
		if (ppg_lines->size() <= 2)
			return;

		// The artifacts in the ppg lines will be of the following types :
		// 1 - Lines before or after a horizontal line will be artifacts
		// 2 - Lines having amplitudes exceeding the limits ThAhigh and ThAlow.
		// 3 - Beats having period less than 240ms.
		// 4 - Beats with less than 50 % of previous valid interbeat - interval

		// Calculating the value of ThAlow and ThAhigh which are amplitude limits.
		int ThAhigh = 1500;
		int ThAlow = 250;

		// Beat thresholds within with the human heart beat remains
		// Heart beat must remain between 230ms to 2400ms.
		double Period_Low_Limit = .23 * this->samp_freq / 2; // 230ms
		double Period_High_Limit = 2.4 * this->samp_freq / 2; // 2400ms

		for (int i = 1, n = ppg_lines->size() - 1; i < n; i++) {
			if ((ppg_lines[0][i].get_slope_dir() * ppg_lines[0][i - 1].get_slope_dir() >= 0) || (ppg_lines[0][i].get_slope_dir() * ppg_lines[0][i + 1].get_slope_dir() >= 0)) {
				continue;
			}
			else if ((ppg_lines[0][i].get_amp() > ThAhigh) || (ppg_lines[0][i].get_amp() < ThAlow)) {
				continue;
			} 
			else if ((ppg_lines[0][i].get_period() < Period_Low_Limit)
				|| (ppg_lines[0][i].get_period() > Period_High_Limit)) {
				continue;
			}
			else {
				to_return->push_back(PPGLines(ppg_lines[0][i]));
			}
		}
	}
	
	void PPGAnalysis::ppg_intensity_waveform(vector<Point> *to_return, vector<PPGLines> *ppg_lines_processed)
	{	
		// intensity of the previous pulse
		double x_prev = -1;
		double x, y;

		for (int i = 0, n = ppg_lines_processed->size(); i < n; i++) {
			if (ppg_lines_processed[0][i].get_slope_dir() == 1) {
				// if slope is greater than 0 we take x2 and y2
				x = ppg_lines_processed[0][i].get_x2();
				y = ppg_lines_processed[0][i].get_y2();

				if (x_prev != x) {
					to_return->push_back(Point(x, y));
					x_prev = x;
				}
			}
			else {
				// if slope is less than 0 we take x1 and y1
				x = ppg_lines_processed[0][i].get_x1();
				y = ppg_lines_processed[0][i].get_y1();
			
				if (x_prev != x) {
					to_return->push_back(Point(x, y));
					x_prev = x;
				}
			}
		}
	}
	
	void PPGAnalysis::ppg_amplitude_waveform(vector<Point> *to_return, vector<PPGLines> *ppg_lines_processed)
	{
		
		for (int i = 0, n = ppg_lines_processed->size(); i < n; i++) {
			if (ppg_lines_processed[0][i].get_slope_dir() == 1)
				to_return->push_back(Point(ppg_lines_processed[0][i].get_x1(), ppg_lines_processed[0][i].get_amp()));
		}
	}

    void PPGAnalysis::ppg_period_waveform(vector<Point> *to_return, vector<PPGLines> *ppg_lines_processed)
    {
        for(int i = 0, n = ppg_lines_processed->size(); i < n - 1; i++) {
            if (ppg_lines_processed[0][i].get_slope_dir() == 1) {
                if (ppg_lines_processed[0][i + 1].get_slope_dir() == -1 && (ppg_lines_processed[0][i + 1].get_x1() - ppg_lines_processed[0][i].get_x2()) < 25)
                    to_return->push_back(Point(ppg_lines_processed[0][i].get_x1(), ppg_lines_processed[0][i].get_period() + ppg_lines_processed[0][i + 1].get_period()));
            }
        }
    }

	void PPGAnalysis::linear_interpolation_with_freq(vector<double> *to_return, vector<Point> *data, int old_freq, int new_freq)
	{
		// This function resamples the data using linear interpolation from old_freq to the new_freq
		
		double freq_ratio = (double)new_freq / old_freq;
		int len = data->size();
		int *new_index = new int[len];
		double *data_slope = new double[len];
		int x_interp = 0; // the x axis value for the interpolated data
		double y_interp;
		int z = 0; // the index of the data set and new_index 

		double x0, y0, x1, y1;

		// creating a new index array for the new_freq
		int prev_index = -1;
		for (int i = 0, n = len; i < n; i++) {
			new_index[i] = round(freq_ratio * data[0][i].get_x());
			// if the previous new index and the current new index are equal add one to the the new index
			if (prev_index == new_index[i])
				new_index[i] += 1;
			prev_index = new_index[i];
		}

		// we need to calculate slope of the data
		for (int i = 0, n = len - 1; i < n; i++) {
			data_slope[i] = (double)(data[0][i + 1].get_y() - data[0][i].get_y()) / (data[0][i + 1].get_x() - data[0][i].get_x());
		}
		
		while (z < len - 1) {
			x0 = new_index[z];
			y0 = data[0][z].get_y();

			x1 = new_index[z + 1];
			y1 = data[0][z + 1].get_y();

			// if the current interpolated point is greater than x1 we increment i and continue
			if (x_interp > x1) {
				z += 1;
				continue;
			}
			else {
				y_interp = ((x_interp - x0) / (x1 - x0)) * (y1 - y0) + y0;

				to_return->push_back(y_interp);

				x_interp += 1;
			}
		}
		
	}

	double PPGAnalysis::ppg_heart_rate(vector<PPGLines> *ppg_lines_processed)
	{
		double sum = 0;
		double j = 0;

		for (int i = 0; i < ppg_lines_processed->size() - 1; i++) {
			if (ppg_lines_processed[0][i].get_slope_dir() == 1) {
				if (ppg_lines_processed[0][i + 1].get_slope_dir() == -1) {
					sum += ppg_lines_processed[0][i].get_period() + ppg_lines_processed[0][i + 1].get_period();
					j++;
				}
			}
		}

		if (j == 0)
			return 0;
		else
			return (double)(j / sum) * 60.0 * this->samp_freq ;
	}

	double PPGAnalysis::ppg_spo2(double *ppg_red_ac, double ppg_red_dc, double *ppg_ir_ac, double ppg_ir_dc)
	{
		double spo2;

		double ppg_ac_rms_red = this->calculate_rms(ppg_red_ac, this->window_size);
		double ppg_ac_rms_ir = this->calculate_rms(ppg_ir_ac, this->window_size);

		double R = (ppg_ac_rms_red / ppg_red_dc) / (ppg_ac_rms_ir / ppg_ir_dc);

		spo2 = -45.060*R*R/10000 + 30.354 *R/100 + 94.845;

		return spo2;
	}

	double PPGAnalysis::calculate_rms(double *ac_signal, int length)
	{
		// calculates rms from the given data
		// rms = sqrt(1 / n * (sum data(i) ^ 2))

		double squared_sum = 0;

		for (int i = 0; i < length; i++) {
			squared_sum += ac_signal[i] * ac_signal[i]; // adds the square of the point
		}

		double rms = sqrt(((double)1.0/(double)length) * (double)squared_sum);

		return rms;
	}

	void PPGAnalysis::ppg_waveform_fft(vector<Point> *to_return, vector<double> *waveform, int freq)
	{
		int n = waveform->size();

		// initializing the fft configuration
		kiss_fft_cfg fft_cfg = kiss_fft_alloc(n, 0, NULL, NULL);
	
		//void kiss_fft(kiss_fft_cfg cfg, const kiss_fft_cpx *fin, kiss_fft_cpx *fout);
		// initializing the input and output buffers
		kiss_fft_cpx *fin = new kiss_fft_cpx[n];
		kiss_fft_cpx *fout = new kiss_fft_cpx[n];

		// preparing the input values
		// setting the fin to the waveform values
		for (int i = 0; i < n; i++) {
			fin[i].r = waveform[0][i];
			fin[i].i = 0;
		}

		// running the fft algorithm
		kiss_fft(fft_cfg, fin, fout);

		/* The following lines of code implement these threee lines
			P2 = abs(Y/L);
			P1 = P2(1:L/2+1);
			P1(2:end-1) = 2*P1(2:end-1);
		*/
		double *P2 = new double[n];
		for (int i = 0; i < n; i++) {
			P2[i] = (double) this->cplx_magnitude(fout[i]) / (double)n;
		}

		double *P1 = new double[n/2 + 1];
		for (int i = 0, tmp = n/2 + 1; i < tmp; i++) {
			P1[i] = P2[i];
		}

		for (int i = 1, tmp = n/2 + 1; i < tmp - 1; i++) {
			P1[i] = 2 * P1[i];
		}

		// f = freq*(0:(L/2))/L;
		for (int i = 0, tmp = n / 2 + 1; i < tmp; i++) {
			double f = (freq * i) / (double)n;
			to_return->push_back(Point(f, P1[i]));
		}

		// freeing the memory
		free(fft_cfg);
		delete fin;
		delete fout;
	}

	Point PPGAnalysis::find_max_in_range(vector<Point> *arr, double low, double high)
	{
		int i = 0; // index of the array
		Point max; // the maximum point in the array
		int tmp;
		int n = arr->size();

		// increment i till low is greater than the required value
		while (low > arr[0][i].get_x() && i < n)
			i++;

		tmp = i;

		while (arr[0][i].get_x() < high && i < n) {
			if (arr[0][i].get_y() > arr[0][tmp].get_y())
				tmp = i;
			i++;
		}

		max = arr[0][tmp];
		
		return max;
	}

    Point PPGAnalysis::ppg_calculate_rr(Point ppg_intensity_fft_max, Point ppg_amplitude_fft_max, Point ppg_period_fft_max)
	{

		double rr_intensity = ppg_intensity_fft_max.get_x() * 60;
		double rr_amplitude = ppg_amplitude_fft_max.get_x() * 60;
        double rr_period = ppg_period_fft_max.get_x() * 60;

		// respiratory rate is the mean of all the rr estimations
        double rr = (rr_intensity + rr_amplitude + rr_period) / 3;

		// to calculate quality we calculate sd of individual rr and see if it is less the 4
		double rr_dev_intensity = rr_intensity - rr;
		double rr_dev_amplitude = rr_amplitude - rr;
        double rr_dev_period = rr_period -rr;

		// the standard deviation
        double sd = sqrt((rr_dev_amplitude * rr_dev_amplitude + rr_dev_intensity * rr_dev_intensity + rr_dev_period * rr_dev_period) / 3);

		return Point(rr, sd);
	}

	double PPGAnalysis::cplx_magnitude(kiss_fft_cpx num)
	{
		double real_sqr = (double)(num.r * num.r);
		double imag_sqr = (double)(num.i * num.i);
		double sum = real_sqr + imag_sqr;

		return pow(sum, 0.5);
	}

	double PPGAnalysis::get_hr()
	{
		return this->hr;
	}

	double PPGAnalysis::get_spo2()
	{
		return this->spo2;
	}

	double PPGAnalysis::get_rr()
	{
		return this->rr;
	}

	double PPGAnalysis::get_rr_std()
	{
		return this->rr_std;
	}

	void PPGAnalysis::run(double *ppg_red_raw, double *ppg_ir_raw, int window_size, int samp_freq) 
	{		
		// cleans the contents of vector before running the complete algorithm
		this->clean_vectors();

		this->window_size = window_size;
		this->samp_freq = samp_freq;

		// passing the data through a high pass filter 
		this->ppg_hpf(ppg_red_raw, this->ppg_red_hpf);
		this->ppg_hpf(ppg_ir_raw, this->ppg_ir_hpf);

		// passing the data through a low pass filter
		this->ppg_dcf(this->ppg_red_hpf, this->ppg_red_ac, this->ppg_red_dc);
		this->ppg_dcf(this->ppg_ir_hpf, this->ppg_ir_ac, this->ppg_ir_dc);

		// converting the ppg_ac into segmented lines
		// this step also does classification assigning max, min, slope and period to lines
		this->ppg_line_segmentation(this->ppg_red_lines, this->ppg_red_ac);
		this->ppg_line_segmentation(this->ppg_ir_lines, this->ppg_ir_ac);

		// removing the artifacts from the newly converted lines
		//this->ppg_artifact_removal(this->ppg_red_lines_processed, this->ppg_red_lines);
		this->ppg_artifact_removal(this->ppg_ir_lines_processed, this->ppg_ir_lines);

		// if the number of segments after artifact removal are less than 3 we just return 0
		if (this->ppg_ir_lines_processed->size() > 3) {
			// extracting the intensity waveform from the processed lines
			//this->ppg_intensity_waveform(this->ppg_red_intensity_waveform, this->ppg_red_lines_processed);
			this->ppg_intensity_waveform(this->ppg_ir_intensity_waveform, this->ppg_ir_lines_processed);

			// extracting the amplitude waveform from the processed lines
			//this->ppg_amplitude_waveform(this->ppg_red_amplitude_waveform, this->ppg_red_lines_processed);
			this->ppg_amplitude_waveform(this->ppg_ir_amplitude_waveform, this->ppg_ir_lines_processed);

            //extracting the period waveform from the processed lines
            // this->ppg_period_waveform(this->ppg_red_period_waveform, this->ppg_red_lines_processed);
            this->ppg_period_waveform(this->ppg_ir_period_waveform, this->ppg_ir_lines_processed);

			// interpolating the waveform data from the processed lines
			//this->linear_interpolation_with_freq(this->ppg_red_intensity_waveform_interp, this->ppg_red_intensity_waveform, this->samp_freq, 4);
			this->linear_interpolation_with_freq(this->ppg_ir_intensity_waveform_interp, this->ppg_ir_intensity_waveform, this->samp_freq, 4);

			//this->linear_interpolation_with_freq(this->ppg_red_amplitude_waveform_interp, this->ppg_red_amplitude_waveform, this->samp_freq, 4);
			this->linear_interpolation_with_freq(this->ppg_ir_amplitude_waveform_interp, this->ppg_ir_amplitude_waveform, this->samp_freq, 4);

            // this->linear_interpolation_with_freq(this->ppg_red_period_waveform_interp, this->ppg_red_period_waveform, this->samp_freq, 4);
            this->linear_interpolation_with_freq(this->ppg_ir_period_waveform_interp, this->ppg_ir_period_waveform, this->samp_freq, 4);

			// fft of the interpolated data
			//this->ppg_waveform_fft(this->ppg_red_intensity_waveform_fft, this->ppg_red_intensity_waveform_interp, 4);
			this->ppg_waveform_fft(this->ppg_ir_intensity_waveform_fft, this->ppg_ir_intensity_waveform_interp, 4);

			//this->ppg_waveform_fft(this->ppg_red_amplitude_waveform_fft, this->ppg_red_amplitude_waveform_interp, 4);
			this->ppg_waveform_fft(this->ppg_ir_amplitude_waveform_fft, this->ppg_ir_amplitude_waveform_interp, 4);

            // this->ppg_waveform(this->ppg_red_period_waveform_fft, this->ppg_red_period_waveform_interp, 4);
            this->ppg_waveform_fft(this->ppg_ir_period_waveform_fft, this->ppg_ir_period_waveform_interp, 4);

			// finding maximum value in fft in range of respiratory freq
			//this->ppg_red_intensity_fft_max = this->find_max_in_range(this->ppg_red_intensity_waveform_fft, 0.067, 1.08);
			this->ppg_ir_intensity_fft_max = this->find_max_in_range(this->ppg_ir_intensity_waveform_fft, 0.067, 1.08);

			//this->ppg_red_amplitude_fft_max = this->find_max_in_range(this->ppg_red_amplitude_waveform_fft, 0.067, 1.08);
			this->ppg_ir_amplitude_fft_max = this->find_max_in_range(this->ppg_ir_amplitude_waveform_fft, 0.067, 1.08);

            //this->ppg_red_period_fft_max = this->find_max_in_range(this->ppg_red_period_waveform_fft, 0.067, 1.08);
            this->ppg_ir_period_fft_max = this->find_max_in_range(this->ppg_ir_period_waveform_fft, 0.067, 1.08);

			// extracting the heart rate from the intensity waveform
			this->hr = this->ppg_heart_rate(this->ppg_ir_lines_processed);

			// extracting spO2 from the singals
			this->spo2 = this->ppg_spo2(this->ppg_red_ac, this->ppg_red_dc, this->ppg_ir_ac, this->ppg_ir_dc);

			// extracting rr from the signals
            this->rr = this->ppg_calculate_rr(this->ppg_ir_intensity_fft_max, this->ppg_ir_amplitude_fft_max, this->ppg_ir_period_fft_max).get_x();

			// setting the standard deviations of different heart rate variations
            this->rr_std = this->ppg_calculate_rr(this->ppg_ir_intensity_fft_max, this->ppg_ir_amplitude_fft_max, this->ppg_ir_period_fft_max).get_y();
		}
		else {
			// extracting the heart rate from the intensity waveform
			this->hr = 0;

			// extracting spO2 from the singals
			this->spo2 = 0;

			// extracting rr from the signals
			this->rr = 0;

			// setting the standard deviations of different heart rate variations
			this->rr_std = 0;
		}
	}
}
