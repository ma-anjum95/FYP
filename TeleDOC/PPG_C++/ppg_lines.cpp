#include "ppg_lines.h"

using namespace std;

namespace PPG 
{
	PPGLines::PPGLines(double x1, double y1, double x2, double y2):
		p1(x1, y1), p2(x2, y2)
	{
		if (y1 > y2) {
			this->max = y1;
			this->min = y2;
		} else {
			this->max = y2;
			this->min = y1;
		}

		this->amp = this->max - this->min;

		this->period = x2 - x1;

		this->slope = double(y2 - y1) / double(x2 - x1);

		if (this->slope > (double)0.0)
			this->slope_dir = 1;
		else if (this->slope < (double)0.0)
			this->slope_dir = -1;
		else
			this->slope_dir = 0;
	}

	PPGLines::PPGLines(Point p1, Point p2):
		PPGLines(p1.get_x(), p1.get_y(), p2.get_x(), p2.get_y()){}

	PPGLines::PPGLines(const PPGLines &to_copy): 
		PPGLines(to_copy.p1, to_copy.p2) {}

		// getter functions
	double PPGLines::get_x1(void) 
	{
		return this->p1.get_x();
	}

	double PPGLines::get_y1(void) 
	{
		return this->p1.get_y();
	}

	double PPGLines::get_x2(void) 
	{
		return this->p2.get_x();
	}

	double PPGLines::get_y2(void) 
	{
		return this->p2.get_y();
	}

	double PPGLines::get_amp(void)
	{
		return this->amp;
	}

	double PPGLines::get_max(void)
	{
		return this->max;
	}

	double PPGLines::get_min(void)
	{
		return this->min;
	}

	double PPGLines::get_period(void)
	{
		return this->period;
	}
	
	double PPGLines::get_slope(void)
	{
		return this->slope;
	}

	double PPGLines::get_slope_dir(void)
	{
		return this->slope_dir;
	}
}