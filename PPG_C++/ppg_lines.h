/*
*	PPG lines class deals with individual pulses that appear in our algorithm
*
*/
#ifndef PPG_LINES_H
#define PPG_LINES_H

#include <iostream>

#include "point.h"

namespace PPG 
{
	class PPGLines 
	{
		Point p1, p2;
		double amp, max, min;
		double period, slope;
		int slope_dir; // 1 for pos, 0 for 0, -1 for neg;

	public:
		PPGLines(double x1, double y1, double x2, double y2);
		PPGLines(Point p1, Point p2);
		PPGLines(const PPGLines &to_copy);

		// getter functions
		double get_x1(void);
		double get_y1(void);
		double get_x2(void);
		double get_y2(void);

		double get_amp(void);
		double get_max(void);
		double get_min(void);

		double get_period(void);
		double get_slope(void);
		double get_slope_dir(void);
	};
}
#endif