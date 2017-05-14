#ifndef POINT_H
#define POINT_H

#include <iostream>

namespace PPG
{
	class Point
	{
		double x;
		double y;

	public:
		Point(double x, double y);
		Point();

		double get_x(void);
		double get_y(void);

		void set(double x, double y);
	};
}

#endif