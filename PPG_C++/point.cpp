#include "point.h"

namespace PPG
{
	Point::Point()
	{
		this->x = 0;
		this->y = 0;
	}

	Point::Point(double x, double y)
	{
		this->x = x;
		this->y = y;
	}

	double Point::get_x(void)
	{
		return this->x;
	}

	double Point::get_y(void)
	{
		return this->y;
	}

	void Point::set(double x, double y)
	{
		this->x = x;
		this->y = y;
	}
}