#ifndef MLX90614_H
#define MLX90614_H

#include <stdio.h>
#include <bcm2835.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <stdint.h>

class MLX90614 
{
private:
	const static int ADDR = 0x5A;

	unsigned int AVG;
	
	const unsigned char ambient_register = 0x06;
	const unsigned char object_register_1 = 0x07;
	const unsigned char object_register_2 = 0x08;
	
	unsigned char buf[2];
	
public:	 
	 
	MLX90614(int AVG = 5);
	
	double ambient_temperature();
	
	double object_temperature_1();
	double object_temperature_2();
	double object_temperature();

};

#endif
