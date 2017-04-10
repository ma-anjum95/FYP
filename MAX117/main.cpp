#include <stdio.h>
//#include <wiringPi.h>
//#include <wiringPiI2C.h>
#include "max30102.h"
#include <time.h>
#include <unistd.h>

// -I/usr/local/include -L/usr/local/lib -lwiringPi

int main(void) 
{
	printf("hello");
	
	// initialze wiringPi to access the GPIO pins on raspberry pi 
	//wiringPiSetupGpio();
	
	
	if (maxim_max30102_init())
		printf("success");
	else
		printf("failure");
	

	uint32_t tmp1, tmp2;
	sleep(2);
	printf("VALUES: \n");
	if (maxim_max30102_read_fifo(&tmp1, &tmp2))
		printf("SUCCESS");
		
	printf("\nred: %08d\n ir: %08d", tmp1, tmp2);

	return true;
	
}
