#include "mlx90614.h"

MLX90614::MLX90614(int AVG) {
	this->AVG = AVG;

	bcm2835_init(); // comment this line if bcm2835 is being initialized somewhere else
	bcm2835_i2c_begin();
	 //bcm2835_i2c_set_baudrate(25000);
    	// set address
    	bcm2835_i2c_setSlaveAddress(MLX90614::ADDR);
}

double MLX90614::ambient_temperature() {
	double calc = 0, temp = 0;

	for(int i = 0; i < AVG; i++) {
		bcm2835_i2c_begin();
		bcm2835_i2c_write ((const char*)&ambient_register, 1);
		bcm2835_i2c_read_register_rs((char*)&ambient_register, (char*)&buf[0], 2);
		temp = (double) (((buf[1]) << 8) + buf[0]);
		temp = (temp * 0.02)-0.01;
		temp = temp - 273.15;
		calc += temp;
	}

	return(calc/AVG);
}

double MLX90614::object_temperature_1() {
	double calc = 0, temp = 0; 

	for(int i = 0; i < AVG; i++) {
		bcm2835_i2c_begin();
		bcm2835_i2c_write ((const char*) &object_register_1, 1);
		bcm2835_i2c_read_register_rs((char*)&object_register_1, (char*)&buf[0], 2);

		temp = (double) (((buf[1]) << 8) + buf[0]);
		temp = (temp * 0.02)-0.01;
		temp = temp - 273.15;
		calc += temp;
	}

	return(calc/AVG);
}

double MLX90614::object_temperature_2() {
	double calc = 0, temp = 0; 

	for(int i = 0; i < AVG; i++) {
		bcm2835_i2c_begin();
		bcm2835_i2c_write ((const char*) &object_register_2, 1);
		bcm2835_i2c_read_register_rs((char*)&object_register_2, (char*)&buf[0], 2);

		temp = (double) (((buf[1]) << 8) + buf[0]);
		temp = (temp * 0.02)-0.01;
		temp = temp - 273.15;
		calc += temp;
	}

	return(calc/AVG);
}

double MLX90614::object_temperature() {
	//return (this->object_temperature_1() + this->object_temperature_2()) / 2;
	return this->object_temperature_1();
}
