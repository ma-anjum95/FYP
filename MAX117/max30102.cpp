/** \file max30102.cpp ******************************************************
*
* Project: MAXREFDES117#
* Filename: max30102.cpp
* Description: This module is an embedded controller driver for the MAX30102
*
* Revision History:
*\n 1-18-2016 Rev 01.00 GL Initial release.
*\n
*
* --------------------------------------------------------------------
*
* This code follows the following naming conventions:
*
* char              ch_pmod_value
* char (array)      s_pmod_s_string[16]
* float             f_pmod_value
* int32_t           n_pmod_value
* int32_t (array)   an_pmod_value[16]
* int16_t           w_pmod_value
* int16_t (array)   aw_pmod_value[16]
* uint16_t          uw_pmod_value
* uint16_t (array)  auw_pmod_value[16]
* uint8_t           uch_pmod_value
* uint8_t (array)   auch_pmod_buffer[16]
* uint32_t          un_pmod_value
* int32_t *         pn_pmod_value
*
* ------------------------------------------------------------------------- */
/*******************************************************************************
* Copyright (C) 2016 Maxim Integrated Products, Inc., All Rights Reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL MAXIM INTEGRATED BE LIABLE FOR ANY CLAIM, DAMAGES
* OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
* OTHER DEALINGS IN THE SOFTWARE.
*
* Except as contained in this notice, the name of Maxim Integrated
* Products, Inc. shall not be used except as stated in the Maxim Integrated
* Products, Inc. Branding Policy.
*
* The mere transfer of this software does not imply any licenses
* of trade secrets, proprietary technology, copyrights, patents,
* trademarks, maskwork rights, or any other form of intellectual
* property whatsoever. Maxim Integrated Products, Inc. retains all
* ownership rights.
*******************************************************************************
*/

#include "max30102.h"

#include <iostream>
using namespace std;

bool maxim_max30102_write_reg(uint8_t uch_addr, uint8_t uch_data)
/**
* \brief        Write a value to a MAX30102 register
* \par          Details
*               This function writes a value to a MAX30102 register
*
* \param[in]    uch_addr    - register address
* \param[in]    uch_data    - register data
*
* \retval       true on success
*/
{
	const uint8_t buf[2] = {uch_addr, uch_data};
	
	if (!bcm2835_i2c_begin())
		return false;
	bcm2835_i2c_setSlaveAddress(0x57);
	
	if (bcm2835_i2c_write((const char*)buf, 2) != BCM2835_I2C_REASON_OK)
		return false;
	bcm2835_i2c_end();
	
	return true;
}

bool maxim_max30102_read_reg(uint8_t uch_addr, uint8_t *puch_data)
/**
* \brief        Read a MAX30102 register
* \par          Details
*               This function reads a MAX30102 register
*
* \param[in]    uch_addr    - register address
* \param[out]   puch_data    - pointer that stores the register data
*
* \retval       true on success
*/
{
  
	if (!bcm2835_i2c_begin())
		return false;
	bcm2835_i2c_setSlaveAddress(0x57);

	if (bcm2835_i2c_read_register_rs((char*)&uch_addr, (char*)puch_data, 1) != BCM2835_I2C_REASON_OK)
		return false;
	bcm2835_i2c_end();
    
	return true;
}

bool maxim_max30102_init()
/**
* \brief        Initialize the MAX30102
* \par          Details
*               This function initializes the MAX30102
*
* \param        None
*
* \retval       true on success
*/
{
	bcm2835_init();
	bcm2835_i2c_begin();
	bcm2835_i2c_setSlaveAddress(0x57);
	
	// resetting the device
	maxim_max30102_reset();
  
	if(!maxim_max30102_write_reg(REG_INTR_ENABLE_1,0xc0)) // INTR setting
		return false;
	if(!maxim_max30102_write_reg(REG_INTR_ENABLE_2,0x00))
		return false;
	if(!maxim_max30102_write_reg(REG_FIFO_WR_PTR,0x00))  //FIFO_WR_PTR[4:0]
		return false;
	if(!maxim_max30102_write_reg(REG_OVF_COUNTER,0x00))  //OVF_COUNTER[4:0]
		return false;
	if(!maxim_max30102_write_reg(REG_FIFO_RD_PTR,0x00))  //FIFO_RD_PTR[4:0]
		return false;
	if(!maxim_max30102_write_reg(REG_FIFO_CONFIG,0x4f))  //sample avg = 4, fifo rollover=false, fifo almost full = 17
		return false;
	if(!maxim_max30102_write_reg(REG_MODE_CONFIG,0x03))   //0x02 for Red only, 0x03 for SpO2 mode 0x07 multimode LED
		return false;
	if(!maxim_max30102_write_reg(REG_SPO2_CONFIG,0x27))  // SPO2_ADC range = 4096nA, SPO2 sample rate (100 Hz), LED pulseWidth (411uS)
		return false;
  
	if(!maxim_max30102_write_reg(REG_LED1_PA,0x24))   //Choose value for ~ 7mA for LED1
		return false;
	if(!maxim_max30102_write_reg(REG_LED2_PA,0x24))   // Choose value for ~ 7mA for LED2
		return false;
	if(!maxim_max30102_write_reg(REG_PILOT_PA,0x7f))   // Choose value for ~ 25mA for Pilot LED
		return false;
	return true;  
}

bool maxim_max30102_read_fifo(uint32_t *pun_red_led, uint32_t *pun_ir_led)
/**
* \brief        Read a set of samples from the MAX30102 FIFO register
* \par          Details
*               This function reads a set of samples from the MAX30102 FIFO register
*
* \param[out]   *pun_red_led   - pointer that stores the red LED reading data
* \param[out]   *pun_ir_led    - pointer that stores the IR LED reading data
*
* \retval       true on success
*/
{
	uint8_t un_temp[6] = {0};
	uint8_t uch_temp, cr1, cr2;
	uint8_t buf[1] = {REG_FIFO_DATA};
	  
	*pun_ir_led = 0;
	*pun_red_led = 0;
	  
	// reading the status of read and write pointer and checking if they are in correct positions
	maxim_max30102_read_reg(REG_FIFO_WR_PTR , &cr2);
	maxim_max30102_read_reg(REG_FIFO_RD_PTR, &cr1);	
	  
	maxim_max30102_read_reg(REG_INTR_STATUS_1, &uch_temp);
	maxim_max30102_read_reg(REG_INTR_STATUS_2, &uch_temp);
	  
	if (!bcm2835_i2c_begin())
		return false;  
	bcm2835_i2c_setSlaveAddress(0x57);
	
	if ((cr2 - cr1 + 32) % 32 <= 0) // 32 is the size of the array so we mod with that r < 0 ? r + b : r
		return false;
	if (bcm2835_i2c_write_read_rs((char*)buf, 1, (char*)un_temp, 6) != BCM2835_I2C_REASON_OK)
		return false;
	bcm2835_i2c_end();

	*pun_red_led = (un_temp[0] << 16) | (un_temp[1] << 8) | un_temp[2];
	*pun_ir_led = (un_temp[3] << 16) | (un_temp[4] << 8) | un_temp[5];
	
	*pun_red_led &= 0x03FFFF;  //Mask MSB [23:18]
	*pun_ir_led &= 0x03FFFF;  //Mask MSB [23:18]
	 
	return true;
}

bool maxim_max30102_reset()
/**
* \brief        Reset the MAX30102
* \par          Details
*               This function resets the MAX30102
*
* \param        None
*
* \retval       true on success
*/
{
    if(!maxim_max30102_write_reg(REG_MODE_CONFIG,0x40))
        return false;
    else
        return true;    
}

