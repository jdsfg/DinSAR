/* Global variable definitions for GMTSAR */
#include "gmtsar.h"
/* verbose 和 debug 现在在 sio_struct.c 中定义 */
/* 定义其他全局变量并初始化 */
int swap = 0;        /* whether to swap bytes 		*/
int quad_pol = 0;    /* quad polarization data 		*/
int force_slope = 0; /* whether to force the slope 		*/
int dopp = 0;        /* whether to calculate doppler 	*/
int roi_flag = 0;    /* whether to write roi.in 		*/
int sio_flag = 0;    /* whether to write PRM file 		*/
int nodata = 0;
int quiet_flag = 0;
int SAR_mode = 0;    /* 0 => high-res                        */
double forced_slope = 0.0; /* value to set chirp_slope to		*/
