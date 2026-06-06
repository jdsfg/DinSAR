/*
 * globals_alos.c
 * Explicit definition of global variables declared extern in image_sio.h
 */
#include "image_sio.h"

int verbose = 0;
int debug = 0;
int roi = 0;
int swap = 0;
int quad_pol = 0;
int ALOS_format = 0;
int force_slope = 0;
int dopp = 0;
int quiet_flag = 0;
int SAR_mode = 0;
int prefix_off = 0;

double forced_slope = 0.0;
double tbias = 0.0;
double rbias = 0.0;
double slc_fact = 1.0;
