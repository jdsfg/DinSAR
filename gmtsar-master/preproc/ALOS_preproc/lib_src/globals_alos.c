/*
 * globals_alos.c
 * 
 * Global variable definitions for ALOS_preproc module.
 * These variables are declared as 'extern' in image_sio.h
 * 
 * Created to fix "multiple definition" linker errors.
 */

/* controls minimal level of output */
int verbose = 0;

/* more output */
int debug = 0;

/* more output */
int roi = 0;

/* whether to swap bytes */
int swap = 0;

/* quad polarization data */
int quad_pol = 0;

/* AUIG: ALOS_format = 0, ERSDAC: ALOS_format = 1, ALOS2: ALOS_format = 2 */
int ALOS_format = 0;

/* whether to set the slope */
int force_slope = 0;

/* whether to calculate doppler */
int dopp = 0;

/* reduce output */
int quiet_flag = 0;

/* SAR_mode: 0 => high-res, 1 => wide obs, 2 => polarimetry */
int SAR_mode = 0;

/* offset needed for ALOS-2 prefix size */
int prefix_off = 0;

/* value to set chirp_slope to */
double forced_slope = 0.0;

/* time bias for clock bias */
double tbias = 0.0;

/* range bias for near range corr */
double rbias = 0.0;

/* factor to convert float to int slc */
double slc_fact = 1.0;
