#ifndef MSREG_H
#define MSREG_H

#include <stddef.h>
#include <stdint.h>

#define MSREG_VERSION "1.0"
#define MSREG_EXIT_OK 0
#define MSREG_EXIT_INPUT 2
#define MSREG_EXIT_GATE 3
#define MSREG_EXIT_INTERNAL 4
#define MSREG_MAX_STAGES 4
#define MSREG_MAX_WARNINGS 8

typedef struct {
    char message[1024];
} MsError;

typedef struct {
    char **lines;
    size_t nlines;
    char *path;
    char **keys;
    char **values;
    size_t nfields;
} MsPrm;

typedef struct {
    char *prm_path;
    char *path;
    MsPrm prm;
    int nx;
    int ny;
    size_t map_size;
    int fd;
    const int16_t *samples;
} MsSlc;

typedef struct {
    int h;
    int w;
    float *v;
} MsImage;

typedef struct {
    double x;
    double y;
} MsXY;

typedef struct {
    double x;
    double y;
    double dr;
    double da;
    double ncc;
    double psr;
    double edge;
    int valid;
    int inlier;
} MsPoint;

typedef struct {
    double c[3];
} MsAxis;

typedef struct {
    MsAxis r;
    MsAxis a;
} MsModel;

typedef struct {
    uint64_t state_hi;
    uint64_t state_lo;
    uint64_t inc_hi;
    uint64_t inc_lo;
    int has_u32;
    uint32_t cached_u32;
} MsRng;

typedef struct {
    int W;
    int R;
    int nxp;
    int nyp;
    double thr;
    int min_inl;
    double min_frac;
    int ovs_a;
    int ovs_r;
} MsSchedule;

typedef struct {
    int ok;
    int looks;
    double ncc;
    double dr;
    double da;
    double psr;
    double second;
    double peak_ratio;
    double time_s;
    double time_read_s;
    char reason[512];
    MsImage texture;
} MsGlobal;

typedef struct {
    int window;
    int search;
    int n_points;
    int n_measured;
    int n_valid;
    int n_inliers;
    double inlier_frac;
    double rms_r;
    double rms_a;
    double sigma;
    double t_in;
    double mad_r;
    double mad_a;
    double robust_sigma_r;
    double robust_sigma_a;
    double rmse_r;
    double rmse_a;
    int metric_n_points;
    double time_s;
    int accepted;
    int has_fit;
    char reason[256];
    double center_dr;
    double center_da;
    MsModel model;
} MsStage;

typedef struct {
    int n_points;
    int n_valid;
    int n_within_0p3;
    double valid_frac;
    double frac_within_0p3_all;
    int has_metrics;
    double median_r;
    double median_a;
    double mad_r;
    double mad_a;
    double robust_sigma_r;
    double robust_sigma_a;
    double rmse_r;
    double rmse_a;
    double max_abs_r;
    double max_abs_a;
    double frac_within_0p3;
    double time_s;
} MsCheck;

typedef struct {
    MsSlc master;
    MsSlc slave;
    MsGlobal global;
    MsStage stages[MSREG_MAX_STAGES];
    int nstages;
    char warnings[MSREG_MAX_WARNINGS][512];
    int nwarnings;
    int status;
    char fail_reason[512];
    int has_model;
    MsModel model;
    char final_stage[16];
    double center_dr;
    double center_da;
    char prm_values[8][64];
    MsPoint *final_points;
    size_t nfinal_points;
    int final_window;
    int coverage_cells;
    MsCheck check;
    int has_check;
    int gate_passed;
    char gate_reasons[8][512];
    int ngate_reasons;
    double time_total_s;
    double check_time_s;
    char output_prm[1024];
} MsEstimate;

typedef struct {
    MsCheck check;
    int passed;
    char reasons[8][512];
    int nreasons;
    double time_s;
} MsVerify;

typedef struct {
    int global_target;
    int seed;
    int threads;
    int write_on_fail;
    int quiet;
    int prep;
    int has_ovs;
    int ovs_a;
    int ovs_r;
    const char *scales;
    MsSchedule schedule[MSREG_MAX_STAGES];
    int nschedule;
    int grid_x;
    int grid_y;
    int grid_is_set;
    const char *master;
    const char *slave;
    const char *slave_slc;
    const char *output;
    const char *report;
    const char *log;
    const char *points;
} MsOptions;

int ms_prm_read(const char *path, MsPrm *p, MsError *err);
void ms_prm_free(MsPrm *p);
const char *ms_prm_get(const MsPrm *p, const char *key);
int ms_prm_write(const MsPrm *src, const char *dst, const MsModel *model,
                 const char *slc_file_override, MsError *err);
int ms_slc_open(MsSlc *s, const char *prm_path, const char *slc_override,
                MsError *err);
void ms_slc_close(MsSlc *s);
void ms_slc_read(const MsSlc *s, int y0, int x0, int h, int w, float *out_ri);
int ms_multilook(const MsSlc *s, int looks, MsImage *out, MsError *err);
void ms_box_mean(const float *src, int h, int w, int k, float *dst);
double ms_median(const double *v, size_t n);

void ms_rng_seed(MsRng *rng, uint64_t seed);
uint64_t ms_rng_next64(MsRng *rng);
uint32_t ms_rng_next32(MsRng *rng);
uint32_t ms_rng_bounded32(MsRng *rng, uint32_t upper_exclusive);
uint64_t ms_rng_bounded64(MsRng *rng, uint64_t upper_exclusive);
int ms_rng_dump_test(uint64_t seed, unsigned long count);
int ms_rng_pattern_test(uint64_t seed, uint32_t first_upper,
                        uint32_t second_upper);
int ms_fit_affine(const MsPoint *points, size_t n, double threshold,
                  MsRng *rng, MsModel *model, unsigned char *inliers,
                  MsStage *stats);
void ms_predict(const MsModel *model, double x, double y, double *dr, double *da);
void ms_prm_values(const MsModel *model, char out[8][64]);
double ms_mad(const double *v, size_t n);
void ms_residual_stats(const double *v, size_t n, double *mad, double *rmse);

int ms_global_stage(const MsSlc *master, const MsSlc *slave, int target,
                    int threads, MsGlobal *out, MsError *err);
void ms_global_free(MsGlobal *g);
int ms_make_lattice(const MsSlc *master, const MsSlc *slave,
                    const MsModel *model, int W, int R, int nxp, int nyp,
                    const MsImage *texture, int looks, double phase_x,
                    double phase_y, MsXY **out, size_t *nout);
int ms_measure_points(const MsSlc *master, const MsSlc *slave,
                      const MsModel *model, const MsXY *xy, size_t nxy,
                      int W, int R, int prep, int ovs_a, int ovs_r,
                      int threads, MsPoint **out, size_t *nout, MsError *err);
int ms_match_patch(const float *tm, const float *ss, int h, int w, int H, int W,
                   int threads, double *dy, double *dx, double *ncc,
                   double *psr, int *edge, MsError *err);
int ms_estimate(const MsOptions *opt, MsEstimate *out, MsError *err);
int ms_check_model(const MsSlc *master, const MsSlc *slave,
                   const MsModel *model, int grid_x, int grid_y, int prep,
                   int ovs_a, int ovs_r, int threads, int use_texture,
                   MsCheck *out, MsError *err);
int ms_verify_registered(const MsSlc *master, const MsSlc *slave,
                         const MsOptions *opt, MsVerify *out, MsError *err);
int ms_write_estimate_report(const char *path, const MsEstimate *e,
                             MsError *err);
int ms_write_verify_report(const char *path, const MsSlc *master,
                           const MsSlc *slave, const MsVerify *v,
                           MsError *err);
int ms_write_log(const char *path, const MsEstimate *e, const char *master,
                 const char *slave, MsError *err);
int ms_write_points(const char *path, const MsEstimate *e, MsError *err);
void ms_log(const MsOptions *opt, const char *fmt, ...);
int ms_estimate_cli(const MsOptions *opt);
int ms_verify_cli(const MsOptions *opt);

#endif
