#include "msreg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint32_t hashmix(uint32_t value, uint32_t *hash_const)
{
    value ^= *hash_const;
    *hash_const *= UINT32_C(0x931e8875);
    value *= *hash_const;
    value ^= value >> 16;
    return value;
}

static uint32_t mix32(uint32_t x, uint32_t y)
{
    uint32_t value = UINT32_C(0xca01f9dd) * x - UINT32_C(0x4973f715) * y;
    value ^= value >> 16;
    return value;
}

void ms_rng_seed(MsRng *rng, uint64_t seed)
{
    uint32_t entropy[2];
    size_t n_entropy = 1;
    uint32_t pool[4];
    uint32_t hc = UINT32_C(0x43b0d7e5);
    uint32_t hash_const;
    uint32_t state_words[8];
    uint64_t seed_words[4];
    size_t i, j;
    __uint128_t state;
    __uint128_t inc;
    const __uint128_t mult =
        (((__uint128_t)UINT64_C(2549297995355413924)) << 64) |
        UINT64_C(4865540595714422341);

    entropy[0] = (uint32_t)seed;
    if (seed > UINT32_MAX) {
        entropy[1] = (uint32_t)(seed >> 32);
        n_entropy = 2;
    }
    for (i = 0; i < 4; i++) {
        pool[i] = hashmix(i < n_entropy ? entropy[i] : 0, &hc);
    }
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            if (i != j) {
                pool[j] = mix32(pool[j], hashmix(pool[i], &hc));
            }
        }
    }
    hash_const = UINT32_C(0x8b51f9dd);
    for (i = 0; i < 8; i++) {
        uint32_t v = pool[i % 4] ^ hash_const;
        hash_const *= UINT32_C(0x58f38ded);
        v *= hash_const;
        v ^= v >> 16;
        state_words[i] = v;
    }
    for (i = 0; i < 4; i++) {
        seed_words[i] = (uint64_t)state_words[2 * i] |
                        ((uint64_t)state_words[2 * i + 1] << 32);
    }
    state = 0;
    inc = (((__uint128_t)seed_words[2]) << 64) | seed_words[3];
    inc = (inc << 1) | 1;
    state = state * mult + inc;
    state += (((__uint128_t)seed_words[0]) << 64) | seed_words[1];
    state = state * mult + inc;
    rng->state_hi = (uint64_t)(state >> 64);
    rng->state_lo = (uint64_t)state;
    rng->inc_hi = (uint64_t)(inc >> 64);
    rng->inc_lo = (uint64_t)inc;
    rng->has_u32 = 0;
    rng->cached_u32 = 0;
}

uint64_t ms_rng_next64(MsRng *rng)
{
    const __uint128_t mult =
        (((__uint128_t)UINT64_C(2549297995355413924)) << 64) |
        UINT64_C(4865540595714422341);
    __uint128_t state =
        (((__uint128_t)rng->state_hi) << 64) | rng->state_lo;
    __uint128_t inc = (((__uint128_t)rng->inc_hi) << 64) | rng->inc_lo;
    uint64_t hi;
    uint64_t lo;
    uint64_t x;
    unsigned int rot;
    state = state * mult + inc;
    rng->state_hi = (uint64_t)(state >> 64);
    rng->state_lo = (uint64_t)state;
    hi = rng->state_hi;
    lo = rng->state_lo;
    x = hi ^ lo;
    rot = (unsigned int)(hi >> 58);
    return (x >> rot) | (x << ((-rot) & 63U));
}

uint32_t ms_rng_next32(MsRng *rng)
{
    uint64_t next;
    if (rng->has_u32) {
        rng->has_u32 = 0;
        return rng->cached_u32;
    }
    next = ms_rng_next64(rng);
    rng->has_u32 = 1;
    rng->cached_u32 = (uint32_t)(next >> 32);
    return (uint32_t)next;
}

uint32_t ms_rng_bounded32(MsRng *rng, uint32_t upper_exclusive)
{
    uint32_t x;
    uint32_t threshold;
    uint64_t m;
    if (upper_exclusive <= 1) {
        return 0;
    }
    threshold = (uint32_t)(-upper_exclusive) % upper_exclusive;
    for (;;) {
        x = ms_rng_next32(rng);
        m = (uint64_t)x * upper_exclusive;
        if ((uint32_t)m >= threshold) {
            return (uint32_t)(m >> 32);
        }
    }
}

uint64_t ms_rng_bounded64(MsRng *rng, uint64_t upper_exclusive)
{
    uint64_t threshold;
    __uint128_t product;
    if (upper_exclusive <= 1) {
        return 0;
    }
    threshold = (uint64_t)(-upper_exclusive) % upper_exclusive;
    for (;;) {
        product = (__uint128_t)ms_rng_next64(rng) * upper_exclusive;
        if ((uint64_t)product >= threshold) {
            return (uint64_t)(product >> 64);
        }
    }
}

int ms_rng_dump_test(uint64_t seed, unsigned long count)
{
    MsRng rng;
    unsigned long i;
    ms_rng_seed(&rng, seed);
    for (i = 0; i < count; i++) {
        if (i != 0) {
            putchar(' ');
        }
        printf("%u", ms_rng_bounded32(&rng, 100000));
    }
    putchar('\n');
    return 0;
}

int ms_rng_pattern_test(uint64_t seed, uint32_t first_upper,
                        uint32_t second_upper)
{
    MsRng rng;
    uint32_t upper[2] = {first_upper, second_upper};
    int stage, i, j;
    ms_rng_seed(&rng, seed);
    for (stage = 0; stage < 2; stage++) {
        for (i = 0; i < 1500; i++) {
            for (j = 0; j < 3; j++) {
                if (stage != 0 || i != 0 || j != 0) {
                    putchar(' ');
                }
                printf("%u", ms_rng_bounded32(&rng, upper[stage]));
            }
        }
    }
    putchar('\n');
    return 0;
}

void ms_predict(const MsModel *model, double x, double y, double *dr, double *da)
{
    *dr = model->r.c[0] + model->r.c[1] * x + model->r.c[2] * y;
    *da = model->a.c[0] + model->a.c[1] * x + model->a.c[2] * y;
}

void ms_prm_values(const MsModel *model, char out[8][64])
{
    double fr = floor(model->r.c[0]);
    double fa = floor(model->a.c[0]);
    snprintf(out[0], 64, "%.0f", fr);
    snprintf(out[1], 64, "%.6f", model->r.c[0] - fr);
    snprintf(out[2], 64, "%.9g", model->r.c[1]);
    snprintf(out[3], 64, "%.9g", model->r.c[2]);
    snprintf(out[4], 64, "%.0f", fa);
    snprintf(out[5], 64, "%.6f", model->a.c[0] - fa);
    snprintf(out[6], 64, "%.9g", model->a.c[1]);
    snprintf(out[7], 64, "%.9g", model->a.c[2]);
}

static int cmp_double_fit(const void *a, const void *b)
{
    double x = *(const double *)a;
    double y = *(const double *)b;
    return (x > y) - (x < y);
}

double ms_mad(const double *v, size_t n)
{
    double med;
    double *dev;
    size_t i;
    if (n == 0) return NAN;
    med = ms_median(v, n);
    dev = malloc(n * sizeof(*dev));
    if (dev == NULL) return NAN;
    for (i = 0; i < n; i++) dev[i] = fabs(v[i] - med);
    med = ms_median(dev, n);
    free(dev);
    return med;
}

void ms_residual_stats(const double *v, size_t n, double *mad, double *rmse)
{
    size_t i;
    double sum = 0.0;
    *mad = n >= 4 ? ms_mad(v, n) : NAN;
    if (n < 4) {
        *rmse = NAN;
        return;
    }
    for (i = 0; i < n; i++) sum += v[i] * v[i];
    *rmse = sqrt(sum / (double)n);
}

static int solve3(double a[3][3], double b[3][2], double out[3][2])
{
    int i, j, k, p;
    double m[3][5];
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) m[i][j] = a[i][j];
        for (j = 0; j < 2; j++) m[i][3 + j] = b[i][j];
    }
    for (k = 0; k < 3; k++) {
        double best = fabs(m[k][k]);
        p = k;
        for (i = k + 1; i < 3; i++) {
            if (fabs(m[i][k]) > best) {
                best = fabs(m[i][k]);
                p = i;
            }
        }
        if (best < 1e-14) return -1;
        if (p != k) {
            for (j = k; j < 5; j++) {
                double t = m[k][j]; m[k][j] = m[p][j]; m[p][j] = t;
            }
        }
        for (i = k + 1; i < 3; i++) {
            double f = m[i][k] / m[k][k];
            for (j = k; j < 5; j++) m[i][j] -= f * m[k][j];
        }
    }
    for (j = 0; j < 2; j++) {
        for (i = 2; i >= 0; i--) {
            double v = m[i][3 + j];
            for (k = i + 1; k < 3; k++) v -= m[i][k] * out[k][j];
            out[i][j] = v / m[i][i];
        }
    }
    return 0;
}

static int least_squares(const double *x, const double *y, const double *dr,
                         const double *da, const double *weights, size_t n,
                         double out[3][2])
{
    double *a = calloc(n * 3, sizeof(*a));
    double *b = calloc(n * 2, sizeof(*b));
    double *v = malloc(n * sizeof(*v));
    size_t i;
    int j, k;
    if (a == NULL || b == NULL || v == NULL) {
        free(a); free(b); free(v); return -1;
    }
    for (i = 0; i < n; i++) {
        double row[3] = {1.0, x[i], y[i]};
        double w = weights == NULL ? 1.0 : weights[i];
        double sw = sqrt(w);
        for (j = 0; j < 3; j++) {
            a[3 * i + j] = sw * row[j];
        }
        b[2 * i] = sw * dr[i];
        b[2 * i + 1] = sw * da[i];
    }
    for (k = 0; k < 3; k++) {
        double norm = 0.0, alpha, vv = 0.0, beta;
        for (i = (size_t)k; i < n; i++) {
            double z = a[3 * i + k];
            norm += z * z;
        }
        norm = sqrt(norm);
        if (norm < 1e-14) { free(a); free(b); free(v); return -1; }
        alpha = -copysign(norm, a[3 * (size_t)k + k]);
        for (i = 0; i < n; i++) v[i] = i < (size_t)k ? 0.0 : a[3 * i + k];
        v[k] -= alpha;
        for (i = (size_t)k; i < n; i++) vv += v[i] * v[i];
        if (vv < 1e-28) { free(a); free(b); free(v); return -1; }
        beta = 2.0 / vv;
        for (j = k; j < 3; j++) {
            double d = 0.0;
            for (i = (size_t)k; i < n; i++) d += v[i] * a[3 * i + j];
            d *= beta;
            for (i = (size_t)k; i < n; i++) a[3 * i + j] -= d * v[i];
        }
        for (j = 0; j < 2; j++) {
            double d = 0.0;
            for (i = (size_t)k; i < n; i++) d += v[i] * b[2 * i + j];
            d *= beta;
            for (i = (size_t)k; i < n; i++) b[2 * i + j] -= d * v[i];
        }
        a[3 * (size_t)k + k] = alpha;
    }
    for (j = 0; j < 2; j++) {
        for (i = 3; i-- > 0;) {
            double z = b[2 * i + (size_t)j];
            size_t r;
            for (r = i + 1; r < 3; r++) z -= a[3 * i + r] * out[r][j];
            out[i][j] = z / a[3 * i + i];
        }
    }
    free(a); free(b); free(v);
    return 0;
}

static double det3(const double a[3][3])
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
           a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
           a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

static double median_copy(const double *v, size_t n)
{
    double *copy = malloc(n * sizeof(*copy));
    double result;
    if (copy == NULL) return NAN;
    memcpy(copy, v, n * sizeof(*copy));
    qsort(copy, n, sizeof(*copy), cmp_double_fit);
    if (n & 1U) result = copy[n / 2];
    else result = 0.5 * (copy[n / 2 - 1] + copy[n / 2]);
    free(copy);
    return result;
}

int ms_fit_affine(const MsPoint *points, size_t n, double threshold,
                  MsRng *rng, MsModel *model, unsigned char *inliers,
                  MsStage *stats)
{
    const int n_iter = 1500;
    double *x, *y, *dr, *da, *res;
    double xc = 0.0, yc = 0.0, scale, coef[3][2] = {{0}};
    double best_score = INFINITY;
    size_t i;
    int iter;
    double *weights;
    if (n < 6) return -1;
    x = malloc(n * sizeof(*x)); y = malloc(n * sizeof(*y));
    dr = malloc(n * sizeof(*dr)); da = malloc(n * sizeof(*da));
    res = malloc(n * sizeof(*res)); weights = malloc(n * sizeof(*weights));
    if (!x || !y || !dr || !da || !res || !weights) {
        free(x); free(y); free(dr); free(da); free(res); free(weights);
        return -1;
    }
    for (i = 0; i < n; i++) {
        x[i] = points[i].x;
        y[i] = points[i].y;
        dr[i] = points[i].dr;
        da[i] = points[i].da;
        xc += x[i];
        yc += y[i];
    }
    xc /= (double)n;
    yc /= (double)n;
    {
        double minx=x[0], maxx=x[0], miny=y[0], maxy=y[0];
        for (i = 1; i < n; i++) {
            if (x[i] < minx) minx=x[i];
            if (x[i] > maxx) maxx=x[i];
            if (y[i] < miny) miny=y[i];
            if (y[i] > maxy) maxy=y[i];
        }
        scale = fmax(fmax(maxx-minx, maxy-miny), 1.0) / 2.0;
    }
    for (i = 0; i < n; i++) {
        x[i] = (x[i] - xc) / scale;
        y[i] = (y[i] - yc) / scale;
    }
    for (iter = 0; iter < n_iter; iter++) {
        size_t idx[3];
        double a[3][3], b[3][2], c[3][2];
        double det;
        double score = 0.0;
        size_t j;
        for (j = 0; j < 3; j++) {
            idx[j] = ms_rng_bounded32(rng, (uint32_t)n);
            a[j][0] = 1.0; a[j][1] = x[idx[j]]; a[j][2] = y[idx[j]];
            b[j][0] = dr[idx[j]]; b[j][1] = da[idx[j]];
        }
        det = det3(a);
        if (fabs(det) <= 0.05 || solve3(a, b, c) != 0) continue;
        for (j = 0; j < n; j++) {
            double pr = c[0][0] + c[1][0] * x[j] + c[2][0] * y[j];
            double pa = c[0][1] + c[1][1] * x[j] + c[2][1] * y[j];
            double rr = dr[j] - pr, ra = da[j] - pa;
            double r2 = rr * rr + ra * ra;
            score += fmin(r2, threshold * threshold);
        }
        if (score < best_score) {
            best_score = score;
            memcpy(coef, c, sizeof(coef));
        }
    }
    if (!isfinite(best_score)) goto fail;
    for (iter = 0; iter < 10; iter++) {
        double *tmp_res = malloc(n * sizeof(*tmp_res));
        double s;
        double next[3][2];
        double max_delta = 0.0;
        if (tmp_res == NULL) goto fail;
        for (i = 0; i < n; i++) {
            double rr = dr[i] - (coef[0][0] + coef[1][0]*x[i] + coef[2][0]*y[i]);
            double ra = da[i] - (coef[0][1] + coef[1][1]*x[i] + coef[2][1]*y[i]);
            res[i] = hypot(rr, ra);
            tmp_res[i] = res[i];
        }
        s = fmax(median_copy(tmp_res, n) / 1.1774, 0.2 / 3.0);
        free(tmp_res);
        for (i = 0; i < n; i++) {
            double u = res[i] / (4.685 * s);
            weights[i] = u < 1.0 ? (1.0-u*u)*(1.0-u*u) : 0.0;
        }
        {
            size_t nw = 0;
            for (i = 0; i < n; i++) if (weights[i] > 0.0) nw++;
            if (nw < 4) break;
        }
        if (least_squares(x, y, dr, da, weights, n, next) != 0) goto fail;
        for (i = 0; i < 3; i++) {
            int j;
            for (j = 0; j < 2; j++) {
                double d = fabs(next[i][j] - coef[i][j]);
                if (d > max_delta) max_delta = d;
            }
        }
        memcpy(coef, next, sizeof(coef));
        if (max_delta < 1e-9) break;
    }
    {
        double *abs_pool = malloc(2 * n * sizeof(*abs_pool));
        double *tmp_res = malloc(n * sizeof(*tmp_res));
        double sigma;
        double tin;
        size_t npool = 0;
        size_t nin = 0;
        double next[3][2];
        if (abs_pool == NULL || tmp_res == NULL) {
            free(abs_pool); free(tmp_res); goto fail;
        }
        for (i = 0; i < n; i++) {
            double rr = dr[i] - (coef[0][0] + coef[1][0]*x[i] + coef[2][0]*y[i]);
            double ra = da[i] - (coef[0][1] + coef[1][1]*x[i] + coef[2][1]*y[i]);
            res[i] = hypot(rr, ra);
            if (res[i] <= threshold) {
                abs_pool[npool++] = fabs(rr);
                abs_pool[npool++] = fabs(ra);
                nin++;
            }
        }
        if (nin < 4) { free(abs_pool); free(tmp_res); goto fail; }
        sigma = 1.4826 * median_copy(abs_pool, npool);
        tin = fmin(threshold, fmax(3.0 * sigma, 0.2));
        free(abs_pool);
        for (i = 0; i < n; i++) weights[i] = res[i] <= tin ? 1.0 : 0.0;
        nin = 0; for (i=0;i<n;i++) if (weights[i] > 0) nin++;
        if (nin < 4 || least_squares(x,y,dr,da,weights,n,next) != 0) { free(tmp_res); goto fail; }
        memcpy(coef,next,sizeof(coef));
        nin = 0;
        for (i = 0; i < n; i++) {
            double rr = dr[i] - (coef[0][0] + coef[1][0]*x[i] + coef[2][0]*y[i]);
            double ra = da[i] - (coef[0][1] + coef[1][1]*x[i] + coef[2][1]*y[i]);
            res[i] = hypot(rr,ra);
            weights[i] = res[i] <= tin ? 1.0 : 0.0;
            if (weights[i] > 0) nin++;
        }
        if (nin < 4 || least_squares(x,y,dr,da,weights,n,next) != 0) { free(tmp_res); goto fail; }
        memcpy(coef,next,sizeof(coef));
        {
            double sse_r=0.0, sse_a=0.0;
            double *rrv=malloc(nin*sizeof(*rrv)), *rav=malloc(nin*sizeof(*rav));
            size_t m=0;
            if (!rrv || !rav) { free(rrv); free(rav); free(tmp_res); goto fail; }
            for (i=0;i<n;i++) {
                double rr=dr[i]-(coef[0][0]+coef[1][0]*x[i]+coef[2][0]*y[i]);
                double ra=da[i]-(coef[0][1]+coef[1][1]*x[i]+coef[2][1]*y[i]);
                if (weights[i] > 0.0) {
                    rrv[m]=rr; rav[m]=ra; m++;
                    sse_r += rr*rr; sse_a += ra*ra;
                }
            }
            model->r.c[1]=coef[1][0]/scale;
            model->r.c[2]=coef[2][0]/scale;
            model->r.c[0]=coef[0][0]-(coef[1][0]*xc+coef[2][0]*yc)/scale;
            model->a.c[1]=coef[1][1]/scale;
            model->a.c[2]=coef[2][1]/scale;
            model->a.c[0]=coef[0][1]-(coef[1][1]*xc+coef[2][1]*yc)/scale;
            for (i=0;i<n;i++) {
                double pr,pa; ms_predict(model,points[i].x,points[i].y,&pr,&pa);
                inliers[i]=(unsigned char)(hypot(points[i].dr-pr,points[i].da-pa)<=tin);
            }
            stats->sigma=sigma; stats->t_in=tin;
            stats->rms_r=sqrt(sse_r/(double)m); stats->rms_a=sqrt(sse_a/(double)m);
            stats->metric_n_points=(int)m;
            stats->mad_r=ms_mad(rrv,m); stats->mad_a=ms_mad(rav,m);
            stats->robust_sigma_r=1.4826*stats->mad_r; stats->robust_sigma_a=1.4826*stats->mad_a;
            stats->rmse_r=stats->rms_r; stats->rmse_a=stats->rms_a;
            free(rrv); free(rav);
        }
        free(tmp_res);
    }
    free(x); free(y); free(dr); free(da); free(res); free(weights);
    return 0;
fail:
    free(x); free(y); free(dr); free(da); free(res); free(weights);
    return -1;
}
