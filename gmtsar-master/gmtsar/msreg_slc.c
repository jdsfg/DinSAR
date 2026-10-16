#define _POSIX_C_SOURCE 200809L
#include "msreg.h"

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

static char *ms_strdup(const char *s)
{
    size_t n = strlen(s) + 1;
    char *p = malloc(n);
    if (p != NULL) {
        memcpy(p, s, n);
    }
    return p;
}

static char *trim(char *s)
{
    char *end;
    while (isspace((unsigned char)*s)) {
        s++;
    }
    end = s + strlen(s);
    while (end > s && isspace((unsigned char)end[-1])) {
        *--end = '\0';
    }
    return s;
}

static int append_string(char ***items, size_t *count, const char *value)
{
    char **next = realloc(*items, (*count + 1) * sizeof(**items));
    if (next == NULL) {
        return -1;
    }
    *items = next;
    (*items)[*count] = ms_strdup(value);
    if ((*items)[*count] == NULL) {
        return -1;
    }
    (*count)++;
    return 0;
}

int ms_prm_read(const char *path, MsPrm *p, MsError *err)
{
    FILE *f;
    char *line = NULL;
    size_t cap = 0;
    ssize_t nread;

    memset(p, 0, sizeof(*p));
    p->path = ms_strdup(path);
    if (p->path == NULL) {
        snprintf(err->message, sizeof(err->message), "out of memory");
        return -1;
    }
    f = fopen(path, "r");
    if (f == NULL) {
        snprintf(err->message, sizeof(err->message), "PRM does not exist: %s: %s",
                 path, strerror(errno));
        return -1;
    }
    while ((nread = getline(&line, &cap, f)) >= 0) {
        char *copy = malloc((size_t)nread + 1);
        char *eq;
        if (copy == NULL) {
            snprintf(err->message, sizeof(err->message), "out of memory");
            free(line);
            fclose(f);
            return -1;
        }
        memcpy(copy, line, (size_t)nread);
        copy[nread] = '\0';
        if (p->nlines == 0) {
            p->lines = malloc(sizeof(*p->lines));
        } else {
            p->lines = realloc(p->lines, (p->nlines + 1) * sizeof(*p->lines));
        }
        if (p->lines == NULL) {
            free(copy);
            free(line);
            fclose(f);
            snprintf(err->message, sizeof(err->message), "out of memory");
            return -1;
        }
        p->lines[p->nlines++] = copy;
        eq = strchr(line, '=');
        if (eq != NULL) {
            char *key;
            char *value;
            *eq = '\0';
            key = trim(line);
            value = trim(eq + 1);
            if (*key != '\0') {
                if (append_string(&p->keys, &p->nfields, key) != 0) {
                    free(line);
                    fclose(f);
                    snprintf(err->message, sizeof(err->message), "out of memory");
                    return -1;
                }
                {
                    char **next = realloc(p->values, p->nfields * sizeof(*p->values));
                    if (next == NULL) {
                        free(line);
                        fclose(f);
                        snprintf(err->message, sizeof(err->message), "out of memory");
                        return -1;
                    }
                    p->values = next;
                    p->values[p->nfields - 1] = ms_strdup(value);
                    if (p->values[p->nfields - 1] == NULL) {
                        free(line);
                        fclose(f);
                        snprintf(err->message, sizeof(err->message), "out of memory");
                        return -1;
                    }
                }
            }
            *eq = '=';
        }
    }
    free(line);
    if (ferror(f)) {
        snprintf(err->message, sizeof(err->message), "failed reading PRM: %s", path);
        fclose(f);
        return -1;
    }
    fclose(f);
    return 0;
}

void ms_prm_free(MsPrm *p)
{
    size_t i;
    for (i = 0; i < p->nlines; i++) {
        free(p->lines[i]);
    }
    for (i = 0; i < p->nfields; i++) {
        free(p->keys[i]);
        free(p->values[i]);
    }
    free(p->lines);
    free(p->keys);
    free(p->values);
    free(p->path);
    memset(p, 0, sizeof(*p));
}

const char *ms_prm_get(const MsPrm *p, const char *key)
{
    size_t i = p->nfields;
    while (i > 0) {
        i--;
        if (strcmp(p->keys[i], key) == 0) {
            return p->values[i];
        }
    }
    return NULL;
}

static int parse_int_field(const MsPrm *p, const char *key, int fallback,
                           int *out)
{
    const char *v = ms_prm_get(p, key);
    char *end;
    double d;
    if (v == NULL) {
        *out = fallback;
        return 0;
    }
    errno = 0;
    d = strtod(v, &end);
    if (errno != 0 || end == v || !isfinite(d) || d < (double)INT32_MIN ||
        d > (double)INT32_MAX) {
        return -1;
    }
    *out = (int)d;
    return 0;
}

static char *resolve_slc_path(const char *prm_path, const char *name)
{
    struct stat st;
    char candidate[4096];
    const char *slash;
    if (name[0] == '/') {
        return stat(name, &st) == 0 ? ms_strdup(name) : NULL;
    }
    slash = strrchr(prm_path, '/');
    if (slash != NULL) {
        size_t prefix = (size_t)(slash - prm_path + 1);
        if (prefix + strlen(name) + 1 < sizeof(candidate)) {
            memcpy(candidate, prm_path, prefix);
            strcpy(candidate + prefix, name);
            if (stat(candidate, &st) == 0) {
                return ms_strdup(candidate);
            }
        }
    }
    if (stat(name, &st) == 0) {
        return ms_strdup(name);
    }
    return NULL;
}

int ms_slc_open(MsSlc *s, const char *prm_path, const char *slc_override,
                MsError *err)
{
    const char *name;
    struct stat st;
    int patches = 1;
    int valid_az = 0;
    size_t row_bytes;
    size_t rows;

    memset(s, 0, sizeof(*s));
    s->fd = -1;
    s->prm_path = ms_strdup(prm_path);
    if (s->prm_path == NULL || ms_prm_read(prm_path, &s->prm, err) != 0) {
        if (s->prm_path == NULL) {
            snprintf(err->message, sizeof(err->message), "out of memory");
        }
        return -1;
    }
    if (parse_int_field(&s->prm, "num_rng_bins", 0, &s->nx) != 0 || s->nx <= 0) {
        snprintf(err->message, sizeof(err->message), "PRM missing valid num_rng_bins: %s",
                 prm_path);
        return -1;
    }
    name = slc_override != NULL ? slc_override : ms_prm_get(&s->prm, "SLC_file");
    if (name == NULL || *name == '\0') {
        snprintf(err->message, sizeof(err->message), "PRM missing SLC_file: %s", prm_path);
        return -1;
    }
    s->path = resolve_slc_path(prm_path, name);
    if (s->path == NULL) {
        snprintf(err->message, sizeof(err->message), "cannot find SLC file: %s (from %s)",
                 name, prm_path);
        return -1;
    }
    if (stat(s->path, &st) != 0 || st.st_size < 0) {
        snprintf(err->message, sizeof(err->message), "cannot stat SLC file: %s", s->path);
        return -1;
    }
    row_bytes = (size_t)4 * (size_t)s->nx;
    if (row_bytes == 0) {
        snprintf(err->message, sizeof(err->message), "invalid SLC row size");
        return -1;
    }
    rows = (size_t)st.st_size / row_bytes;
    if (parse_int_field(&s->prm, "num_patches", 1, &patches) != 0 ||
        parse_int_field(&s->prm, "num_valid_az", 0, &valid_az) != 0) {
        snprintf(err->message, sizeof(err->message), "invalid PRM azimuth dimensions");
        return -1;
    }
    if (valid_az > 0 && patches > 0 && (size_t)patches <= SIZE_MAX / (size_t)valid_az) {
        size_t prm_rows = (size_t)patches * (size_t)valid_az;
        if (rows > prm_rows) {
            rows = prm_rows;
        }
    }
    if (rows > INT32_MAX || rows < 64 || s->nx < 64) {
        snprintf(err->message, sizeof(err->message), "SLC dimensions too small or too large: %s",
                 s->path);
        return -1;
    }
    s->ny = (int)rows;
    s->map_size = (size_t)st.st_size;
    s->fd = open(s->path, O_RDONLY);
    if (s->fd < 0) {
        snprintf(err->message, sizeof(err->message), "cannot open SLC: %s: %s",
                 s->path, strerror(errno));
        return -1;
    }
    s->samples = mmap(NULL, s->map_size, PROT_READ, MAP_PRIVATE, s->fd, 0);
    if (s->samples == MAP_FAILED) {
        s->samples = NULL;
        snprintf(err->message, sizeof(err->message), "cannot mmap SLC: %s: %s",
                 s->path, strerror(errno));
        close(s->fd);
        s->fd = -1;
        return -1;
    }
    return 0;
}

void ms_slc_close(MsSlc *s)
{
    if (s->samples != NULL) {
        munmap((void *)s->samples, s->map_size);
    }
    if (s->fd >= 0) {
        close(s->fd);
    }
    free(s->prm_path);
    free(s->path);
    ms_prm_free(&s->prm);
    memset(s, 0, sizeof(*s));
    s->fd = -1;
}

void ms_slc_read(const MsSlc *s, int y0, int x0, int h, int w, float *out_ri)
{
    int y;
    memset(out_ri, 0, (size_t)h * (size_t)w * 2 * sizeof(*out_ri));
    for (y = 0; y < h; y++) {
        int sy = y0 + y;
        int xa;
        int xb;
        int x;
        const int16_t *row;
        float *dst;
        if (sy < 0 || sy >= s->ny) {
            continue;
        }
        xa = x0 < 0 ? -x0 : 0;
        xb = w;
        if (x0 + xb > s->nx) {
            xb = s->nx - x0;
        }
        if (xa >= xb) {
            continue;
        }
        row = s->samples + (size_t)sy * (size_t)(2 * s->nx);
        dst = out_ri + (size_t)y * (size_t)(2 * w);
        for (x = xa; x < xb; x++) {
            int sx = x0 + x;
            dst[2 * x] = (float)row[2 * sx];
            dst[2 * x + 1] = (float)row[2 * sx + 1];
        }
    }
}

int ms_multilook(const MsSlc *s, int looks, MsImage *out, MsError *err)
{
    int y, x, ay, ax;
    int h;
    int w;
    if (looks < 1) {
        snprintf(err->message, sizeof(err->message), "looks must be positive");
        return -1;
    }
    h = s->ny / looks;
    w = s->nx / looks;
    if (h <= 0 || w <= 0 || (size_t)h > SIZE_MAX / (size_t)w) {
        snprintf(err->message, sizeof(err->message), "invalid multilook dimensions");
        return -1;
    }
    out->h = h;
    out->w = w;
    out->v = malloc((size_t)h * (size_t)w * sizeof(*out->v));
    if (out->v == NULL) {
        snprintf(err->message, sizeof(err->message), "out of memory for multilook image");
        return -1;
    }
    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            float sum = 0.0f;
            for (ay = 0; ay < looks; ay++) {
                const int16_t *row = s->samples +
                    (size_t)(y * looks + ay) * (size_t)(2 * s->nx);
                float rsum = 0.0f;
                for (ax = 0; ax < looks; ax++) {
                    float ri = (float)row[2 * (x * looks + ax)];
                    float rq = (float)row[2 * (x * looks + ax) + 1];
                    float vi = ri * ri;
                    float vq = rq * rq;
                    rsum += vi + vq;
                }
                sum += rsum;
            }
            out->v[(size_t)y * (size_t)w + (size_t)x] =
                sum / (float)(looks * looks);
        }
    }
    return 0;
}

void ms_box_mean(const float *src, int h, int w, int k, float *dst)
{
    int lo;
    int hi;
    int ph;
    int pw;
    int y, x;
    double *padded;
    double *ycum;
    double *vsum;
    double *xcum;
    if (k <= 1) {
        memcpy(dst, src, (size_t)h * (size_t)w * sizeof(*dst));
        return;
    }
    lo = k / 2;
    hi = k - 1 - lo;
    ph = h + lo + hi;
    pw = w + lo + hi;
    padded = malloc((size_t)ph * (size_t)pw * sizeof(*padded));
    ycum = calloc((size_t)(ph + 1) * (size_t)pw, sizeof(*ycum));
    vsum = malloc((size_t)h * (size_t)pw * sizeof(*vsum));
    xcum = malloc((size_t)h * (size_t)(pw + 1) * sizeof(*xcum));
    if (padded == NULL || ycum == NULL || vsum == NULL || xcum == NULL) {
        for (y = 0; y < h; y++) {
            for (x = 0; x < w; x++) {
                double sum = 0.0;
                int dy, dx;
                for (dy = -lo; dy <= hi; dy++) {
                    int yy = y + dy;
                    if (yy < 0) yy = 0;
                    if (yy >= h) yy = h - 1;
                    for (dx = -lo; dx <= hi; dx++) {
                        int xx = x + dx;
                        if (xx < 0) xx = 0;
                        if (xx >= w) xx = w - 1;
                        sum += src[(size_t)yy * (size_t)w + (size_t)xx];
                    }
                }
                dst[(size_t)y * (size_t)w + (size_t)x] =
                    (float)(sum / (double)(k * k));
            }
        }
        free(padded);
        free(ycum);
        free(vsum);
        free(xcum);
        return;
    }
    for (y = 0; y < ph; y++) {
        int yy = y - lo;
        if (yy < 0) yy = 0;
        if (yy >= h) yy = h - 1;
        for (x = 0; x < pw; x++) {
            int xx = x - lo;
            if (xx < 0) xx = 0;
            if (xx >= w) xx = w - 1;
            padded[(size_t)y * (size_t)pw + (size_t)x] =
                (double)src[(size_t)yy * (size_t)w + (size_t)xx];
            ycum[(size_t)(y + 1) * (size_t)pw + (size_t)x] =
                ycum[(size_t)y * (size_t)pw + (size_t)x] +
                padded[(size_t)y * (size_t)pw + (size_t)x];
        }
    }
    for (y = 0; y < h; y++) {
        for (x = 0; x < pw; x++) {
            vsum[(size_t)y * (size_t)pw + (size_t)x] =
                ycum[(size_t)(y + k) * (size_t)pw + (size_t)x] -
                ycum[(size_t)y * (size_t)pw + (size_t)x];
        }
    }
    for (y = 0; y < h; y++) {
        double *row = xcum + (size_t)y * (size_t)(pw + 1);
        row[0] = 0.0;
        for (x = 0; x < pw; x++) {
            row[x + 1] = row[x] + vsum[(size_t)y * (size_t)pw + (size_t)x];
        }
        for (x = 0; x < w; x++) {
            dst[(size_t)y * (size_t)w + (size_t)x] =
                (float)((row[x + k] - row[x]) / (double)(k * k));
        }
    }
    free(padded);
    free(ycum);
    free(vsum);
    free(xcum);
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a;
    double y = *(const double *)b;
    return (x > y) - (x < y);
}

double ms_median(const double *v, size_t n)
{
    double *copy;
    double result;
    if (n == 0) {
        return NAN;
    }
    copy = malloc(n * sizeof(*copy));
    if (copy == NULL) {
        return NAN;
    }
    memcpy(copy, v, n * sizeof(*copy));
    qsort(copy, n, sizeof(*copy), cmp_double);
    if ((n & 1U) != 0U) {
        result = copy[n / 2];
    } else {
        result = 0.5 * (copy[n / 2 - 1] + copy[n / 2]);
    }
    free(copy);
    return result;
}

static void format_model_values(const MsModel *model, char out[8][64])
{
    ms_prm_values(model, out);
}

int ms_prm_write(const MsPrm *src, const char *dst, const MsModel *model,
                 const char *slc_file_override, MsError *err)
{
    static const char *keys[9] = {
        "rshift", "sub_int_r", "stretch_r", "a_stretch_r",
        "ashift", "sub_int_a", "stretch_a", "a_stretch_a", "SLC_file"
    };
    char values[8][64];
    char tmp[4096];
    FILE *f;
    size_t i;
    int has_override = slc_file_override != NULL;
    const char *dst_slash;
    const char *src_slash;

    if (model != NULL) {
        format_model_values(model, values);
    }
    dst_slash = strrchr(dst, '/');
    src_slash = strrchr(src->path, '/');
    if (slc_file_override == NULL && dst_slash != NULL && src_slash != NULL &&
        (size_t)(dst_slash - dst) != (size_t)(src_slash - src->path)) {
        has_override = 1;
        slc_file_override = ms_prm_get(src, "SLC_file");
    }
    if (snprintf(tmp, sizeof(tmp), "%s.tmp", dst) >= (int)sizeof(tmp)) {
        snprintf(err->message, sizeof(err->message), "output PRM path too long");
        return -1;
    }
    f = fopen(tmp, "w");
    if (f == NULL) {
        snprintf(err->message, sizeof(err->message), "cannot write PRM: %.850s: %.120s",
                 tmp, strerror(errno));
        return -1;
    }
    for (i = 0; i < src->nlines; i++) {
        char *line = ms_strdup(src->lines[i]);
        char *eq;
        char *key;
        size_t k;
        int skip = 0;
        if (line == NULL) {
            fclose(f);
            unlink(tmp);
            snprintf(err->message, sizeof(err->message), "out of memory");
            return -1;
        }
        eq = strchr(line, '=');
        if (eq != NULL) {
            *eq = '\0';
            key = trim(line);
            for (k = 0; k < 8; k++) {
                if (model != NULL && strcmp(key, keys[k]) == 0) {
                    skip = 1;
                }
            }
            if (has_override && strcmp(key, "SLC_file") == 0) {
                skip = 1;
            }
        }
        free(line);
        if (!skip && fputs(src->lines[i], f) == EOF) {
            fclose(f);
            unlink(tmp);
            snprintf(err->message, sizeof(err->message), "failed writing PRM: %.950s", tmp);
            return -1;
        }
    }
    if (model != NULL) {
        for (i = 0; i < 8; i++) {
            if (fprintf(f, "%s = %s\n", keys[i], values[i]) < 0) {
                fclose(f);
                unlink(tmp);
                snprintf(err->message, sizeof(err->message), "failed writing PRM: %.950s", tmp);
                return -1;
            }
        }
    }
    if (has_override && slc_file_override != NULL &&
        fprintf(f, "SLC_file = %s\n", slc_file_override) < 0) {
        fclose(f);
        unlink(tmp);
        snprintf(err->message, sizeof(err->message), "failed writing PRM: %.950s", tmp);
        return -1;
    }
    if (fclose(f) != 0 || rename(tmp, dst) != 0) {
        snprintf(err->message, sizeof(err->message), "cannot finalize PRM: %s: %s",
                 dst, strerror(errno));
        unlink(tmp);
        return -1;
    }
    return 0;
}
