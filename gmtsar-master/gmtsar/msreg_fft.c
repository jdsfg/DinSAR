#define _POSIX_C_SOURCE 200809L
#include "msreg.h"

#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef MS_HAVE_FFTW_THREADS
#include <fftw3.h>
#endif

static int fast_len(int n)
{
    int m;
    for (m = n; ; m++) {
        int k = m;
        while (k % 2 == 0) k /= 2;
        while (k % 3 == 0) k /= 3;
        while (k % 5 == 0) k /= 5;
        if (k == 1) return m;
    }
}

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void setup_threads(int n)
{
#ifdef MS_HAVE_FFTW_THREADS
    static int initialized = 0;
    if (!initialized) {
        fftwf_init_threads();
        initialized = 1;
    }
    fftwf_plan_with_nthreads(n > 0 ? n : 1);
#else
    (void)n;
#endif
}

static size_t csize(int h, int w)
{
    return (size_t)h * (size_t)(w / 2 + 1);
}

static fftwf_complex *make_spec(const float *src, int h, int w, int fh, int fw)
{
    float *pad = fftwf_alloc_real((size_t)fh * (size_t)fw);
    fftwf_complex *out = fftwf_alloc_complex(csize(fh, fw));
    fftwf_plan p;
    int y;
    if (pad == NULL || out == NULL) {
        fftwf_free(pad); fftwf_free(out); return NULL;
    }
    memset(pad, 0, (size_t)fh * (size_t)fw * sizeof(*pad));
    for (y = 0; y < h; y++) {
        memcpy(pad + (size_t)y * (size_t)fw,
               src + (size_t)y * (size_t)w, (size_t)w * sizeof(*src));
    }
    p = fftwf_plan_dft_r2c_2d(fh, fw, pad, out, FFTW_ESTIMATE);
    if (p == NULL) {
        fftwf_free(pad); fftwf_free(out); return NULL;
    }
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    fftwf_free(pad);
    return out;
}

static int inverse_product(const fftwf_complex *a, const fftwf_complex *b,
                           int fh, int fw, float *out)
{
    size_t i, nc = csize(fh, fw);
    fftwf_complex *p = fftwf_alloc_complex(nc);
    fftwf_plan plan;
    if (p == NULL) return -1;
    for (i = 0; i < nc; i++) {
        p[i][0] = a[i][0] * b[i][0] + a[i][1] * b[i][1];
        p[i][1] = a[i][0] * b[i][1] - a[i][1] * b[i][0];
    }
    plan = fftwf_plan_dft_c2r_2d(fh, fw, p, out, FFTW_ESTIMATE);
    if (plan == NULL) {
        fftwf_free(p); return -1;
    }
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    {
        float scale = 1.0f / ((float)fh * (float)fw);
        size_t n = (size_t)fh * (size_t)fw;
        for (i = 0; i < n; i++) out[i] *= scale;
    }
    fftwf_free(p);
    return 0;
}

static float *corr(const fftwf_complex *a, const fftwf_complex *b, int h, int w)
{
    float *out = fftwf_alloc_real((size_t)h * (size_t)w);
    if (out == NULL || inverse_product(a, b, h, w, out) != 0) {
        fftwf_free(out); return NULL;
    }
    return out;
}

static double parabola(double l, double c, double r)
{
    double d = l - 2.0 * c + r;
    return d == 0.0 ? 0.0 : 0.5 * (l - r) / d;
}

static int round_even(double x)
{
    double lo=floor(x),f=x-lo;
    if(f<0.5)return (int)lo;
    if(f>0.5)return (int)(lo+1.0);
    return ((long long)lo&1LL)?(int)(lo+1.0):(int)lo;
}

static int wrap(int i, int n)
{
    i %= n;
    return i < 0 ? i + n : i;
}

static double second_peak(const float *ncc, int h, int w, int iy, int ix,
                          int radius)
{
    double best = -INFINITY;
    int y, x, dy, dx;
    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            float v = ncc[(size_t)y * (size_t)w + (size_t)x];
            int local = v > -1.5f;
            int yy = abs(((y - iy + h / 2) % h + h) % h - h / 2);
            int xx = abs(((x - ix + w / 2) % w + w) % w - w / 2);
            if (!local || (yy <= radius && xx <= radius)) continue;
            for (dy = -1; dy <= 1 && local; dy++) {
                for (dx = -1; dx <= 1; dx++) {
                    float q;
                    if (dy == 0 && dx == 0) continue;
                    q = ncc[(size_t)wrap(y + dy, h) * (size_t)w +
                            (size_t)wrap(x + dx, w)];
                    if (v < q) { local = 0; break; }
                }
            }
            if (local && v > best) best = v;
        }
    }
    return best;
}

static int global_xcorr(const MsImage *a, const MsImage *ma,
                       const MsImage *b, const MsImage *mb,
                       int fh, int fw, float min_pixels,
                       float **ncc_out)
{
    fftwf_complex *fam=NULL,*fbm=NULL,*fa=NULL,*fb=NULL;
    float *nraw=NULL,*sa=NULL,*sb=NULL,*sab=NULL,*saa=NULL,*sbb=NULL,*ncc=NULL;
    size_t i, n=(size_t)fh*(size_t)fw;
    int ok=-1;
    fam=make_spec(ma->v,ma->h,ma->w,fh,fw);
    fbm=make_spec(mb->v,mb->h,mb->w,fh,fw);
    fa=make_spec(a->v,a->h,a->w,fh,fw);
    fb=make_spec(b->v,b->h,b->w,fh,fw);
    if (!fam||!fbm||!fa||!fb) goto done;
    nraw=corr(fam,fbm,fh,fw);
    sa=corr(fa,fbm,fh,fw);
    sb=corr(fam,fb,fh,fw);
    sab=corr(fa,fb,fh,fw);
    if (!nraw||!sa||!sb||!sab) goto done;
    /* Squared images use the same validity mask. */
    {
        float *aa=malloc((size_t)a->h*(size_t)a->w*sizeof(*aa));
        float *bb=malloc((size_t)b->h*(size_t)b->w*sizeof(*bb));
        fftwf_complex *fa_sq,*fb_sq;
        if (!aa||!bb) { free(aa); free(bb); goto done; }
        for (i=0;i<(size_t)a->h*(size_t)a->w;i++) aa[i]=a->v[i]*a->v[i]*ma->v[i];
        for (i=0;i<(size_t)b->h*(size_t)b->w;i++) bb[i]=b->v[i]*b->v[i]*mb->v[i];
        fa_sq=make_spec(aa,a->h,a->w,fh,fw);
        fb_sq=make_spec(bb,b->h,b->w,fh,fw);
        free(aa); free(bb);
        if (!fa_sq||!fb_sq) { fftwf_free(fa_sq); fftwf_free(fb_sq); goto done; }
        saa=corr(fa_sq,fbm,fh,fw);
        sbb=corr(fam,fb_sq,fh,fw);
        fftwf_free(fa_sq); fftwf_free(fb_sq);
    }
    if (!saa||!sbb) goto done;
    ncc=malloc(n*sizeof(*ncc));
    if (!ncc) goto done;
    {
        int nb=a->h*a->w < b->h*b->w ? a->h*a->w : b->h*b->w;
        for (i=0;i<n;i++) {
            float count=fmaxf(nraw[i],1.0f);
            float num=sab[i]-sa[i]*sb[i]/count;
            float va=fmaxf(saa[i]-sa[i]*sa[i]/count,0.0f);
            float vb=fmaxf(sbb[i]-sb[i]*sb[i]/count,0.0f);
            float den=sqrtf(va*vb);
            ncc[i]=(nraw[i]>=min_pixels && den>1e-6f) ? num/den : -2.0f;
        }
        (void)nb;
    }
    *ncc_out=ncc; ncc=NULL; ok=0;
done:
    fftwf_free(fam); fftwf_free(fbm); fftwf_free(fa); fftwf_free(fb);
    fftwf_free(nraw); fftwf_free(sa); fftwf_free(sb); fftwf_free(sab);
    fftwf_free(saa); fftwf_free(sbb); free(ncc);
    return ok;
}

static int log_highpass(const MsImage *power, int hpbox, MsImage *hp,
                        MsImage *mask, MsError *err)
{
    size_t n=(size_t)power->h*(size_t)power->w,i,nvalid=0;
    double *v=malloc(n*sizeof(*v));
    float *lg=malloc(n*sizeof(*lg));
    float *mean=malloc(n*sizeof(*mean));
    double med,floorv;
    if (!v||!lg||!mean) { free(v);free(lg);free(mean);snprintf(err->message,sizeof(err->message),"out of memory");return -1; }
    for(i=0;i<n;i++) if(power->v[i]>0.0f) v[nvalid++]=power->v[i];
    if(nvalid==0){free(v);free(lg);free(mean);snprintf(err->message,sizeof(err->message),"image is all zero");return -1;}
    med=ms_median(v,nvalid); floorv=1e-3*med; free(v);
    hp->h=mask->h=power->h; hp->w=mask->w=power->w;
    hp->v=malloc(n*sizeof(*hp->v));mask->v=malloc(n*sizeof(*mask->v));
    if(!hp->v||!mask->v){free(lg);free(mean);free(hp->v);free(mask->v);hp->v=NULL;mask->v=NULL;snprintf(err->message,sizeof(err->message),"out of memory");return -1;}
    for(i=0;i<n;i++){float z=power->v[i];mask->v[i]=z>0?1.0f:0.0f;lg[i]=z>0?logf(fmaxf(z,(float)floorv)):0.0f;}
    ms_box_mean(lg,power->h,power->w,hpbox,mean);
    for(i=0;i<n;i++) hp->v[i]=mask->v[i]>0?lg[i]-mean[i]:0.0f;
    free(lg);free(mean);return 0;
}

int ms_global_stage(const MsSlc *master,const MsSlc *slave,int target,
                    int threads,MsGlobal *out,MsError *err)
{
    MsImage pa={0},pb={0},a={0},b={0},ma={0},mb={0};
    float *ncc=NULL;
    int big,L,fh,fw,iy=0,ix=0,y,x;
    double t0,tread,peak=-INFINITY,side_sum=0.0,side_sq=0.0,psr;
    size_t side_n=0;
    int ok=-1;
    memset(out,0,sizeof(*out));
    t0=now_s();
    big=master->nx;if(master->ny>big)big=master->ny;if(slave->nx>big)big=slave->nx;if(slave->ny>big)big=slave->ny;
    L=(int)ceil((double)big/(double)target);if(L<1)L=1;
    if(ms_multilook(master,L,&pa,err)||ms_multilook(slave,L,&pb,err))goto done;
    if(log_highpass(&pa,15,&a,&ma,err)||log_highpass(&pb,15,&b,&mb,err))goto done;
    tread=now_s()-t0;
    fh=fast_len(a.h+b.h-1);fw=fast_len(a.w+b.w-1);
    setup_threads(threads);
    {
        float minpix=0.25f*(float)((size_t)a.h*(size_t)a.w<(size_t)b.h*(size_t)b.w?(size_t)a.h*(size_t)a.w:(size_t)b.h*(size_t)b.w);
        if(global_xcorr(&a,&ma,&b,&mb,fh,fw,minpix,&ncc)!=0)goto done;
    }
    for(y=0;y<fh;y++)for(x=0;x<fw;x++){size_t k=(size_t)y*(size_t)fw+(size_t)x;if(ncc[k]>peak){peak=ncc[k];iy=y;ix=x;}}
    out->ok=0;out->looks=L;out->ncc=peak;
    if(peak<=-1.5){snprintf(out->reason,sizeof(out->reason),"insufficient overlap");goto success;}
    {
        int ly=iy<=b.h-1?iy:iy-fh;
        int lx=ix<=b.w-1?ix:ix-fw;
        double sy=0,sx=0;
        int valid=1,dy,dx;
        for(dy=-1;dy<=1;dy++)for(dx=-1;dx<=1;dx++)if(ncc[(size_t)wrap(iy+dy,fh)*(size_t)fw+(size_t)wrap(ix+dx,fw)]<=-1.5f)valid=0;
        if(valid){
            sy=parabola(ncc[(size_t)wrap(iy-1,fh)*fw+ix],ncc[(size_t)iy*fw+ix],ncc[(size_t)wrap(iy+1,fh)*fw+ix]);
            sx=parabola(ncc[(size_t)iy*fw+wrap(ix-1,fw)],ncc[(size_t)iy*fw+ix],ncc[(size_t)iy*fw+wrap(ix+1,fw)]);
        }
        for(y=0;y<fh;y++)for(x=0;x<fw;x++){
            int yy=abs(((y-iy+fh/2)%fh+fh)%fh-fh/2);
            int xx=abs(((x-ix+fw/2)%fw+fw)%fw-fw/2);
            if(ncc[(size_t)y*fw+x]>-1.5f && !(yy<=6&&xx<=6)){
                double v=ncc[(size_t)y*fw+x];side_sum+=v;side_sq+=v*v;side_n++;
            }
        }
        if(side_n<50){snprintf(out->reason,sizeof(out->reason),"too few sidelobe samples");goto success;}
        {
            double mean=side_sum/(double)side_n;
            double sd=sqrt(fmax(side_sq/(double)side_n-mean*mean,0.0));
            psr=(peak-mean)/(sd+1e-12);
        }
        out->dr=(lx+sx)*L;out->da=(ly+sy)*L;
        out->psr=psr;out->second=second_peak(ncc,fh,fw,iy,ix,10);
        out->peak_ratio=out->second>0?peak/out->second:INFINITY;
        out->time_read_s=tread;
        out->texture=a;memset(&a,0,sizeof(a));
        if(peak<0.10||psr<8.0){
            snprintf(out->reason,sizeof(out->reason),"global peak untrusted (ncc=%.3f or psr=%.1f)",peak,psr);
            goto success;
        }
        out->ok=1;
    }
success:
    out->time_s=now_s()-t0;
    ok=0;
done:
    free(pa.v);free(pb.v);free(a.v);free(b.v);free(ma.v);free(mb.v);free(ncc);
    if(ok!=0)ms_global_free(out);
    return ok;
}

void ms_global_free(MsGlobal *g){free(g->texture.v);memset(g,0,sizeof(*g));}

static void prep_image(float *a,size_t n,int mode)
{
    size_t i;double mean=0.0;
    if(mode==1)return;
    if(mode==2){for(i=0;i<n;i++)mean+=a[i];mean/=n;mean+=1e-6;}
    for(i=0;i<n;i++){
        if(mode==0)a[i]=sqrtf(a[i]);
        else if(mode==2)a[i]=log1pf(a[i]/(float)mean);
    }
}

static void abs_complex(float *ri,size_t n,float *out)
{
    size_t i;for(i=0;i<n;i++)out[i]=hypotf(ri[2*i],ri[2*i+1]);
}

static int complex_upsample(const float *in,int h,int w,int osa,int osr,int crop,
                            int threads,float **out,int *oh,int *ow)
{
    int fh=h*osa,fw=w*osr,hh=h/2,ww=w/2,y,x;
    fftwf_complex *full=fftwf_alloc_complex((size_t)h*(size_t)w);
    fftwf_complex *pad=fftwf_alloc_complex((size_t)fh*(size_t)fw);
    fftwf_complex *zf=fftwf_alloc_complex((size_t)fh*(size_t)fw);
    fftwf_plan pf,pi;
    float *amp;
    setup_threads(threads);
    if(!full||!pad||!zf){fftwf_free(full);fftwf_free(pad);fftwf_free(zf);return -1;}
    /* FFTW complex-to-complex input uses interleaved real/imag values. */
    pf=fftwf_plan_dft_2d(h,w,(fftwf_complex *)in,full,FFTW_FORWARD,FFTW_ESTIMATE);
    if(!pf)goto fail;
    fftwf_execute(pf);fftwf_destroy_plan(pf);
    memset(pad,0,(size_t)fh*(size_t)fw*sizeof(*pad));
    for(y=0;y<h;y++)for(x=0;x<w;x++){
        int yy=y<hh?y:fh-(h-y);
        int xx=x<ww?x:fw-(w-x);
        pad[(size_t)yy*fw+xx][0]=full[(size_t)y*w+x][0];
        pad[(size_t)yy*fw+xx][1]=full[(size_t)y*w+x][1];
    }
    pi=fftwf_plan_dft_2d(fh,fw,pad,zf,FFTW_BACKWARD,FFTW_ESTIMATE);
    if(!pi)goto fail;
    fftwf_execute(pi);fftwf_destroy_plan(pi);
    *oh=(h-2*crop)*osa;*ow=(w-2*crop)*osr;
    amp=malloc((size_t)*oh*(size_t)*ow*sizeof(*amp));
    if(!amp)goto fail;
    {
        float scale=(float)(osa*osr)/((float)fh*(float)fw);
        for(y=0;y<*oh;y++)for(x=0;x<*ow;x++){
            int yy=(y+crop*osa),xx=(x+crop*osr);
            float zr=zf[(size_t)yy*fw+xx][0],zi=zf[(size_t)yy*fw+xx][1];
            amp[(size_t)y*(*ow)+x]=hypotf(zr*scale,zi*scale);
        }
    }
    *out=amp;fftwf_free(full);fftwf_free(pad);fftwf_free(zf);return 0;
fail:
    fftwf_free(full);fftwf_free(pad);fftwf_free(zf);return -1;
}

int ms_match_patch(const float *tm,const float *ss,int h,int w,int H,int W,
                   int threads,double *dy,double *dx,double *ncc_out,
                   double *psr_out,int *edge,MsError *err)
{
    int ly=H-h+1,lx=W-w+1,fh,fw,y,x,iy=0,ix=0;
    float mt=0,ms=0;
    double tt=0,peak=-INFINITY,side_sum=0,side_sq=0;
    size_t nside=0;
    float *tz=NULL,*ss0=NULL,*num=NULL;
    double *den=NULL,*ncc=NULL,*var=NULL,*s1=NULL,*s2=NULL;
    fftwf_complex *ft=NULL,*fs=NULL,*g=NULL,*g_inverse=NULL;
    fftwf_plan p1=NULL,p2=NULL,pi=NULL;
    size_t i;
    if(ly<3||lx<3)return 1;
    fh=fast_len(H);fw=fast_len(W);setup_threads(threads);
    tz=malloc((size_t)h*(size_t)w*sizeof(*tz));
    ss0=malloc((size_t)H*(size_t)W*sizeof(*ss0));
    num=calloc((size_t)fh*(size_t)fw,sizeof(*num));
    den=malloc((size_t)ly*(size_t)lx*sizeof(*den));
    ncc=malloc((size_t)ly*(size_t)lx*sizeof(*ncc));
    var=malloc((size_t)ly*(size_t)lx*sizeof(*var));
    s1=calloc((size_t)(H+1)*(size_t)(W+1),sizeof(*s1));
    s2=calloc((size_t)(H+1)*(size_t)(W+1),sizeof(*s2));
    if(!tz||!ss0||!num||!den||!ncc||!var||!s1||!s2)goto oom;
    for(i=0;i<(size_t)h*w;i++)mt+=tm[i];
    mt/=(float)(h*w);
    for(i=0;i<(size_t)H*W;i++)ms+=ss[i];
    ms/=(float)(H*W);
    for(i=0;i<(size_t)h*w;i++){tz[i]=tm[i]-mt;tt+=(double)tz[i]*(double)tz[i];}
    if(tt<1e-9)goto invalid;
    for(i=0;i<(size_t)H*W;i++)ss0[i]=ss[i]-ms;
    ft=fftwf_alloc_complex(csize(fh,fw));fs=fftwf_alloc_complex(csize(fh,fw));g=fftwf_alloc_complex(csize(fh,fw));g_inverse=fftwf_alloc_complex(csize(fh,fw));
    if(!ft||!fs||!g||!g_inverse)goto oom;
    {
        float *pt=fftwf_alloc_real((size_t)fh*fw),*ps=fftwf_alloc_real((size_t)fh*fw);
        if(!pt||!ps){fftwf_free(pt);fftwf_free(ps);goto oom;}
        memset(pt,0,(size_t)fh*fw*sizeof(*pt));memset(ps,0,(size_t)fh*fw*sizeof(*ps));
        for(y=0;y<h;y++)memcpy(pt+(size_t)y*fw,tz+(size_t)y*w,(size_t)w*sizeof(*tz));
        for(y=0;y<H;y++)memcpy(ps+(size_t)y*fw,ss0+(size_t)y*W,(size_t)W*sizeof(*ss0));
        p1=fftwf_plan_dft_r2c_2d(fh,fw,pt,ft,FFTW_ESTIMATE);
        p2=fftwf_plan_dft_r2c_2d(fh,fw,ps,fs,FFTW_ESTIMATE);
        if(!p1||!p2){fftwf_free(pt);fftwf_free(ps);goto oom;}
        fftwf_execute(p1);fftwf_execute(p2);
        for(i=0;i<csize(fh,fw);i++){
            g[i][0]=ft[i][0]*fs[i][0]+ft[i][1]*fs[i][1];
            g[i][1]=ft[i][0]*fs[i][1]-ft[i][1]*fs[i][0];
        }
        memcpy(g_inverse,g,csize(fh,fw)*sizeof(*g_inverse));
        pi=fftwf_plan_dft_c2r_2d(fh,fw,g_inverse,num,FFTW_ESTIMATE);
        if(!pi){fftwf_free(pt);fftwf_free(ps);goto oom;}
        fftwf_execute(pi);
        for(i=0;i<(size_t)fh*fw;i++)num[i]/=(float)(fh*fw);
        fftwf_free(pt);fftwf_free(ps);
    }
    fftwf_destroy_plan(p1);p1=NULL;fftwf_destroy_plan(p2);p2=NULL;fftwf_destroy_plan(pi);pi=NULL;
    for(y=0;y<H;y++)for(x=0;x<W;x++){
        size_t k=(size_t)(y+1)*(W+1)+(x+1);
        double v=ss0[(size_t)y*W+x];
        s1[k]=v+s1[k-(W+1)]+s1[k-1]-s1[k-(W+2)];
        s2[k]=v*v+s2[k-(W+1)]+s2[k-1]-s2[k-(W+2)];
    }
    {
        double *positive=malloc((size_t)ly*lx*sizeof(*positive));size_t nv=0;
        double med,vfloor;
        if(!positive)goto oom;
        for(y=0;y<ly;y++)for(x=0;x<lx;x++){
            double a=s1[(size_t)y*(W+1)+x],b=s1[(size_t)(y+h)*(W+1)+x];
            double c=s1[(size_t)y*(W+1)+x+w],d=s1[(size_t)(y+h)*(W+1)+x+w];
            double su=d-b-c+a;
            a=s2[(size_t)y*(W+1)+x];b=s2[(size_t)(y+h)*(W+1)+x];
            c=s2[(size_t)y*(W+1)+x+w];d=s2[(size_t)(y+h)*(W+1)+x+w];
            var[(size_t)y*lx+x]=d-b-c+a-su*su/(double)(h*w);
            if(var[(size_t)y*lx+x]>0)positive[nv++]=var[(size_t)y*lx+x];
        }
        med=ms_median(positive,nv);free(positive);
        vfloor=fmax(1e-9,1e-4*med);
        for(y=0;y<ly;y++)for(x=0;x<lx;x++){
            size_t k=(size_t)y*lx+x;
            double vv=fmax(var[k],vfloor);
            den[k]=sqrt(tt*vv);
            ncc[k]=(double)num[(size_t)y*fw+x]/den[k];
            if(ncc[k]>peak){peak=ncc[k];iy=y;ix=x;}
        }
    }
    {
        unsigned char *excluded=calloc((size_t)ly*lx,1);
        if(!excluded)goto oom;
        for(y=iy-2;y<=iy+2;y++)for(x=ix-2;x<=ix+2;x++){
            if(y>=0&&y<ly&&x>=0&&x<lx)excluded[(size_t)y*lx+x]=1;
        }
        for(i=0;i<(size_t)ly*lx;i++)if(!excluded[i]){side_sum+=ncc[i];side_sq+=ncc[i]*ncc[i];nside++;}
        free(excluded);
    }
    *psr_out=nside>8?(peak-side_sum/nside)/(sqrt(fmax(side_sq/nside-(side_sum/nside)*(side_sum/nside),0.0))+1e-12):0.0;
    *edge=(iy==0||iy==ly-1||ix==0||ix==lx-1);
    *dy=iy-(ly-1)/2.0;*dx=ix-(lx-1)/2.0;*ncc_out=peak;
    if(*edge)goto done;
    {
        int up=16,ny=2*up+1,nx=2*up+1;
        int nh=fh,nxh=fw/2+1;
        double *eyr=malloc((size_t)ny*nh*sizeof(*eyr)),*eyi=malloc((size_t)ny*nh*sizeof(*eyi));
        double *exr=malloc((size_t)nx*nxh*sizeof(*exr)),*exi=malloc((size_t)nx*nxh*sizeof(*exi));
        double *cup=malloc((size_t)ny*nx*sizeof(*cup));
        double *tmp_r=malloc((size_t)ny*nxh*sizeof(*tmp_r));
        double *tmp_i=malloc((size_t)ny*nxh*sizeof(*tmp_i));
        if(!eyr||!eyi||!exr||!exi||!cup||!tmp_r||!tmp_i){
            free(eyr);free(eyi);free(exr);free(exi);free(cup);free(tmp_r);free(tmp_i);goto oom;
        }
        for(y=0;y<ny;y++)for(x=0;x<nh;x++){
            double u=iy+(y-up)/(double)up;
            double f=(x<=nh/2?x:x-nh)/(double)nh;
            double a=2*M_PI*u*f;eyr[(size_t)y*nh+x]=cos(a);eyi[(size_t)y*nh+x]=sin(a);
        }
        for(y=0;y<nx;y++)for(x=0;x<nxh;x++){
            double v=ix+(y-up)/(double)up;
            double f=(double)x/(double)fw;
            double a=2*M_PI*f*v;exr[(size_t)y*nxh+x]=cos(a);exi[(size_t)y*nxh+x]=sin(a);
        }
        for(y=0;y<ny;y++)for(x=0;x<nxh;x++){
            double re=0.0,im=0.0;
            double wx=(x==0||(fw%2==0&&x==fw/2))?1.0:2.0;
            int ky;
            for(ky=0;ky<nh;ky++){
                size_t k=(size_t)ky*nxh+x;
                double gr=g[k][0]*wx,gi=g[k][1]*wx;
                double er=eyr[(size_t)y*nh+ky],ei=eyi[(size_t)y*nh+ky];
                double tr=er*gr-ei*gi,ti=er*gi+ei*gr;
                re+=tr;im+=ti;
            }
            tmp_r[(size_t)y*nxh+x]=re;
            tmp_i[(size_t)y*nxh+x]=im;
        }
        for(y=0;y<ny;y++)for(x=0;x<nx;x++){
            double re=0.0;
            int kx;
            for(kx=0;kx<nxh;kx++){
                double tr=tmp_r[(size_t)y*nxh+kx],ti=tmp_i[(size_t)y*nxh+kx];
                double xr=exr[(size_t)x*nxh+kx],xi=exi[(size_t)x*nxh+kx];
                re+=tr*xr-ti*xi;
            }
            cup[(size_t)y*nx+x]=re/((double)fh*fw);
        }
        {
            double best=-INFINITY;int by=0,bx=0;
            for(y=0;y<ny;y++)for(x=0;x<nx;x++){
                double u=iy+(y-up)/(double)up,v=ix+(x-up)/(double)up;
                int u0=(int)floor(u),v0=(int)floor(v);
                double tu=u-u0,tv=v-v0,d;
                if(u0<0)u0=0;
                if(u0>ly-2)u0=ly-2;
                if(v0<0)v0=0;
                if(v0>lx-2)v0=lx-2;
                d=(1-tu)*(1-tv)*den[(size_t)u0*lx+v0]+(1-tu)*tv*den[(size_t)u0*lx+v0+1]+tu*(1-tv)*den[(size_t)(u0+1)*lx+v0]+tu*tv*den[(size_t)(u0+1)*lx+v0+1];
                cup[(size_t)y*nx+x]/=d;
                if(cup[(size_t)y*nx+x]>best){best=cup[(size_t)y*nx+x];by=y;bx=x;}
            }
            {
                double fu=(by>0&&by<ny-1)?parabola(cup[(size_t)(by-1)*nx+bx],cup[(size_t)by*nx+bx],cup[(size_t)(by+1)*nx+bx])/up:0;
                double fv=(bx>0&&bx<nx-1)?parabola(cup[(size_t)by*nx+bx-1],cup[(size_t)by*nx+bx],cup[(size_t)by*nx+bx+1])/up:0;
                *dy=iy+(by-up)/(double)up+fu-(ly-1)/2.0;
                *dx=ix+(bx-up)/(double)up+fv-(lx-1)/2.0;
                *ncc_out=best;
            }
        }
        free(eyr);free(eyi);free(exr);free(exi);free(cup);free(tmp_r);free(tmp_i);
    }
done:
    free(tz);free(ss0);free(num);free(den);free(ncc);free(var);free(s1);free(s2);
    fftwf_free(ft);fftwf_free(fs);fftwf_free(g);fftwf_free(g_inverse);
    return 0;
invalid:
    free(tz);free(ss0);free(num);free(den);free(ncc);free(var);free(s1);free(s2);
    fftwf_free(ft);fftwf_free(fs);fftwf_free(g);fftwf_free(g_inverse);
    return 1;
oom:
    snprintf(err->message,sizeof(err->message),"out of memory or FFTW plan failure");
    free(tz);free(ss0);free(num);free(den);free(ncc);free(var);free(s1);free(s2);
    fftwf_destroy_plan(p1);fftwf_destroy_plan(p2);fftwf_destroy_plan(pi);
    fftwf_free(ft);fftwf_free(fs);fftwf_free(g);fftwf_free(g_inverse);
    return -1;
}

int ms_measure_points(const MsSlc *master,const MsSlc *slave,const MsModel *model,
                      const MsXY *xy,size_t nxy,int W,int R,int prep,int osa,int osr,
                      int threads,MsPoint **out,size_t *nout,MsError *err)
{
    size_t i,n=0;
    MsPoint *rows=calloc(nxy,sizeof(*rows));
    int crop=(osa==1&&osr==1)?0:3;
    if(!rows&&nxy){snprintf(err->message,sizeof(err->message),"out of memory");return -1;}
    for(i=0;i<nxy;i++){
        double pdr,pda,dy,dx,ncc,psr;
        int ox,oy,mx0,my0,edge;
        int th=W+2*crop,sh=W+2*R+2*crop,at=th*osa,ar=th*osr,st=sh*osa,sr=sh*osr;
        size_t nt=(size_t)th*th,ns=(size_t)sh*sh,j,zeros;
        float *zt=malloc(nt*2*sizeof(*zt)),*zs=malloc(ns*2*sizeof(*zs));
        float *tm=NULL,*ss=NULL;
        int tmh,tmw,ssh,ssw;
        ms_predict(model,xy[i].x,xy[i].y,&pdr,&pda);
        ox=round_even(pdr);oy=round_even(pda);
        mx0=round_even(xy[i].x)-W/2;my0=round_even(xy[i].y)-W/2;
        if(!zt||!zs){free(zt);free(zs);free(rows);snprintf(err->message,sizeof(err->message),"out of memory");return -1;}
        ms_slc_read(master,my0-crop,mx0-crop,th,th,zt);
        zeros=0;for(j=0;j<nt;j++)if(zt[2*j]==0&&zt[2*j+1]==0)zeros++;
        if(zeros*5>nt){free(zt);free(zs);continue;}
        ms_slc_read(slave,my0+oy-R-crop,mx0+ox-R-crop,sh,sh,zs);
        zeros=0;for(j=0;j<ns;j++)if(zs[2*j]==0&&zs[2*j+1]==0)zeros++;
        if(zeros*5>ns){free(zt);free(zs);continue;}
        if(osa==1&&osr==1){
            tmh=th;tmw=th;ssh=sh;ssw=sh;
            tm=malloc((size_t)tmh*tmw*sizeof(*tm));ss=malloc((size_t)ssh*ssw*sizeof(*ss));
            if(!tm||!ss){free(zt);free(zs);free(tm);free(ss);free(rows);return -1;}
            abs_complex(zt,nt,tm);abs_complex(zs,ns,ss);
        }else{
            tmh=at;tmw=ar;ssh=st;ssw=sr;
            if(complex_upsample(zt,th,th,osa,osr,crop,threads,&tm,&tmh,&tmw)!=0||
               complex_upsample(zs,sh,sh,osa,osr,crop,threads,&ss,&ssh,&ssw)!=0){
                free(zt);free(zs);free(tm);free(ss);free(rows);snprintf(err->message,sizeof(err->message),"complex FFT oversampling failed");return -1;
            }
        }
        free(zt);free(zs);
        prep_image(tm,(size_t)tmh*tmw,prep);prep_image(ss,(size_t)ssh*ssw,prep);
        {
        int match_rc=ms_match_patch(tm,ss,tmh,tmw,ssh,ssw,threads,&dy,&dx,&ncc,&psr,&edge,err);
        if(match_rc<0){
            free(tm);free(ss);free(rows);return -1;
        }
        if(match_rc>0){free(tm);free(ss);continue;}
        }
        free(tm);free(ss);
        rows[n].x=mx0+(W-1)/2.0;rows[n].y=my0+(W-1)/2.0;
        rows[n].dr=ox+dx/(double)osr;rows[n].da=oy+dy/(double)osa;
        rows[n].ncc=ncc;rows[n].psr=psr;rows[n].edge=edge;rows[n].valid=0;rows[n].inlier=0;n++;
    }
    *out=rows;*nout=n;return 0;
}
