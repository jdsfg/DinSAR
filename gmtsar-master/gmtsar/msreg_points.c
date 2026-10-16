#define _POSIX_C_SOURCE 200809L
#include "msreg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int ms_make_lattice(const MsSlc *master,const MsSlc *slave,const MsModel *model,
                    int W,int R,int nxp,int nyp,const MsImage *texture,int looks,
                    double phase_x,double phase_y,MsXY **out,size_t *nout)
{
    int half=W/2+R+4;
    double span_x=(double)master->nx-2.0*half;
    double span_y=(double)master->ny-2.0*half;
    double cw,ch;
    int cx,cy;
    size_t n=0,cap=(size_t)nxp*(size_t)nyp;
    MsXY *pts=malloc(cap*sizeof(*pts));
    double positive_med=NAN;
    if(!pts)return -1;
    if(span_x<=0.0||span_y<=0.0){free(pts);*out=NULL;*nout=0;return 0;}
    cw=span_x/(double)nxp;
    ch=span_y/(double)nyp;
    if(texture&&texture->v){
        size_t nn=(size_t)texture->h*(size_t)texture->w;
        double *positive=malloc(nn*sizeof(*positive));size_t np=0,i;
        if(positive){
            for(i=0;i<nn;i++)if(texture->v[i]>0.0f)positive[np++]=texture->v[i];
            if(np)positive_med=ms_median(positive,np);
            free(positive);
        }
    }
    for(cy=0;cy<nyp;cy++)for(cx=0;cx<nxp;cx++){
        double sx=half+(cx+phase_x)*cw;
        double sy=half+(cy+phase_y)*ch;
        double dr=0.0,da=0.0;
        if(texture&&texture->v){
            double tx0=sx-.35*cw,tx1=sx+.35*cw;
            double ty0=sy-.35*ch,ty1=sy+.35*ch;
            int i0=(int)(tx0/looks),i1=(int)(tx1/looks)+1;
            int j0=(int)(ty0/looks),j1=(int)(ty1/looks)+1;
            int x,y,bestx=-1,besty=-1;
            float best=-INFINITY;
            if(i0<0)i0=0;
            if(j0<0)j0=0;
            if(i1>texture->w)i1=texture->w;
            if(j1>texture->h)j1=texture->h;
            for(y=j0;y<j1;y++)for(x=i0;x<i1;x++){
                float v=texture->v[(size_t)y*(size_t)texture->w+(size_t)x];
                if(v>best){best=v;bestx=x;besty=y;}
            }
            if(j1>j0&&i1>i0){
                if(isfinite(positive_med)&&best<0.2*positive_med)continue;
                sx=(bestx+0.5)*looks;
                sy=(besty+0.5)*looks;
            }
        }
        sx=fmin(fmax(sx,(double)half),(double)master->nx-half);
        sy=fmin(fmax(sy,(double)half),(double)master->ny-half);
        if(model){
            ms_predict(model,sx,sy,&dr,&da);
            if(!isfinite(dr)||!isfinite(da))continue;
            if(!(half<=sx+dr&&sx+dr<slave->nx-half&&
                 half<=sy+da&&sy+da<slave->ny-half))continue;
        }
        pts[n].x=sx;pts[n].y=sy;n++;
    }
    *out=pts;*nout=n;return 0;
}

static void calc_metrics(const double *r,const double *a,size_t n,MsCheck *c)
{
    size_t i;
    double *rad=malloc(n*sizeof(*rad)),*absr=malloc(n*sizeof(*absr)),*absa=malloc(n*sizeof(*absa));
    double maxr=0,maxa=0,median_r,median_a;
    if(!rad||!absr||!absa){free(rad);free(absr);free(absa);return;}
    median_r=ms_median(r,n);median_a=ms_median(a,n);
    for(i=0;i<n;i++){
        rad[i]=hypot(r[i],a[i]);absr[i]=fabs(r[i]);absa[i]=fabs(a[i]);
        if(fabs(r[i]-median_r)>maxr)maxr=fabs(r[i]-median_r);
        if(fabs(a[i]-median_a)>maxa)maxa=fabs(a[i]-median_a);
        if(rad[i]<=0.3)c->n_within_0p3++;
    }
    c->median_r=median_r;c->median_a=median_a;
    c->mad_r=ms_mad(r,n);c->mad_a=ms_mad(a,n);
    c->robust_sigma_r=1.4826*c->mad_r;c->robust_sigma_a=1.4826*c->mad_a;
    c->rmse_r=sqrt(0.0);c->rmse_a=sqrt(0.0);
    for(i=0;i<n;i++){c->rmse_r+=r[i]*r[i];c->rmse_a+=a[i]*a[i];}
    c->rmse_r=sqrt(c->rmse_r/(double)n);c->rmse_a=sqrt(c->rmse_a/(double)n);
    c->max_abs_r=maxr;c->max_abs_a=maxa;c->frac_within_0p3=(double)c->n_within_0p3/(double)n;
    c->has_metrics=1;free(rad);free(absr);free(absa);
}

int ms_check_model(const MsSlc *master,const MsSlc *slave,const MsModel *model,
                   int grid_x,int grid_y,int prep,int ovs_a,int ovs_r,int threads,
                   int use_texture,MsCheck *out,MsError *err)
{
    MsGlobal global={0};MsXY *xy=NULL;size_t nxy=0;
    MsPoint *points=NULL;size_t nmeas=0,i,nvalid=0;
    double *rr=NULL,*ra=NULL,t0=now_s();
    MsModel m=*model;
    int L=1,target=1500;
    MsImage texture={0};
    memset(out,0,sizeof(*out));
    if(use_texture){
        if(ms_global_stage(master,slave,target,threads,&global,err)!=0)return -1;
        L=global.looks;
        texture.h=global.texture.h;texture.w=global.texture.w;
        texture.v=malloc((size_t)texture.h*(size_t)texture.w*sizeof(*texture.v));
        if(!texture.v){ms_global_free(&global);return -1;}
        {
            size_t k,n=(size_t)texture.h*(size_t)texture.w;
            for(k=0;k<n;k++)texture.v[k]=global.texture.v[k]*global.texture.v[k];
            {
                int box=(int)nearbyint(128.0/L);if(box<3)box=3;
                ms_box_mean(texture.v,texture.h,texture.w,box,global.texture.v);
            }
        }
    } else {
        global.texture.h=0;
        global.texture.w=0;
        global.texture.v=NULL;
    }
    if(ms_make_lattice(master,slave,&m,128,8,grid_x,grid_y,
                       use_texture?&global.texture:NULL,L,.25,.75,&xy,&nxy)!=0)goto fail;
    if(ms_measure_points(master,slave,&m,xy,nxy,128,8,prep,ovs_a,ovs_r,threads,&points,&nmeas,err)!=0)goto fail;
    rr=malloc(nmeas*sizeof(*rr));ra=malloc(nmeas*sizeof(*ra));
    if(!rr||!ra)goto fail;
    out->n_points=(int)nxy;
    for(i=0;i<nmeas;i++){
        double pr,pa;
        if(points[i].ncc<0.08||points[i].edge>=0.5)continue;
        ms_predict(&m,points[i].x,points[i].y,&pr,&pa);
        rr[nvalid]=points[i].dr-pr;ra[nvalid]=points[i].da-pa;nvalid++;
    }
    out->n_valid=(int)nvalid;
    out->valid_frac=nxy?(double)nvalid/(double)nxy:0.0;
    if(nvalid>=4)calc_metrics(rr,ra,nvalid,out);
    if(nxy){out->frac_within_0p3_all=(double)out->n_within_0p3/(double)nxy;}
    out->time_s=now_s()-t0;
    free(texture.v);free(xy);free(points);free(rr);free(ra);ms_global_free(&global);return 0;
fail:
    free(texture.v);free(xy);free(points);free(rr);free(ra);ms_global_free(&global);
    return -1;
}
