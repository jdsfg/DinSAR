#define _POSIX_C_SOURCE 200809L
#include "msreg.h"

#include <getopt.h>
#include <math.h>
#include <stdarg.h>
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
void ms_log(const MsOptions *opt,const char *fmt,...)
{
    va_list ap;
    if(opt->quiet)return;
    va_start(ap,fmt);
    vfprintf(stderr,fmt,ap);
    va_end(ap);
}
static void debug_stage(int window,const MsPoint *points,size_t n,const MsModel *model)
{
    const char *dir=getenv("MSREG_DEBUG_DIR");
    char path[2048];
    FILE *f;
    size_t i;
    int written;
    if(!dir||!*dir)return;
    written=snprintf(path,sizeof(path),"%s/stage_S%d.tsv",dir,window);
    if(written<0||(size_t)written>=sizeof(path))return;
    f=fopen(path,"w");
    if(!f){
        fprintf(stderr,"[ms_register] cannot write debug dump %s\n",path);
        return;
    }
    fputs("index\tx\ty\tdr\tda\tncc\tpsr\tedge\tinlier\tresidual_r\tresidual_a\n",f);
    for(i=0;i<n;i++){
        double pr,pa;
        ms_predict(model,points[i].x,points[i].y,&pr,&pa);
        fprintf(f,"%zu\t%.6f\t%.6f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%d\t%.12g\t%.12g\n",
        i,points[i].x,points[i].y,points[i].dr,points[i].da,points[i].ncc,
        points[i].psr,points[i].edge,points[i].inlier,
        points[i].dr-pr,points[i].da-pa);
    }
    fclose(f);
}
static void default_schedule(MsOptions *o)
{
    static const MsSchedule d[4]={
        {
            512,64,7,7,2.5,10,.20,1,1
        },
        {
            256,16,12,12,1.2,15,.25,1,1
        },
        {
            128,8,20,20,.7,25,.25,2,2
        },
        {
            64,6,28,28,.5,40,.25,2,2
        }
    };
    memcpy(o->schedule,d,sizeof(d));
    o->nschedule=4;
}
static int count_coverage(const MsPoint *p,size_t n,int nx,int ny)
{
    int cells[9]={
        0
    },i,count=0;
    for(i=0;i<(int)n;i++)if(p[i].inlier){
        int x=(int)(3.0*p[i].x/(double)nx),y=(int)(3.0*p[i].y/(double)ny);
        if(x<0)x=0;
        if(x>2)x=2;
        if(y<0)y=0;
        if(y>2)y=2;
        cells[y*3+x]=1;
    }
    for(i=0;i<9;i++)count+=cells[i];
    return count;
}
static void add_reason(MsEstimate *e,const char *reason)
{
    if(e->ngate_reasons<8)snprintf(e->gate_reasons[e->ngate_reasons++],512,"%s",reason);
}
int ms_estimate(const MsOptions *opt,MsEstimate *out,MsError *err)
{
    MsRng rng;
    MsImage power={
        0
    };
    double t0=now_s();
    int i,L,target=opt->global_target>0?opt->global_target:1500;
    MsModel current={
        {
            {
                0,0,0
            }
        },{
            {
                0,0,0
            }
        }
    };
    memset(out,0,sizeof(*out));
    out->master.fd=out->slave.fd=-1;
    if(ms_slc_open(&out->master,opt->master,NULL,err)!=0)return -1;
    if(ms_slc_open(&out->slave,opt->slave,NULL,err)!=0){
        ms_slc_close(&out->master);
        return -1;
    }
    if(ms_global_stage(&out->master,&out->slave,target,opt->threads,&out->global,err)!=0)goto fail;
    if(!out->global.ok){
        out->status=MSREG_EXIT_GATE;
        snprintf(out->fail_reason,sizeof(out->fail_reason),"%s",out->global.reason);
        out->time_total_s=now_s()-t0;
        return 0;
    }
    L=out->global.looks;
    current.r.c[0]=out->global.dr;
    current.a.c[0]=out->global.da;
    power.h=out->global.texture.h;
    power.w=out->global.texture.w;
    power.v=malloc((size_t)power.h*(size_t)power.w*sizeof(*power.v));
    if(!power.v){
        snprintf(err->message,sizeof(err->message),"out of memory");
        goto fail;
    }
    for(i=0;i<power.h*power.w;i++)power.v[i]=out->global.texture.v[i]*out->global.texture.v[i];
    ms_rng_seed(&rng,(uint64_t)(opt->seed<0?0:opt->seed));
    for(i=0;i<opt->nschedule;i++){
        MsSchedule sc=opt->schedule[i];
        MsImage texture={
            0
        };
        MsXY *xy=NULL;
        size_t nxy=0;
        MsPoint *measured=NULL;
        size_t nmeasured=0;
        MsPoint *valid=NULL;
        size_t nvalid=0,j;
        unsigned char *mask=NULL;
        MsStage *st=&out->stages[out->nstages];
        double fit_start=now_s();
        int box=(int)nearbyint((double)sc.W/(double)L);
        if(box<3)box=3;
        texture.h=power.h;
        texture.w=power.w;
        texture.v=malloc((size_t)power.h*(size_t)power.w*sizeof(*texture.v));
        if(!texture.v){
            snprintf(err->message,sizeof(err->message),"out of memory");
            free(xy);
            goto fail;
        }
        ms_box_mean(power.v,power.h,power.w,box,texture.v);
        if(ms_make_lattice(&out->master,&out->slave,&current,sc.W,sc.R,sc.nxp,sc.nyp,
        &texture,L,.5,.5,&xy,&nxy)!=0){
            free(texture.v);
            snprintf(err->message,sizeof(err->message),"out of memory");
            goto fail;
        }
        if(ms_measure_points(&out->master,&out->slave,&current,xy,nxy,sc.W,sc.R,
        opt->prep,sc.ovs_a,sc.ovs_r,opt->threads,
        &measured,&nmeasured,err)!=0){
            free(texture.v);
            free(xy);
            goto fail;
        }
        valid=malloc(nmeasured*sizeof(*valid));
        mask=calloc(nmeasured,sizeof(*mask));
        if((nmeasured&&(!valid||!mask))){
            free(texture.v);
            free(xy);
            free(measured);
            free(valid);
            free(mask);
            snprintf(err->message,sizeof(err->message),"out of memory");
            goto fail;
        }
        for(j=0;j<nmeasured;j++)if(measured[j].ncc>=.08&&measured[j].edge<.5)valid[nvalid++]=measured[j];
        memset(st,0,sizeof(*st));
        st->window=sc.W;
        st->search=sc.R;
        st->n_points=(int)nxy;
        st->n_measured=(int)nmeasured;
        st->n_valid=(int)nvalid;
        st->time_s=now_s()-fit_start;
        st->has_fit=0;
        if(nvalid>=6&&ms_fit_affine(valid,nvalid,sc.thr,&rng,&current,mask,st)==0){
            size_t in=0;
            for(j=0;j<nvalid;j++)if(mask[j])in++;
            st->has_fit=1;
            st->n_inliers=(int)in;
            st->inlier_frac=nvalid?(double)in/(double)nvalid:0.0;
            st->accepted=(in>=(size_t)sc.min_inl&&st->inlier_frac>=sc.min_frac);
            st->model=current;
            ms_predict(&current,out->master.nx/2.0,out->master.ny/2.0,&st->center_dr,&st->center_da);
            for(j=0;j<nmeasured;j++)measured[j].inlier=0;
            {
                size_t k=0;
                for(j=0;j<nmeasured;j++)if(measured[j].ncc>=.08&&measured[j].edge<.5){
                    measured[j].inlier=mask[k];
                    k++;
                }
            }
            if(!st->accepted)snprintf(st->reason,sizeof(st->reason),"stage gate rejected");
        }
        else{
            snprintf(st->reason,sizeof(st->reason),"not enough valid points or affine fit failed");
        }
        if(st->has_fit){
            for(j=0;j<nvalid;j++)valid[j].inlier=mask[j];
            debug_stage(sc.W,valid,nvalid,&current);
        }
        st->time_s=now_s()-fit_start;
        if(st->has_fit&&st->accepted){
            MsPoint *keep=malloc(nvalid*sizeof(*keep));
            if(nvalid&&!keep){
                free(texture.v);
                free(xy);
                free(measured);
                free(valid);
                free(mask);
                snprintf(err->message,sizeof(err->message),"out of memory");
                goto fail;
            }
            for(j=0;j<nvalid;j++)keep[j]=valid[j];
            free(out->final_points);
            out->final_points=keep;
            out->nfinal_points=nvalid;
            out->final_window=sc.W;
            out->model=current;
            out->has_model=1;
            snprintf(out->final_stage,sizeof(out->final_stage),"S%d",sc.W);
            out->center_dr=st->center_dr;
            out->center_da=st->center_da;
        }
        out->nstages++;
        free(texture.v);
        free(xy);
        free(measured);
        free(valid);
        free(mask);
        if(!st->accepted)break;
    }
    free(power.v);
    power.v=NULL;
    out->time_total_s=now_s()-t0;
    if(!out->has_model){
        out->status=MSREG_EXIT_GATE;
        snprintf(out->fail_reason,sizeof(out->fail_reason),"no stage produced a valid affine model");
        add_reason(out,out->fail_reason);
        goto gate_done;
    }
    out->check_time_s=now_s();
    if(ms_check_model(&out->master,&out->slave,&out->model,14,14,opt->prep,
    opt->has_ovs?opt->ovs_a:2,opt->has_ovs?opt->ovs_r:2,
    opt->threads,0,&out->check,err)!=0)goto fail;
    out->check_time_s=now_s()-out->check_time_s;
    out->has_check=1;
    out->coverage_cells=count_coverage(out->final_points,out->nfinal_points,
    out->master.nx,out->master.ny);
    {
        int fi=out->nstages-1;
        while(fi>=0&&out->stages[fi].window!=out->final_window)fi--;
        if(fi<0)fi=out->nstages-1;
        {
            MsStage *last=&out->stages[fi];
            if(out->final_window>256)add_reason(out,"minimum_stage: expected S256 or finer");
            if(last->n_inliers<20)add_reason(out,"final_inliers below 20");
            if(last->inlier_frac<.30)add_reason(out,"final inlier fraction below 0.30");
            if(last->rms_r>.35||last->rms_a>.35)add_reason(out,"final fit RMSE exceeds 0.35 px");
            if(out->coverage_cells<5)add_reason(out,"inliers cover fewer than 5 of 9 image cells");
            if(out->check.n_valid<8)add_reason(out,"check valid points below 8");
            if(out->check.valid_frac<.30)add_reason(out,"check valid fraction below 0.30");
            if(fmax(fabs(out->check.median_r),fabs(out->check.median_a))>.12)add_reason(out,"check median residual exceeds 0.12 px");
            if(out->check.robust_sigma_r>.35||out->check.robust_sigma_a>.35)add_reason(out,"check robust sigma exceeds 0.35 px");
        }
    }
    out->gate_passed=out->ngate_reasons==0;
    out->status=out->gate_passed?MSREG_EXIT_OK:MSREG_EXIT_GATE;
    gate_done:
    out->time_total_s=now_s()-t0;
    if(out->has_model)ms_prm_values(&out->model,out->prm_values);
    if(opt->output&&(out->gate_passed||opt->write_on_fail)&&out->has_model){
        const char *override=NULL;
        char *odir=strdup(opt->output),*sdir=strdup(out->slave.prm_path);
        if(odir&&sdir){
            char *a=strrchr(odir,'/'),*b=strrchr(sdir,'/');
            if(a&&b){
                *a='\0';
                *b='\0';
                if(strcmp(odir,sdir)!=0)override=out->slave.path;
            }
            else if((a!=NULL)!=(b!=NULL))override=out->slave.path;
        }
        if(ms_prm_write(&out->slave.prm,opt->output,&out->model,override,err)!=0){
            free(odir);
            free(sdir);
            goto fail;
        }
        free(odir);
        free(sdir);
        snprintf(out->output_prm,sizeof(out->output_prm),"%s",opt->output);
    }
    return 0;
    fail:
    free(power.v);
    ms_slc_close(&out->master);
    ms_slc_close(&out->slave);
    return -1;
}
int ms_verify_registered(const MsSlc *master,const MsSlc *slave,const MsOptions *opt,
MsVerify *out,MsError *err)
{
    MsModel zero={
        {
            {
                0,0,0
            }
        },{
            {
                0,0,0
            }
        }
    };
    double t=now_s();
    memset(out,0,sizeof(*out));
    if(master->nx!=slave->nx||master->ny!=slave->ny){
        snprintf(err->message,sizeof(err->message),"master/slave dimensions differ");
        return -1;
    }
    if(ms_check_model(master,slave,&zero,opt->grid_x,opt->grid_y,opt->prep,
    opt->has_ovs?opt->ovs_a:2,opt->has_ovs?opt->ovs_r:2,
    opt->threads,1,&out->check,err)!=0)return -1;
    if(out->check.n_valid<8)snprintf(out->reasons[out->nreasons++],512,"check valid points below 8");
    if(out->check.valid_frac<.30)snprintf(out->reasons[out->nreasons++],512,"check valid fraction below 0.30");
    if(fmax(fabs(out->check.median_r),fabs(out->check.median_a))>.12)snprintf(out->reasons[out->nreasons++],512,"check median residual exceeds 0.12 px");
    if(out->check.robust_sigma_r>.35||out->check.robust_sigma_a>.35)snprintf(out->reasons[out->nreasons++],512,"check robust sigma exceeds 0.35 px");
    out->passed=out->nreasons==0;
    out->time_s=now_s()-t;
    return 0;
}
static int parse_ovs(const char *s,int *a,int *r)
{
    char *end;
    long x=strtol(s,&end,10);
    if(*end!='x'&&*end!='X')return -1;
    long y=strtol(end+1,&end,10);
    if(*end||x<1||y<1||x>16||y>16)return -1;
    *a=(int)x;
    *r=(int)y;
    return 0;
}
static int parse_scales(const char *s,MsOptions *o)
{
    char *copy=strdup(s),*save=NULL,*tok;
    MsSchedule base[4];
    int n=0,i;
    if(!copy)return -1;
    default_schedule(o);
    memcpy(base,o->schedule,sizeof(base));
    o->nschedule=0;
    for(tok=strtok_r(copy,",",&save);tok;tok=strtok_r(NULL,",",&save)){
        int w=atoi(tok),found=0;
        if(n>=MSREG_MAX_STAGES){
            free(copy);
            return -1;
        }
        for(i=0;i<4;i++)if(base[i].W==w){
            MsSchedule sc=base[i];
            if(o->has_ovs){
                sc.ovs_a=o->ovs_a;
                sc.ovs_r=o->ovs_r;
            }
            o->schedule[n++]=sc;
            found=1;
            break;
        }
        if(!found||n>MSREG_MAX_STAGES){
            free(copy);
            return -1;
        }
    }
    free(copy);
    o->nschedule=n;
    return n?0:-1;
}
static int parse_cli(int argc,char **argv,MsOptions *o)
{
    memset(o,0,sizeof(*o));
    o->threads=1;
    o->global_target=1500;
    o->seed=0;
    o->prep=0;
    default_schedule(o);
    if(argc<2)return -1;
    if(strcmp(argv[1],"--version")==0){
        puts(MSREG_VERSION);
        exit(0);
    }
    if(strcmp(argv[1],"estimate")==0){
        static const struct option opts[]={
            {
                "output",required_argument,0,'o'
            },{
                "report",required_argument,0,1
            },
            {
                "log",required_argument,0,2
            },{
                "points",required_argument,0,3
            },
            {
                "threads",required_argument,0,4
            },{
                "seed",required_argument,0,5
            },
            {
                "global-target",required_argument,0,6
            },{
                "prep",required_argument,0,7
            },
            {
                "ovs",required_argument,0,8
            },{
                "write-on-fail",no_argument,0,9
            },
            {
                "scales",required_argument,0,10
            },
            {
                "quiet",no_argument,0,'q'
            },{
                0,0,0,0
            }
        };
        int c;
        optind=2;
        while((c=getopt_long(argc,argv,"o:q",opts,NULL))!=-1){
            switch(c){
                case 'o':o->output=optarg;
                break;
                case 'q':o->quiet=1;
                break;
                case 1:o->report=optarg;
                break;
                case 2:o->log=optarg;
                break;
                case 3:o->points=optarg;
                break;
                case 4:o->threads=atoi(optarg);
                break;
                case 5:o->seed=atoi(optarg);
                break;
                case 6:o->global_target=atoi(optarg);
                break;
                case 7:
                if(!strcmp(optarg,"sqrt"))o->prep=0;
                else if(!strcmp(optarg,"amp"))o->prep=1;
                else if(!strcmp(optarg,"log"))o->prep=2;
                else return -1;
                break;
                case 8:o->has_ovs=1;
                if(parse_ovs(optarg,&o->ovs_a,&o->ovs_r))return -1;
                break;
                case 9:o->write_on_fail=1;
                break;
                case 10:o->scales=optarg;
                break;
                default:return -1;
            }
        }
        if(o->scales&&parse_scales(o->scales,o)!=0)return -1;
        if(o->has_ovs){
            int j;
            for(j=0;j<o->nschedule;j++){
                o->schedule[j].ovs_a=o->ovs_a;
                o->schedule[j].ovs_r=o->ovs_r;
            }
        }
        if(argc-optind!=2)return -1;
        o->master=argv[optind];
        o->slave=argv[optind+1];
    }
    else if(strcmp(argv[1],"verify")==0){
        static const struct option opts[]={
            {
                "slave-slc",required_argument,0,1
            },{
                "report",required_argument,0,2
            },
            {
                "grid",required_argument,0,3
            },{
                "prep",required_argument,0,4
            },
            {
                "threads",required_argument,0,5
            },{
                "quiet",no_argument,0,'q'
            },{
                0,0,0,0
            }
        };
        int c;
        optind=2;
        o->grid_x=14;
        o->grid_y=14;
        while((c=getopt_long(argc,argv,"q",opts,NULL))!=-1)switch(c){
            case 'q':o->quiet=1;
            break;
            case 1:o->slave_slc=optarg;
            break;
            case 2:o->report=optarg;
            break;
            case 3:if(sscanf(optarg,"%dx%d",&o->grid_x,&o->grid_y)!=2)return -1;
            o->grid_is_set=1;
            break;
            case 4:if(!strcmp(optarg,"sqrt"))o->prep=0;
            else if(!strcmp(optarg,"amp"))o->prep=1;
            else if(!strcmp(optarg,"log"))o->prep=2;
            else return -1;
            break;
            case 5:o->threads=atoi(optarg);
            break;
            default:return -1;
        }
        if(argc-optind!=2)return -1;
        o->master=argv[optind];
        o->slave=argv[optind+1];
    }
    else return -1;
    if(o->threads<1)o->threads=1;
    return 0;
}
int ms_estimate_cli(const MsOptions *opt)
{
    MsEstimate e;
    MsError err={
        {
            0
        }
    };
    int rc;
    MsOptions run=*opt;
    if(run.output==NULL)run.output=run.slave;
    rc=ms_estimate(&run,&e,&err);
    if(rc!=0){
        fprintf(stderr,"ms_register: %s\n",err.message);
        return MSREG_EXIT_INPUT;
    }
    if(opt->report&&ms_write_estimate_report(opt->report,&e,&err)!=0){
        fprintf(stderr,"%s\n",err.message);
        rc=MSREG_EXIT_INPUT;
    }
    if(opt->log&&ms_write_log(opt->log,&e,opt->master,opt->slave,&err)!=0){
        fprintf(stderr,"%s\n",err.message);
        rc=MSREG_EXIT_INPUT;
    }
    if(opt->points&&ms_write_points(opt->points,&e,&err)!=0){
        fprintf(stderr,"%s\n",err.message);
        rc=MSREG_EXIT_INPUT;
    }
    if(!opt->quiet)fprintf(stderr,"ms_register estimate: %s\n",e.gate_passed?"passed":"failed");
    ms_slc_close(&e.master);
    ms_slc_close(&e.slave);
    free(e.final_points);
    return rc?rc:e.status;
}
int ms_verify_cli(const MsOptions *opt)
{
    MsSlc m,s;
    MsVerify v;
    MsError err={
        {
            0
        }
    };
    if(ms_slc_open(&m,opt->master,NULL,&err)!=0){
        fprintf(stderr,"%s\n",err.message);
        return MSREG_EXIT_INPUT;
    }
    if(ms_slc_open(&s,opt->slave,opt->slave_slc,&err)!=0){
        fprintf(stderr,"%s\n",err.message);
        ms_slc_close(&m);
        return MSREG_EXIT_INPUT;
    }
    if(ms_verify_registered(&m,&s,opt,&v,&err)!=0){
        fprintf(stderr,"%s\n",err.message);
        ms_slc_close(&m);
        ms_slc_close(&s);
        return MSREG_EXIT_INPUT;
    }
    if(opt->report&&ms_write_verify_report(opt->report,&m,&s,&v,&err)!=0){
        fprintf(stderr,"%s\n",err.message);
        ms_slc_close(&m);
        ms_slc_close(&s);
        return MSREG_EXIT_INPUT;
    }
    if(!opt->quiet)fprintf(stderr,"ms_register verify: %s\n",v.passed?"passed":"failed");
    ms_slc_close(&m);
    ms_slc_close(&s);
    return v.passed?0:MSREG_EXIT_GATE;
}
int main(int argc,char **argv)
{
    MsOptions o;
    if(argc==4&&strcmp(argv[1],"--rng-dump")==0)
    return ms_rng_dump_test((uint64_t)strtoull(argv[2],NULL,10),
    strtoul(argv[3],NULL,10));
    if(argc==5&&strcmp(argv[1],"--rng-pattern")==0)
    return ms_rng_pattern_test((uint64_t)strtoull(argv[2],NULL,10),
    (uint32_t)strtoul(argv[3],NULL,10),
    (uint32_t)strtoul(argv[4],NULL,10));
    if(argc==10&&strcmp(argv[1],"--prm-roundtrip")==0){
        MsPrm p;
        MsError err={
            {
                0
            }
        };
        MsModel model={
            {
                {
                    0,0,0
                }
            },{
                {
                    0,0,0
                }
            }
        };
        int i;
        if(ms_prm_read(argv[2],&p,&err)!=0){
            fprintf(stderr,"%s\n",err.message);
            return MSREG_EXIT_INPUT;
        }
        for(i=0;i<3;i++)model.r.c[i]=strtod(argv[4+i],NULL);
        for(i=0;i<3;i++)model.a.c[i]=strtod(argv[7+i],NULL);
        if(ms_prm_write(&p,argv[3],&model,NULL,&err)!=0){
            fprintf(stderr,"%s\n",err.message);
            ms_prm_free(&p);
            return MSREG_EXIT_INPUT;
        }
        ms_prm_free(&p);
        return MSREG_EXIT_OK;
    }
    if(parse_cli(argc,argv,&o)!=0){
        fprintf(stderr,"usage: ms_register estimate [options] master.prm slave.prm | ms_register verify [options] master.prm slave.prm\n");
        return MSREG_EXIT_INPUT;
    }
    return strcmp(argv[1],"estimate")==0?ms_estimate_cli(&o):ms_verify_cli(&o);
}
