#include "msreg.h"

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void quote(FILE *f,const char *s)
{
    const unsigned char *p=(const unsigned char *)(s?s:"");
    fputc('"',f);
    for(;*p;p++){
        if(*p=='"'||*p=='\\'){fputc('\\',f);fputc(*p,f);}
        else if(*p=='\n')fputs("\\n",f);
        else if(*p=='\r')fputs("\\r",f);
        else if(*p=='\t')fputs("\\t",f);
        else if(*p<32)fprintf(f,"\\u%04x",*p);
        else fputc(*p,f);
    }
    fputc('"',f);
}

static void number(FILE *f,double x)
{
    if(!isfinite(x))fputs("null",f);
    else fprintf(f,"%.10g",x);
}

static void log_metric(FILE *f,double x)
{
    if(isfinite(x))fprintf(f,"%.3f",x);
    else fputc('-',f);
}

static void check_json(FILE *f,const MsCheck *c)
{
    fprintf(f,"{\n    \"n_points\": %d,\n    \"n_valid\": %d,\n    \"valid_frac\": ",c->n_points,c->n_valid);
    number(f,c->valid_frac);
    fprintf(f,",\n    \"n_within_0p3\": %d,\n    \"frac_within_0p3_all_candidates\": ",c->n_within_0p3);
    number(f,c->frac_within_0p3_all);
    fputs(",\n    \"within_0p3_definition\": \"hypot(residual_r, residual_a) <= 0.3 pixel; no centering\",\n    \"frac_within_0p3_definition\": \"n_within_0p3 / n_valid (legacy; requires >=4 valid points)\",\n    \"frac_within_0p3_all_candidates_definition\": \"n_within_0p3 / n_points; null if no candidates\",\n    \"metric_population\": \"valid independent check points (ncc >= ncc_min and edge < 0.5)\",\n    \"metric_n_points\": ",f);
    fprintf(f,"%d",c->n_valid);
    fputs(",\n    \"mad_definition\": \"median(abs(residual - median(residual)))\",\n    \"robust_sigma_definition\": \"1.4826 * MAD (median-centered)\",\n    \"rmse_definition\": \"sqrt(mean(residual**2)); no mean subtraction\",\n    \"metric_units\": \"pixel\",\n    \"metric_min_points\": 4,\n    \"mad_r\": ",f);
    number(f,c->mad_r);
    fputs(",\n    \"robust_sigma_r\": ",f);number(f,c->robust_sigma_r);
    fputs(",\n    \"rmse_r\": ",f);number(f,c->rmse_r);
    fputs(",\n    \"mad_a\": ",f);number(f,c->mad_a);
    fputs(",\n    \"robust_sigma_a\": ",f);number(f,c->robust_sigma_a);
    fputs(",\n    \"rmse_a\": ",f);number(f,c->rmse_a);
    fputs(",\n    \"rms_definition\": \"1.4826 * MAD (median-centered); legacy rms_r/rms_a are not RMSE\"",f);
    if(c->has_metrics){
        fputs(",\n    \"median_r\": ",f);number(f,c->median_r);
        fputs(",\n    \"median_a\": ",f);number(f,c->median_a);
        fputs(",\n    \"rms_r\": ",f);number(f,c->robust_sigma_r);
        fputs(",\n    \"rms_a\": ",f);number(f,c->robust_sigma_a);
        fputs(",\n    \"max_abs_r\": ",f);number(f,c->max_abs_r);
        fputs(",\n    \"max_abs_a\": ",f);number(f,c->max_abs_a);
        fputs(",\n    \"frac_within_0p3\": ",f);number(f,c->frac_within_0p3);
    }
    fputs("\n  }",f);
}

static void gate_json(FILE *f,int pass,const char reasons[][512],int nreasons)
{
    int i;
    fprintf(f,"{\n    \"passed\": %s,\n    \"reasons\": [",pass?"true":"false");
    for(i=0;i<nreasons;i++){if(i)fputs(", ",f);quote(f,reasons[i]);}
    fputs("],\n    \"thresholds\": {\"min_stage\": \"S256\", \"min_inliers\": 20, \"min_inlier_frac\": 0.30, \"max_rms\": 0.35, \"min_coverage\": 5, \"check_median\": 0.12, \"check_rms\": 0.35, \"check_valid_frac\": 0.30}\n  }",f);
}

static void check_json_compact(FILE *f,const MsCheck *c)
{
    fprintf(f,"{\"n_points\": %d, \"n_valid\": %d, \"valid_frac\": ",c->n_points,c->n_valid);
    number(f,c->valid_frac);
    fprintf(f,", \"n_within_0p3\": %d, \"frac_within_0p3_all_candidates\": ",c->n_within_0p3);
    number(f,c->frac_within_0p3_all);
    fputs(", \"within_0p3_definition\": \"hypot(residual_r, residual_a) <= 0.3 pixel; no centering\", \"frac_within_0p3_definition\": \"n_within_0p3 / n_valid (legacy; requires >=4 valid points)\", \"frac_within_0p3_all_candidates_definition\": \"n_within_0p3 / n_points; null if no candidates\", \"metric_population\": \"valid independent check points (ncc >= ncc_min and edge < 0.5)\", \"metric_n_points\": ",f);
    fprintf(f,"%d",c->n_valid);
    fputs(", \"mad_definition\": \"median(abs(residual - median(residual)))\", \"robust_sigma_definition\": \"1.4826 * MAD (median-centered)\", \"rmse_definition\": \"sqrt(mean(residual**2)); no mean subtraction\", \"metric_units\": \"pixel\", \"metric_min_points\": 4, \"mad_r\": ",f);
    number(f,c->mad_r);
    fputs(", \"robust_sigma_r\": ",f);number(f,c->robust_sigma_r);
    fputs(", \"rmse_r\": ",f);number(f,c->rmse_r);
    fputs(", \"mad_a\": ",f);number(f,c->mad_a);
    fputs(", \"robust_sigma_a\": ",f);number(f,c->robust_sigma_a);
    fputs(", \"rmse_a\": ",f);number(f,c->rmse_a);
    fputs(", \"rms_definition\": \"1.4826 * MAD (median-centered); legacy rms_r/rms_a are not RMSE\"",f);
    if(c->has_metrics){
        fputs(", \"median_r\": ",f);number(f,c->median_r);
        fputs(", \"median_a\": ",f);number(f,c->median_a);
        fputs(", \"rms_r\": ",f);number(f,c->robust_sigma_r);
        fputs(", \"rms_a\": ",f);number(f,c->robust_sigma_a);
        fputs(", \"max_abs_r\": ",f);number(f,c->max_abs_r);
        fputs(", \"max_abs_a\": ",f);number(f,c->max_abs_a);
        fputs(", \"frac_within_0p3\": ",f);number(f,c->frac_within_0p3);
    }
    fputc('}',f);
}

static void gate_json_compact(FILE *f,int pass,const char reasons[][512],int nreasons)
{
    int i;
    fprintf(f,"{\"passed\": %s, \"reasons\": [",pass?"true":"false");
    for(i=0;i<nreasons;i++){if(i)fputs(", ",f);quote(f,reasons[i]);}
    fputs("], \"thresholds\": {\"min_stage\": \"S256\", \"min_inliers\": 20, \"min_inlier_frac\": 0.30, \"max_rms\": 0.35, \"min_coverage\": 5, \"check_median\": 0.12, \"check_rms\": 0.35, \"check_valid_frac\": 0.30}}",f);
}

static void model_json(FILE *f,const MsModel *m)
{
    int i;
    fputs("{\"cr\": [",f);
    for(i=0;i<3;i++){if(i)fputs(", ",f);number(f,m->r.c[i]);}
    fputs("], \"ca\": [",f);
    for(i=0;i<3;i++){if(i)fputs(", ",f);number(f,m->a.c[i]);}
    fputs("]}",f);
}

int ms_write_estimate_report(const char *path,const MsEstimate *e,MsError *err)
{
    FILE *f=fopen(path,"w");
    int i;
    if(!f){snprintf(err->message,sizeof(err->message),"cannot write report %s: %s",path,strerror(errno));return -1;}
    fputs("{\n  \"module\": \"ms_register\",\n  \"version\": ",f);quote(f,MSREG_VERSION);
    fputs(",\n  \"master\": {\"prm\": ",f);quote(f,e->master.prm_path);
    fputs(", \"slc\": ",f);quote(f,e->master.path);
    fprintf(f,", \"nx\": %d, \"ny\": %d},\n  \"slave\": {\"prm\": ",e->master.nx,e->master.ny);
    quote(f,e->slave.prm_path);fputs(", \"slc\": ",f);quote(f,e->slave.path);
    fprintf(f,", \"nx\": %d, \"ny\": %d},\n  \"stages\": [",e->slave.nx,e->slave.ny);
    for(i=0;i<e->nstages;i++){
        const MsStage *s=&e->stages[i];
        if(i)fputc(',',f);
        fprintf(f,"\n    {\"name\": \"S%d\", \"window\": %d, \"search\": %d, \"n_points\": %d, \"n_measured\": %d, \"n_valid\": %d, \"n_inliers\": %d, \"inlier_frac\": ",
                s->window,s->window,s->search,s->n_points,s->n_measured,s->n_valid,s->n_inliers);
        number(f,s->inlier_frac);
        fputs(", \"rms_r\": ",f);number(f,s->rms_r);
        fputs(", \"rms_a\": ",f);number(f,s->rms_a);
        fprintf(f,", \"time_s\": ");number(f,s->time_s);
        fputs(", \"mad_definition\": \"median(abs(residual - median(residual)))\", \"robust_sigma_definition\": \"1.4826 * MAD (median-centered)\", \"rmse_definition\": \"sqrt(mean(residual**2)); no mean subtraction\", \"metric_units\": \"pixel\", \"metric_min_points\": 4, \"mad_r\": ",f);number(f,s->mad_r);
        fputs(", \"robust_sigma_r\": ",f);number(f,s->robust_sigma_r);
        fputs(", \"rmse_r\": ",f);number(f,s->rmse_r);
        fputs(", \"mad_a\": ",f);number(f,s->mad_a);
        fputs(", \"robust_sigma_a\": ",f);number(f,s->robust_sigma_a);
        fputs(", \"rmse_a\": ",f);number(f,s->rmse_a);
        fputs(", \"rms_definition\": \"sqrt(mean(residual**2)); no mean subtraction\", \"metric_population\": \"final least-squares input inliers; before returned mask recomputation\", \"metric_n_points\": ",f);fprintf(f,"%d",s->metric_n_points);
        if(s->has_fit){
            fputs(", \"sigma\": ",f);number(f,s->sigma);
            fputs(", \"t_in\": ",f);number(f,s->t_in);
            fputs(", \"sigma_definition\": \"1.4826 * median(abs(residual)); pooled r/a on inl0 before least-squares refinement; not median-centered MAD\"",f);
            fprintf(f,", \"accepted\": %s, \"center_dr\": ",s->accepted?"true":"false");number(f,s->center_dr);
            fputs(", \"center_da\": ",f);number(f,s->center_da);fputs(", \"model\": ",f);model_json(f,&s->model);
        }else{fputs(", \"accepted\": false, \"reason\": ",f);quote(f,s->reason);}
        fputc('}',f);
    }
    fputs("\n  ],\n  \"warnings\": [],\n  \"global\": {\"ok\": ",f);fputs(e->global.ok?"true":"false",f);
    fprintf(f,", \"looks\": %d, \"ncc\": ",e->global.looks);number(f,e->global.ncc);
    fputs(", \"time_s\": ",f);number(f,e->global.time_s);
    fputs(", \"dr\": ",f);number(f,e->global.dr);fputs(", \"da\": ",f);number(f,e->global.da);
    fputs(", \"psr\": ",f);number(f,e->global.psr);fputs(", \"second\": ",f);number(f,e->global.second);
    fputs(", \"peak_ratio\": ",f);number(f,e->global.peak_ratio);fputs(", \"time_read_s\": ",f);number(f,e->global.time_read_s);
    fputs("},\n  \"time_total_s\": ",f);number(f,e->time_total_s);
    fprintf(f,",\n  \"status\": \"%s\"",e->gate_passed?"passed":"failed");
    if(e->has_model){
        fputs(",\n  \"final_stage\": ",f);quote(f,e->final_stage);
        fputs(",\n  \"model\": ",f);model_json(f,&e->model);
        fprintf(f,",\n  \"prm\": {\"rshift\": %s, \"sub_int_r\": ",e->prm_values[0]);quote(f,e->prm_values[1]);
        fputs(", \"stretch_r\": ",f);quote(f,e->prm_values[2]);
        fputs(", \"a_stretch_r\": ",f);quote(f,e->prm_values[3]);
        fprintf(f,", \"ashift\": %s, \"sub_int_a\": ",e->prm_values[4]);quote(f,e->prm_values[5]);
        fputs(", \"stretch_a\": ",f);quote(f,e->prm_values[6]);
        fputs(", \"a_stretch_a\": ",f);quote(f,e->prm_values[7]);fputc('}',f);
        fputs(",\n  \"center_offset\": {\"range\": ",f);number(f,e->center_dr);fputs(", \"azimuth\": ",f);number(f,e->center_da);fputc('}',f);
    }
    if(e->has_check){
        fputs(",\n  \"check\": ",f);check_json(f,&e->check);
        fputs(",\n  \"check_time_s\": ",f);number(f,e->check_time_s);
        fprintf(f,",\n  \"coverage_cells\": %d",e->coverage_cells);
    }
    fputs(",\n  \"gate\": ",f);gate_json(f,e->gate_passed,(const char (*)[512])e->gate_reasons,e->ngate_reasons);
    if(e->output_prm[0]){fputs(",\n  \"output_prm\": ",f);quote(f,e->output_prm);}
    if(e->fail_reason[0]){fputs(",\n  \"fail_reason\": ",f);quote(f,e->fail_reason);}
    fputs("\n}\n",f);
    if(fclose(f)!=0){snprintf(err->message,sizeof(err->message),"cannot close report: %s",path);return -1;}
    return 0;
}

int ms_write_verify_report(const char *path,const MsSlc *m,const MsSlc *s,const MsVerify *v,MsError *err)
{
    FILE *f=fopen(path,"w");if(!f){snprintf(err->message,sizeof(err->message),"cannot write report: %s",path);return -1;}
    fputs("{\n  \"module\": \"ms_register\",\n  \"version\": ",f);quote(f,MSREG_VERSION);
    fputs(",\n  \"mode\": \"verify\",\n  \"master\": ",f);quote(f,m->path);fputs(",\n  \"slave\": ",f);quote(f,s->path);
    fputs(",\n  \"check\": ",f);check_json(f,&v->check);fputs(",\n  \"gate\": ",f);gate_json(f,v->passed,(const char (*)[512])v->reasons,v->nreasons);
    fputs(",\n  \"time_s\": ",f);number(f,v->time_s);fputs("\n}\n",f);
    if(fclose(f)!=0){snprintf(err->message,sizeof(err->message),"cannot close report: %s",path);return -1;}return 0;
}

int ms_write_log(const char *path,const MsEstimate *e,const char *master,const char *slave,MsError *err)
{
    FILE *f=fopen(path,"w");int i;
    if(!f){snprintf(err->message,sizeof(err->message),"cannot write log: %s",path);return -1;}
    fprintf(f,"# ms_register %s  多尺度配准结果\nmaster_prm = %s\nslave_prm = %s\nstatus = %s\n",MSREG_VERSION,master,slave,e->gate_passed?"passed":"failed");
    if(e->has_model){
        fprintf(f,"range_offset = %.3f pixel\nazimuth_offset = %.3f pixel\n",e->center_dr,e->center_da);
        fputs("# 偏移取影像中心处, 从影像坐标 = 主影像坐标 + 偏移\n",f);
        fprintf(f,"rshift = %s\nsub_int_r = %s\nstretch_r = %s\na_stretch_r = %s\nashift = %s\nsub_int_a = %s\nstretch_a = %s\na_stretch_a = %s\n",
                e->prm_values[0],e->prm_values[1],e->prm_values[2],e->prm_values[3],e->prm_values[4],e->prm_values[5],e->prm_values[6],e->prm_values[7]);
    }
    if(e->global.ok)fprintf(f,"global: dr=%.2f da=%.2f ncc=%.3f psr=%.1f\n",e->global.dr,e->global.da,e->global.ncc,e->global.psr);
    for(i=0;i<e->nstages;i++){
        const MsStage *s=&e->stages[i];
        fprintf(f,"stage S%d: points=%d valid=%d inliers=%d (%.0f%%) rmse=(",
                s->window,s->n_points,s->n_valid,s->n_inliers,100*s->inlier_frac);
        log_metric(f,s->rms_r);fputc(',',f);log_metric(f,s->rms_a);
        fputs(") robust_sigma=(",f);log_metric(f,s->robust_sigma_r);fputc(',',f);log_metric(f,s->robust_sigma_a);
        fprintf(f,") accepted=%s\n",s->accepted?"True":"False");
    }
    if(e->has_check){fputs("check = ",f);check_json_compact(f,&e->check);fputc('\n',f);}
    fputs("gate = ",f);gate_json_compact(f,e->gate_passed,(const char (*)[512])e->gate_reasons,e->ngate_reasons);fputc('\n',f);
    if(e->fail_reason[0])fprintf(f,"fail_reason = %s\n",e->fail_reason);
    if(fclose(f)!=0){snprintf(err->message,sizeof(err->message),"cannot close log: %s",path);return -1;}return 0;
}

int ms_write_points(const char *path,const MsEstimate *e,MsError *err)
{
    FILE *f=fopen(path,"w");size_t i;
    if(!f){snprintf(err->message,sizeof(err->message),"cannot write points: %s",path);return -1;}
    fprintf(f,"# x_master dr y_master da ncc*100 inlier  (window=%d)  从影像坐标 = 主影像坐标 + 偏移\n",e->final_window);
    for(i=0;i<e->nfinal_points;i++){
        const MsPoint *p=&e->final_points[i];
        fprintf(f," %.1f %.3f %.1f %.3f %.2f %d\n",p->x,p->dr,p->y,p->da,100*p->ncc,p->inlier?1:0);
    }
    if(fclose(f)!=0){snprintf(err->message,sizeof(err->message),"cannot close points: %s",path);return -1;}return 0;
}
