// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include <signal.h>
#include "common.h"
#include "meschach/matrix.h"

Real iota1=5.2;
Real iota2=0.619;
Real phi=17;
Real eta1=1.703205e-5;
Real eta2=2.9526;
int interactive_mode_requested = 0;

void request_interactive_mode(int s) {
  if(interactive_mode_requested == 0) {
    printf("\nInteractive mode requested. Please wait.\n");
    interactive_mode_requested = 1;
    signal(SIGINT, request_interactive_mode);
  }
  else {
    // User is probably trying to exit the app - they've hit ctrl+c twice
    exit(0);
  }
}

void spade_v_output(VEC* vec) {
  printf("\n\n");
  for(int i = 0; i < vec->dim; i++) {
    printf("%d:\t%Lf\n", i, vec->ve[i]);
  }
}

void data_read_ce(char * data_file_name, Data * data, int * N, Real k) {
  char buffer[30];
  sprintf(buffer,"%s-ce.dat",data_file_name);

  int Nce;
  FILE *fp = fopen(buffer,"r");
  fscanf(fp,"%d",&Nce);
  Real *ct = (Real *) calloc(Nce,sizeof(Real));
  Real *ti = (Real *) calloc(Nce,sizeof(Real));

  for (int i=0;i<Nce;i++) {
  #if REAL == DOUBLE
    // lf
    fscanf(fp,"%lf %lf", &ct[i],&ti[i]);
  #elif REAL == FLOAT
    // f
    fscanf(fp,"%f %f", &ct[i],&ti[i]);
  #elif REAL == LONGDOUBLE
    // Lf
    fscanf(fp,"%Lf %Lf", &ct[i],&ti[i]);
  #endif
  }

  VEC *vti = v_get(Nce);
  for (int i=0;i<Nce;i++) {
    vti->ve[i] = ti[i];
  }

  PERM *order = px_get(vti->dim);
  v_sort(vti,order);

  // Data file year is zero based but does not necessarily start at zero
  // we want to take floor of vti->ve[0] and subtract it from the ceil of vti->ve[Nce-1]
  // to ensure we are getting the correct Y (total number of years).
  data->Y = ceil(vti->ve[Nce-1]) - floor(vti->ve[0]);
  *N = round(data->Y/k);

  PX_FREE(order);
  V_FREE(vti);

  data->cat = v_get((*N)+1);

  Real cek;
  
  if (QUARTER){
    cek = k/4;
    data->eff = v_get(4*(*N)+1);
  }
  else {
    cek = k/2;
    data->eff = v_get(2*(*N)+1);
  }

  for (int i=0;i<Nce;i++) {
    int idx_e = floor((ti[i]+cek/2)/cek);
    int idx_c = floor((ti[i]+k/2)/k);
    data->cat->ve[idx_c] += (ct[i]/k)/1e3;
    data->eff->ve[idx_e] += 1.0/cek;
  }

  free(ct);
  free(ti);

  fclose(fp);
}

void data_read_ce_new(char * data_file_name, NewData * newdata) {
  char buffer[30];
  sprintf(buffer,"%s-ce-new.dat",data_file_name);

  int Nk,B;
  FILE *fp = fopen(buffer,"r");
  fscanf(fp,"%d",&Nk);
  Real *knots = (Real *) calloc(Nk,sizeof(Real));

  for (int i=0;i<Nk;i++) {
  #if REAL == DOUBLE
    // lf
    fscanf(fp,"%lf", &knots[i]);
  #elif REAL == FLOAT
    // f
    fscanf(fp,"%f", &knots[i]);
  #elif REAL == LONGDOUBLE
    // Lf
    fscanf(fp,"%Lf", &knots[i]);
  #endif
  }

  newdata->nK = Nk;
  newdata->knots_c = v_get(Nk);

  for (int i=0;i<Nk;i++)
    newdata->knots_c->ve[i] = knots[i];
  
  free(knots);

  Real *knots_e = (Real *) calloc(Nk,sizeof(Real));

  for (int i=0;i<Nk;i++) {
  #if REAL == DOUBLE
    // lf
    fscanf(fp,"%lf", &knots_e[i]);
  #elif REAL == FLOAT
    // f
    fscanf(fp,"%f", &knots_e[i]);
  #elif REAL == LONGDOUBLE
    // Lf
    fscanf(fp,"%Lf", &knots_e[i]);
  #endif
  }

  newdata->knots_e = v_get(Nk);

  for (int i=0;i<Nk;i++)
    newdata->knots_e->ve[i] = knots_e[i];
  
  free(knots_e);
  
  fscanf(fp,"%d",&B);
  newdata->B = B;
  
  Real *splcoef = (Real *) calloc(B,sizeof(Real));

  for (int i=0;i<B;i++) {
  #if REAL == DOUBLE
    // lf
    fscanf(fp,"%lf", &splcoef[i]);
  #elif REAL == FLOAT
    // f
    fscanf(fp,"%f", &splcoef[i]);
  #elif REAL == LONGDOUBLE
    // Lf
    fscanf(fp,"%Lf", &splcoef[i]);
  #endif
  }
  
  newdata->splcoef_c = v_get(B);
  for (int i=0;i<B;i++)
    newdata->splcoef_c->ve[i] = splcoef[i];

  free(splcoef);


  Real *splcoef_e = (Real *) calloc(B,sizeof(Real));

  for (int i=0;i<B;i++) {
  #if REAL == DOUBLE
    // lf
    fscanf(fp,"%lf", &splcoef_e[i]);
  #elif REAL == FLOAT
    // f
    fscanf(fp,"%f", &splcoef_e[i]);
  #elif REAL == LONGDOUBLE
    // Lf
    fscanf(fp,"%Lf", &splcoef_e[i]);
  #endif
  }
  
  newdata->splcoef_e = v_get(B);
  for (int i=0;i<B;i++)
    newdata->splcoef_e->ve[i] = splcoef_e[i];

  free(splcoef_e);

  fclose(fp);
}


void data_read_lf_new(char * data_file_name, NewData * newdata) {
  char buffer[30];
  sprintf(buffer,"%s-lf-new.dat",data_file_name);

  int Nlf;
  FILE *fp = fopen(buffer,"r");
  fscanf(fp,"%d",&Nlf);
  Real *ln = (Real *) calloc(Nlf,sizeof(Real));

  for (int i=0;i<Nlf;i++) {
  #if REAL == DOUBLE
    // lf
    fscanf(fp,"%lf", &ln[i]);
  #elif REAL == FLOAT
    // f
    fscanf(fp,"%f", &ln[i]);
  #elif REAL == LONGDOUBLE
    // Lf
    fscanf(fp,"%Lf", &ln[i]);
  #endif
  }

  newdata->Nlf = Nlf;
  newdata->ln = v_get(Nlf);

  for (int i=0;i<Nlf;i++)
    newdata->ln->ve[i] = ln[i];
  
  free(ln);

  Real *tl = (Real *) calloc(Nlf,sizeof(Real));

  for (int i=0;i<Nlf;i++) {
  #if REAL == DOUBLE
    // lf
    fscanf(fp,"%lf", &tl[i]);
  #elif REAL == FLOAT
    // f
    fscanf(fp,"%f", &tl[i]);
  #elif REAL == LONGDOUBLE
    // Lf
    fscanf(fp,"%Lf", &tl[i]);
  #endif
  }

  newdata->tl = v_get(Nlf);

  for (int i=0;i<Nlf;i++)
    newdata->tl->ve[i] = tl[i];
  
  free(tl);
  
  fclose(fp);
}


void data_read_lf(char * data_file_name, Data * data, int N, Real k, int minfish) {
  char buffer[30];
  sprintf(buffer,"%s-lf.dat",data_file_name);

  int Nlf;
  FILE * fp = fopen(buffer,"r");
  fscanf(fp,"%d",&Nlf);

  Real *ln = (Real *) calloc(Nlf,sizeof(Real));
  Real *tl = (Real *) calloc(Nlf,sizeof(Real));

  for (int i=0;i<Nlf;i++) {
    #if REAL == DOUBLE
        // lf
        fscanf(fp,"%lf %lf", &ln[i],&tl[i]);
    #elif REAL == FLOAT
        // f
        fscanf(fp,"%f %f", &ln[i],&tl[i]);
    #elif REAL == LONGDOUBLE
        // Lf
        fscanf(fp,"%Lf %Lf", &ln[i],&tl[i]);
    #endif
  }

  fclose(fp);

  VEC *lfv = v_get(Nlf);
  VEC *ilv = v_get(Nlf);
  VEC *cnt = v_get((int)2*N);

  for (int i=0;i<Nlf;i++) {
    lfv->ve[i] = (Real)ln[i];
    ilv->ve[i] = (int)(N + floor((tl[i]+k/2)/k));
    cnt->ve[(int)(N + floor((tl[i]+k/2)/k))] += 1;
  }

  data->n = 0;
  for (int i=0;i<cnt->dim;i++)  {    
    if (cnt->ve[i] > minfish) {
      data->n += 1;
    }
  }

  data->lf = (Real **) calloc(data->n,sizeof(Real *));
  data->t_id = (int *) calloc(data->n,sizeof(int));
  data->t_sz = (int *) calloc(data->n,sizeof(int));

  int kk=0;
  for (int i=0;i<cnt->dim;i++) {
    if (cnt->ve[i] > minfish) {
      data->t_id[kk] = i;
      data->t_sz[kk] = cnt->ve[i];   
      data->lf[kk] = calloc(data->t_sz[kk],sizeof(Real));

      int jj=0;
      for (int j=0;j<Nlf;j++) {
        if (ilv->ve[j]==i) {
          data->lf[kk][jj] = lfv->ve[j];
          jj += 1;
        }
      }

      kk+=1;
    }
  }

  V_FREE(lfv);
  V_FREE(ilv);
  V_FREE(cnt);

  free(ln);
  free(tl);
}

void optim_control_read(char * optim_file_name, OptimControl * optim) {
  #if REAL == DOUBLE
    // lf
    char * format = "%lf\n";
  #elif REAL == FLOAT
    // f
    char * format = "%f\n";
  #elif REAL == LONGDOUBLE
    // Lf
    char * format = "%Lf\n";
  #endif

  FILE * fp;
  fp = fopen("control.optim","r");
  fscanf(fp,format,&optim->stp);
  fscanf(fp,format,&optim->ftol);
  fscanf(fp,format,&optim->gtol);
  fscanf(fp,format,&optim->stpmin);
  fscanf(fp,format,&optim->stpmax);
  fscanf(fp,"%d",&optim->maxfev);
  fclose(fp);
  optim->xtol = DBL_EPSILON;
}
