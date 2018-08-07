// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"
#include "../machinery/spade_solve.h"
#include "../model/fishing/selectivity.h"
#include "../model/fishing/effort.h"
#include "../model/fishing/catch.h"
#include "../model/biology/weight.h"
#include "../machinery/Q.h"
#include "../util/util.h"

void plot( 

	  Parameters *parameters,
	  Data *d,
	  char * label

	   )
{

  Solve_Core_Args core_args;
  
  int I,J;
  I = d->I;
  J = d->J;

  core_args.x = m_get(I+1,J+1);
  core_args.u = m_get(I+1,J+1);

  // todo: Review this. Condition has been flipped to maintain same functionality as last commit
  // in non-bigmatrices mode but may not be correct for bigmatrices mode. Affects xh in particular
  // as it was of dimension 401 here but was of dimension 402 in working copy.
  if (SGNM)
    {
       core_args.xh = m_get(I+1,J+2);
       core_args.uh = m_get(I+1,J+2);
       core_args.xn = m_get(I+1,J+2);
       core_args.xhh = m_get(I+1,J+2);
       core_args.un = m_get(I+1,J+2);
     }
  else
    {
       core_args.xh = m_get(I+1,J+1);
       core_args.uh = m_get(I+1,J+1);
       core_args.xn = m_get(I+1,J+1);
       core_args.xhh = m_get(I+1,J+1);
       core_args.un = m_get(I+1,J+1);
    }

  core_args.Ui = v_get(I+1);
  core_args.Uh = v_get(I+1);
  core_args.Uhh = v_get(I+1);
  core_args.idxi = iv_get(I);

  solve(parameters,d->eff,d->k,d->S,d->Y,&core_args);

  MAT *x = core_args.x;
  MAT *u = core_args.u;

  M_FREE(core_args.xh);
  M_FREE(core_args.uh);
  M_FREE(core_args.xn);
  M_FREE(core_args.xhh);
  M_FREE(core_args.un);
  V_FREE(core_args.Ui);
  V_FREE(core_args.Uh);
  V_FREE(core_args.Uhh);
  IV_FREE(core_args.idxi);

  Real iota = parameters->iota.value;
  iota *= 1e-3;

  int bigJ;
  if (!SGNM){
    bigJ = x->n - 1;
    J = x->n - x->m;

  } else {
    J = x->n - 1;   
  }
  
  FILE *sp1 = fopen("plot1.txt","w");

  for (int i=0;i<x->m;i++)
    {

      int terminator;
      if(!SGNM)
	terminator = J+i;
      else
	terminator = J;

      VEC *ctt = v_get(terminator+1);
      VEC *xt = v_get(terminator+1);      
      
      xt = get_row(x,i,xt);

      if (!SGNM)
	xt = v_resize(xt,terminator+1);
      
      for (int j=0;j<=terminator;j++)
	ctt->ve[j] = s(x->me[i][j])*iota*e(d->eff,d->k,d->k*(i-d->S),d->Y)*w(x->me[i][j])*u->me[i][j];
      fprintf(sp1,"%lf %lf\n",d->k*(i-d->S),Q(xt,ctt));

      V_FREE(ctt);
      V_FREE(xt);

    }

  fclose(sp1);

  FILE *sp2 = fopen("plot2.txt","w");

  for (int i=0;i<x->m;i++) 
    fprintf(sp2,"%lf %lf\n",d->k*(i-d->S),1e3*c(d->cat,d->k,d->k*(i-d->S)));

  fclose(sp2);

  char sbuffer[100];
  sprintf(sbuffer,"./plo > plotc_%s.pdf",label);
  system(sbuffer); 
  //free(sbuffer);

  int S = d->S;
  Real k = d->k;

  int lfi=0;

  for (int i=S;i<x->m;i++)
    {

      int terminator;
      if(!SGNM)
	terminator = J+i;
      else
	terminator = J;


      VEC *xt = v_get(terminator+1);
      get_row(x,i,xt);      
      VEC *ut = v_get(terminator+1);
      get_row(u,i,ut);
      
      if (!SGNM)
	{
	  xt = v_resize(xt,terminator+1);
	  ut = v_resize(ut,terminator+1);
	}
      
      VEC *v = v_get(terminator+1);

      for (int j=0;j<=terminator;j++)
	v->ve[j] = iota*e(d->eff,k,k*(i-S),d->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < d->n && d->t_id[lfi]==i) 
	{
      
	  VEC *dt = v_get(d->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = d->lf[lfi][j];

	  Real bw = get_bw(dt);

	  VEC *l = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(d->cat,k,k*(i - S)) / Q(xt,l); 

	  FILE *p1 = fopen("plot1.txt","w");

	  for (int j=0;j<=terminator;j++)
	    fprintf(p1,"%lf %lf\n",xt->ve[j],al*l->ve[j]);

	  fclose(p1);

	  FILE *p2 = fopen("plot2.txt","w");

	  for (int j=0;j<=terminator;j++)
	    fprintf(p2,"%lf %lf\n",xt->ve[j],v->ve[j]);

	  fclose(p2);

	  char buffer[100];

	  sprintf(buffer,"./plo > plotl%.3lf_%s.pdf",k*(d->t_id[lfi] - S),label);	

	  system(buffer);

	  //free(buffer);

	  V_FREE(dt);

	  lfi += 1;

	}

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(v);

    }

  M_FREE(core_args.x);
  M_FREE(core_args.u);

}
