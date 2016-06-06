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

  int I = d->I+1;
  int J = d->J+1;

  Solve_Core_Args core_args;
  
  core_args.x = m_get(I,J);
  core_args.u = m_get(I,J);
  core_args.xh = m_get(I,J+1);
  core_args.uh = m_get(I,J+1);
  core_args.xn = m_get(I,J+1);
  core_args.xhh = m_get(I,J+1);
  core_args.un = m_get(I,J+1);
  core_args.Ui = v_get(I);
  core_args.Uh = v_get(I);
  core_args.Uhh = v_get(I);
  core_args.idxi = iv_get(I-1);  

  solve(parameters,d->eff,d->k,d->S,&core_args);

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

  double iota = parameters->iota.value;
  iota *= 1e-3;

  VEC *ctt = v_get(x->n);
  VEC *xt = v_get(x->n);

  FILE *sp1 = fopen("plot1.txt","w");

  for (int i=0;i<x->m;i++)
    {

      xt = get_row(x,i,xt);
      for (int j=0;j<x->n;j++)
	ctt->ve[j] = s(x->me[i][j])*iota*e(d->eff,d->k,d->k*(i-d->S))*w(x->me[i][j])*u->me[i][j];
      fprintf(sp1,"%f %f\n",d->k*(i-d->S),Q(xt,ctt));

    }

  fclose(sp1);

  FILE *sp2 = fopen("plot2.txt","w");

  for (int i=0;i<x->m;i++) 
    fprintf(sp2,"%f %f\n",d->k*(i-d->S),c(d->cat,d->k,d->k*(i-d->S)));

  fclose(sp2);

  char sbuffer[100];
  sprintf(sbuffer,"./plo > plotc_%s.pdf",label);
  system(sbuffer); 
  //free(sbuffer);

  V_FREE(ctt);
  V_FREE(xt);

  int S = d->S;
  double k = d->k;

  int lfi=0;

  for (int i=S;i<x->m;i++)
    {

      VEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      VEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      VEC *v = v_get(x->n);

      for (int j=0;j<x->n;j++)
	v->ve[j] = iota*e(d->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < d->n && d->t_id[lfi]==i) 
	{
      
	  VEC *dt = v_get(d->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = d->lf[lfi][j];

	  double bw = get_bw(dt);

	  VEC *l = v_get(xt->dim);

	  for (int j=0;j<xt->dim;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  double al = c(d->cat,k,k*(i - S)) / Q(xt,l); 

	  FILE *p1 = fopen("plot1.txt","w");

	  for (int j=0;j<x->n;j++)
	    fprintf(p1,"%f %f\n",xt->ve[j],al*l->ve[j]);

	  fclose(p1);

	  FILE *p2 = fopen("plot2.txt","w");

	  for (int j=0;j<x->n;j++)
	    fprintf(p2,"%f %f\n",xt->ve[j],v->ve[j]);

	  fclose(p2);

	  char buffer[100];

	  sprintf(buffer,"./plo > plotl%.3f_%s.pdf",k*(d->t_id[lfi] - S),label);	

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
