#include <math.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "../VMGMM/solvers/spade_solve.h"
#include "../socbio/fixed/selectivity.h"
#include "../socbio/fixed/weight.h"
#include "../socbio/variable/effort.h"
#include "../socbio/variable/catch.h"
#include "../VMGMM/solvers/Q.h"
#include "../util/util.h"

void plot( 

	  VEC *theta,
	  struct DATA *d,
	  char * label

	   )
{

  int I = d->I+1;
  int J = d->J+1;
  MAT *x = m_get(I,J);
  MAT *u = m_get(I,J);
  MAT *xh = m_get(I,J+1);
  MAT *uh = m_get(I,J+1);
  MAT *xn = m_get(I,J+1);
  MAT *xhh = m_get(I,J+1);
  MAT *un = m_get(I,J+1);
  VEC *Ui = v_get(I);
  VEC *Uh = v_get(I);
  VEC *Uhh = v_get(I);
  IVEC *idxi = iv_get(I-1);

  solve(theta,x,u,xhh,xh,xn,uh,un,Ui,Uh,Uhh,idxi,d->eff,d->k,d->S);

  M_FREE(xh);
  M_FREE(uh);
  M_FREE(xn);
  M_FREE(xhh);
  M_FREE(un);
  V_FREE(Ui);
  V_FREE(Uh);
  V_FREE(Uhh);
  IV_FREE(idxi);

  double iota = theta->ve[3];
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

  M_FREE(x);
  M_FREE(u);

}
