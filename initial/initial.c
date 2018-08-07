// Copyright 2016 State of Queensland
// This file is part of SPADE
// See spade.c, COPYING, COPYING.LESSER

#include <math.h>
#include "../meschach/matrix.h"
#include "../model/fishing/selectivity.h"
#include "../common.h"
#include "../machinery/Qn.h"
#include "../machinery/Q.h"
#include "../util/util.h"

VEC * calc_alpha(

			Real a,
			Real k,
			Real w,
			Real bt,
			Real f

			)
{

  VEC *retur = v_get(3);

  VEC *x = v_get(500);
  VEC *v = v_get(x->dim);
  VEC *t = v_get(x->dim);

  for (int j=0;j<x->dim;j++)
    x->ve[j] = j*w/500.;

  for (int j=0;j<t->dim;j++)
    t->ve[j] = (bt + f*s(x->ve[j])) / (k*(w-x->ve[j]));

  for (int j=0;j<v->dim;j++)
    v->ve[j] = a*(A1*x->ve[j] + A2*pow(x->ve[j],2.))*exp(-Qn(x,t,j+1))/(k*(w-x->ve[j]));

  Real qv = Q(x,v);

  retur->ve[0] = pow(qv-1,2.);

  VEC *z = v_get(x->dim);
  VEC *q = v_get(x->dim);

  for (int j=0;j<v->dim;j++)
    z->ve[j] = pow(x->ve[j],2.)*exp(-Qn(x,t,j+1))/(k*(w-x->ve[j]));

  retur->ve[1] = 2*Q(x,z)*(qv-1); 

  return retur;

}

VEC * VMGMM_eq(

		  VEC *x,
		  Data *d,
		  VEC *grad,
		  Real *f

		  )
{

  Real kappa = .1; // FIX THIS
  Real omega = 160; // FIX THIS

  VEC *xx = v_get(d->J);
  for (int j=0;j<xx->dim;j++)
    xx->ve[j] = h*j;

  int nd=0;
  for (int i=0;i<d->n;i++)
    nd += d->t_sz[i];

  VEC *dt = v_get(nd);
  int ii=0;
  for (int i=0;i<d->n;i++)
    for (int jj=0;jj<d->t_sz[i];jj++)
      {
	dt->ve[ii] = d->lf[i][jj];
	ii++;
      }

  Real bw = get_bw(dt);

  VEC *l = v_get(d->J);

  for (int j=0;j<l->dim;j++)
    for (int jj=0;jj<dt->dim;jj++)
      l->ve[j] += exp( -pow((xx->ve[j] - dt->ve[jj])/bw,2.) );

  Real Ql = Q(xx,l);

  for (int j=0;j<l->dim;j++)
    l->ve[j] /= Ql;

  VEC *v = v_get(xx->dim);
  VEC *w = v_get(xx->dim);

  for (int j=0;j<w->dim;j++)
    w->ve[j] = (x->ve[0] + x->ve[1]*s(xx->ve[j])) / (kappa*(omega-xx->ve[j]));

  for (int j=1;j<v->dim;j++)
    v->ve[j] = s(xx->ve[j])*exp(-Qn(xx,w,j+1))/(omega-xx->ve[j]);
      
  Real B = Q(xx,v);

  VEC * vn = v_get(v->dim);

  for (int j=0;j<v->dim;j++)
    vn->ve[j] = v->ve[j]/B;
  /*
  printf("\n");
  for (int j=0;j<w->dim;j++)
    printf("%f %f\n",xx->ve[j],vn->ve[j]);
  printf("e\n\n");
  for (int j=0;j<w->dim;j++)
    printf("%f %f\n",xx->ve[j],l->ve[j]);
  exit(1);
  */
  VEC *objf = v_get(v->dim);

  for (int j=0;j<v->dim;j++)
    objf->ve[j]=l->ve[j]*log(l->ve[j]/(vn->ve[j]+1e-12) + 1e-12);
    
  *f = Q(xx,objf);

  VEC *w2 = v_get(w->dim);

  for (int j=0;j<w->dim-1;j++)
    w2->ve[j] = 1/(kappa*(omega-xx->ve[j]));

  VEC *w3 = v_get(w->dim);
  for (int j=1;j<w->dim;j++)
    w3->ve[j] = s(xx->ve[j])*exp(-Qn(xx,w,j+1))*Qn(xx,w2,j+1)/(omega-xx->ve[j]);

  Real C = Q(xx,w3);

  VEC * fc = v_get(w->dim);

  for (int j=1;j<w->dim;j++)
    fc->ve[j] = (B*Qn(xx,w2,j+1)-C ) / B;
    
  VEC *w4 = v_get(w->dim);
  
  for (int j=0;j<w->dim;j++)
    w4->ve[j] =  s(xx->ve[j])/(kappa*(omega-xx->ve[j]));

  VEC *w5 = v_get(w->dim);

  for (int j=0;j<w->dim;j++)
    w5->ve[j] = s(xx->ve[j])*exp(-Qn(xx,w,j+1))*Qn(xx,w4,j+1)/(omega-xx->ve[j]);

  Real D = Q(xx,w5);

  VEC * fc2 = v_get(w->dim);

  for (int j=1;j<w->dim;j++)
    fc2->ve[j] = (B*Qn(xx,w4,j+1)-D) / B;

  VEC *dGdb = v_get(w->dim);
  VEC *dGdf = v_get(w->dim);
    
  for (int j=0;j<w->dim;j++)
    {
      dGdb->ve[j]=l->ve[j]*fc->ve[j];
      dGdf->ve[j]=l->ve[j]*fc2->ve[j];
    }

  grad->ve[0] = Q(xx,dGdb);
  grad->ve[1] = Q(xx,dGdf);

  V_FREE(v);
  V_FREE(w);
  V_FREE(w2);
  V_FREE(w3);
  V_FREE(w4);
  V_FREE(w5);
  V_FREE(dGdb);
  V_FREE(dGdf);

  return grad;

}
