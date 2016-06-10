#include <math.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"
#include "spade_solve.h"
#include "../model/fishing/selectivity.h"
#include "../model/biology/weight.h"
#include "../model/fishing/effort.h"
#include "../model/fishing/catch.h"
#include "Q.h"
#include "../util/util.h"

Real K_dr(

  Parameters * parameters,
  Data *d
  
  )
{

  Solve_Core_Args core_args;
  int I = d->I;
  int J = d->J;
    
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

  Real iota = parameters->iota.value;

  iota *= 1e-3;
  int S = d->S;
  Real k = d->k;

  int lfi=0;

  VEC *ht = v_get(x->m);
  VEC *tt = v_get(x->m);

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

    Real bw = get_bw(dt);

    VEC *l = v_get(xt->dim);

    for (int j=0;j<xt->dim;j++)
      for (int jj=0;jj<dt->dim;jj++)
        l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

    Real al = 1e3*c(d->cat,k,k*(i - S)) / Q(xt,l); 

    VEC *ld = v_get(x->n);

    for (int j=0;j<xt->dim;j++)
      ld->ve[j] = pow(v->ve[j] - al*l->ve[j],2.);

    ht->ve[i] = Q(xt,ld);

    lfi += 1;

    V_FREE(dt);
    V_FREE(l);
    V_FREE(ld);

  } 
      else 
  {
    ht->ve[i] = pow(Q(xt,v)-1e3*c(d->cat,k,k*(i - S)),2.);
  }

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(v);

    }

  Real blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  M_FREE(core_args.x);
  M_FREE(core_args.u);
  M_FREE(core_args.xh);
  M_FREE(core_args.uh);
  M_FREE(core_args.xn);
  M_FREE(core_args.xhh);
  M_FREE(core_args.un);
  V_FREE(core_args.Ui);
  V_FREE(core_args.Uh);
  V_FREE(core_args.Uhh);
  IV_FREE(core_args.idxi);

  return blah; 

}


Real K(

     Parameters *parameters,
		 Data *d,
		 Solve_Core_Args *core_args
		 
	
		 )
{


  solve(parameters,d->eff,d->k,d->S,core_args);

  MAT *x = core_args->x;
  MAT *u = core_args->u;

  Real iota = parameters->iota.value;

  iota *= 1e-3;
  int S = d->S;
  Real k = d->k;

  int lfi=0;

  VEC *ht = v_get(x->m);
  VEC *tt = v_get(x->m);

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

	  Real bw = get_bw(dt);

	  VEC *l = v_get(xt->dim);

	  for (int j=0;j<xt->dim;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(d->cat,k,k*(i - S)) / Q(xt,l); //1e3*

	  VEC *ld = v_get(x->n);

	  for (int j=0;j<xt->dim;j++)
	    ld->ve[j] = pow(v->ve[j] - al*l->ve[j],2.);

	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	} 
      else 
	{
	  ht->ve[i] = pow(Q(xt,v)-1e3*c(d->cat,k,k*(i - S)),2.);
	}

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(v);

    }

  Real blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah; 

}

Real G(

	 MAT *p,
	 MAT *x,
	 MAT *u,
	 Data *data,
	 Real iota
	
	 )
{

  iota *=1e-3;
  int lfi=0;

  int S = data->S;
  Real k = data->k;

  VEC *ht = v_get(x->m);
  VEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      VEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      VEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      VEC *pt = v_get(x->n);
      get_row(p,i,pt); 

      VEC *v = v_get(x->n);
      VEC *pv = v_get(x->n);

      for (int j=0;j<x->n;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j] + 1e-3*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      //      if(data->t_id[lfi]==i) 
      if(lfi < data->n && data->t_id[lfi]==i) 
	{
      
	  VEC *dt = v_get(data->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = data->lf[lfi][j];

	  Real bw = get_bw(dt);

	  VEC *l = v_get(xt->dim);

	  for (int j=0;j<xt->dim;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(data->cat,k,k*(i - S)) / Q(xt,l);

	  VEC *ld = v_get(x->n);

	  for (int j=0;j<xt->dim;j++)
	    ld->ve[j] = 2*(v->ve[j] - al*l->ve[j])*pv->ve[j];

	  /*
	    if (v->ve[j] < al*l->ve[j])
	      ld->ve[j] = -pv->ve[j];
	    else
	    ld->ve[j] = pv->ve[j];*/

	  //	    ld->ve[j] = pv->ve[j]*(v->ve[j] - al*l->ve[j]) / fabs(v->ve[j] - al*l->ve[j]);

	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	} 
      else 
	{ 
          //printf("%d\n",i);
          //printf("%f\n",Q(xt,pt));
          //printf("%f\n",Q(xt,v));
          //printf("%f\n",c(k*(i-tmi.I)));

	  ht->ve[i] = 2*(Q(xt,v)-1e3*c(data->cat,k,k*(i-S)))*Q(xt,pv);
	  /*    
	  if (Q(xt,v)<1e3*c(data->cat,k,k*(i-S)))
	    ht->ve[i] = -Q(xt,pv);//*(Q(xt,v)-1e3*c(k*(i - tmi.I)))/fabs(Q(xt,v)-1e3*c(k*(i - tmi.I)));
	  else
	  ht->ve[i]=Q(xt,pv);*/
	}

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(pt);
      V_FREE(v);
      V_FREE(pv);

    }

  Real blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah;

}


/*Real G_ni_for_condition_number(

	    MAT *p,
	    MAT *x,
	    MAT *u,
	    Data *data,
	    Real iota

	    )
{

  int S = data->S;
  Real k = data->k;

  int lfi=0;

  VEC *ht = v_get(x->m);
  VEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      VEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      VEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      VEC *pt = v_get(x->n);
      get_row(p,i,pt); 

      VEC *v = v_get(x->n);
      VEC *pv = v_get(x->n);

      for (int j=0;j<x->n;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      if(data->t_id[lfi]==i) 
	{
      
	  VEC *l = v_get(xt->dim);
	  for (int j=0;j<xt->dim;j++)
	    l->ve[j] = ut->ve[j]; 

	  for (int j=0;j<xt->dim;j++)
	    if (v->ve[j] < al*l->ve[j])
	      ld->ve[j] = -pv->ve[j];
	    else
	      ld->ve[j] = pv->ve[j];

	  ht->ve[i] = Q(xt,ld);

          if (lfi<data->n)
	    lfi += 1;

	} 
      else 
	{ 
          //printf("%d\n",i);
          //printf("%f\n",Q(xt,pt));
          //printf("%f\n",Q(xt,v));
          //printf("%f\n",c(k*(i-tmi.I)));
        
	  if (Q(xt,v)<1e3*c(data->cat,k,k*(i-S)))
	    ht->ve[i] = -Q(xt,pv);//*(Q(xt,v)-1e3*c(k*(i - tmi.I)))/fabs(Q(xt,v)-1e3*c(k*(i - tmi.I)));
	  else
	    ht->ve[i]=Q(xt,pv);
	}

      tt->ve[i] = k*(i-S);

    }


  if (BLAH) {

    FILE *p1 = fopen("plot.txt","w");

    for (int i=S;i<x->m;i++)
      fprintf(p1,"%f %f\n",tt->ve[i],ht->ve[i]);

    fclose(p1);

    char buffer[100];

    sprintf(buffer,"./plo1 > plotht.pdf");

    system(buffer);

  }

  return Q(tt,ht);

}
*/

Real G_ni(

	    MAT *p,
	    MAT *x,
	    MAT *u,
	    Data *data,
	    Real iota

	    )
{

  iota *=1e-3;
  int S = data->S;
  Real k = data->k;

  int lfi=0;

  VEC *ht = v_get(x->m);
  VEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      VEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      VEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      VEC *pt = v_get(x->n);
      get_row(p,i,pt); 

      VEC *v = v_get(x->n);
      VEC *pv = v_get(x->n);

      for (int j=0;j<x->n;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      if(lfi < data->n && data->t_id[lfi]==i) 
	{
      
	  VEC *dt = v_get(data->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = data->lf[lfi][j];

	  Real bw = get_bw(dt);

	  VEC *l = v_get(xt->dim);

	  for (int j=0;j<xt->dim;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(data->cat,k,k*(i - S)) / Q(xt,l);

	  VEC *ld = v_get(x->n);

	  for (int j=0;j<xt->dim;j++)
	    ld->ve[j] = 2*(v->ve[j] - al*l->ve[j])*pv->ve[j];

	  /*
	  for (int j=0;j<xt->dim;j++)
	    if (v->ve[j] < al*l->ve[j])
	      ld->ve[j] = -pv->ve[j];
	    else
	    ld->ve[j] = pv->ve[j];*/

	  // Ld->ve[j] = pv->ve[j]*(v->ve[j] - al*l->ve[j]) / fabs(v->ve[j] - al*l->ve[j]);

	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	} 
      else 
	{ 
          //printf("%d\n",i);
          //printf("%f\n",Q(xt,pt));
          //printf("%f\n",Q(xt,v));
          //printf("%f\n",c(k*(i-tmi.I)));

	  ht->ve[i] = 2*(Q(xt,v)-1e3*c(data->cat,k,k*(i-S)))*Q(xt,pv);
	  /*        
	  if (Q(xt,v)<1e3*c(data->cat,k,k*(i-S)))
	    ht->ve[i] = -Q(xt,pv);//*(Q(xt,v)-1e3*c(k*(i - tmi.I)))/fabs(Q(xt,v)-1e3*c(k*(i - tmi.I)));
	  else
	  ht->ve[i]=Q(xt,pv);*/
	}

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(pt);
      V_FREE(v);
      V_FREE(pv);

    }

  Real blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah;

}
