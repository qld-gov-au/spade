#include <math.h>
#include "../../meschach/matrix.h"
#include "../../common.h"
#include "../../socbio/fixed/s.h"
#include "../../socbio/fixed/w.h"
#include "../../socbio/variable/effort.h"
#include "../../socbio/variable/catch.h"
#include "../solvers/Q.h"
#include "../../util/util.h"

double H(

		 MAT *x,
		 MAT *u,
		 struct DATA *data,
		 double iota
	
		 )
{

  iota *= 1e-3;
  int S = data->S;
  double k = data->k;

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
	v->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < data->n && data->t_id[lfi]==i) 
	{
      
	  VEC *dt = v_get(data->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = data->lf[lfi][j];

	  double bw = get_bw(dt);

	  VEC *l = v_get(xt->dim);

	  for (int j=0;j<xt->dim;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  double al = c(data->cat,k,k*(i - S)) / Q(xt,l); //1e3*

	  if (PLOT) {

	    FILE *p1 = fopen("plot1.txt","w");

	    for (int j=0;j<x->n;j++)
	      fprintf(p1,"%f %f\n",xt->ve[j],al*l->ve[j]);

	    fclose(p1);

	    FILE *p2 = fopen("plot2.txt","w");

	    for (int j=0;j<x->n;j++)
	      fprintf(p2,"%f %f\n",xt->ve[j],v->ve[j]);

	    fclose(p2);

	    char buffer[100];

	    sprintf(buffer,"./plo > plotl%.3f.pdf",k*(data->t_id[lfi] - S));

	    system(buffer); 

	  }

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
	  ht->ve[i] = pow(Q(xt,v)-c(data->cat,k,k*(i - S)),2.);
	}

      tt->ve[i] = k*(i-S);

      V_FREE(xt);
      V_FREE(ut);
      V_FREE(v);

    }

  if (PLOT) {

    FILE *p1 = fopen("plot.txt","w");

    for (int i=S;i<x->m;i++)
      fprintf(p1,"%f %f\n",tt->ve[i],ht->ve[i]);

    fclose(p1);

    char buffer[100];

    sprintf(buffer,"./plo1 > plotht.pdf");

    system(buffer);

  }

  double blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah; 

}

double G(

	 MAT *p,
	 MAT *x,
	 MAT *u,
	 struct DATA *data,
	 double iota
	
	 )
{

  iota *=1e-3;
  int lfi=0;

  int S = data->S;
  double k = data->k;

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

	  double bw = get_bw(dt);

	  VEC *l = v_get(xt->dim);

	  for (int j=0;j<xt->dim;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  double al = c(data->cat,k,k*(i - S)) / Q(xt,l);

	  if (PLOTDERIV) {

	    FILE *p1 = fopen("plot1.txt","w");

	    for (int j=0;j<x->n;j++)
	      fprintf(p1,"%f %f\n",xt->ve[j],al*l->ve[j]);

	    fclose(p1);

	    FILE *p2 = fopen("plot2.txt","w");

	    for (int j=0;j<x->n;j++)
	      fprintf(p2,"%f %f\n",xt->ve[j],pv->ve[j]);

	    fclose(p2);

	    char buffer[100];

	    sprintf(buffer,"./plo > plotl%.3f.pdf",k*(data->t_id[lfi] - S));

	    system(buffer); 

	  }

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

	  ht->ve[i] = 2*(Q(xt,v)-c(data->cat,k,k*(i-S)))*Q(xt,pv);
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


  if (PLOTDERIV) {

    FILE *p1 = fopen("plot.txt","w");

    for (int i=S;i<x->m;i++)
      fprintf(p1,"%f %f\n",tt->ve[i],ht->ve[i]);

    fclose(p1);

    char buffer[100];

    sprintf(buffer,"./plo1 > plotht.pdf");

    system(buffer);

  }

  double blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah;

}


/*double G_ni_for_condition_number(

	    MAT *p,
	    MAT *x,
	    MAT *u,
	    struct DATA *data,
	    double iota

	    )
{

  int S = data->S;
  double k = data->k;

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

double G_ni(

	    MAT *p,
	    MAT *x,
	    MAT *u,
	    struct DATA *data,
	    double iota

	    )
{

  iota *=1e-3;
  int S = data->S;
  double k = data->k;

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

	  double bw = get_bw(dt);

	  VEC *l = v_get(xt->dim);

	  for (int j=0;j<xt->dim;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  double al = c(data->cat,k,k*(i - S)) / Q(xt,l);

	  if (PLOTDERIV) 
	    {

	      FILE *p1 = fopen("plot1.txt","w");

	      for (int j=0;j<x->n;j++)
		fprintf(p1,"%f %f\n",xt->ve[j],al*l->ve[j]);

	      fclose(p1);

	      FILE *p2 = fopen("plot2.txt","w");

	      for (int j=0;j<x->n;j++)
		fprintf(p2,"%f %f\n",xt->ve[j],pv->ve[j]);

	      fclose(p2);

	      char buffer[100];

	      sprintf(buffer,"./plo > plotl%.3f.pdf",k*(data->t_id[lfi] - S));

	      system(buffer);

	    }

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

	  ht->ve[i] = 2*(Q(xt,v)-c(data->cat,k,k*(i-S)))*Q(xt,pv);
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

  if (PLOTDERIV) {

    FILE *p1 = fopen("plot.txt","w");

    for (int i=S;i<x->m;i++)
      fprintf(p1,"%f %f\n",tt->ve[i],ht->ve[i]);

    fclose(p1);

    char buffer[100];

    sprintf(buffer,"./plo1 > plotht.pdf");

    system(buffer);

  }


  double blah = Q(tt,ht);

  V_FREE(ht);
  V_FREE(tt);

  return blah;

}