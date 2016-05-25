﻿#include <math.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "optim.h"

VEC * bfgs(

	   VEC * (*model)(VEC *,struct DATA *,VEC *,double *),
	   VEC *x,
	   struct DATA *data

	   )   
{

  int n = x->dim;
  int nfev = 0;

  double f;

  VEC *s = v_get(n);
  VEC *y = v_get(n);
  VEC *g = v_get(n);
  VEC *p = v_get(n);
  VEC *u = v_get(n);

  VEC *oldx = v_get(n);
  VEC *oldg = v_get(n);

  MAT *B = m_get(n,n);

  m_ident(B);

  g = (*model)(x,data,g,&f); 
  nfev += 1;

  while (1)
    {

      printf("v_norm_inf(g): %f\n",v_norm_inf(g));

      if (v_norm_inf(g) < 1e-4)
	break;

      mv_mlt(B,g,p);      
      sv_mlt(-1./v_norm2(p),p,p);

      v_copy(x,oldx);
      v_copy(g,oldg);

      // More-Thuente line search
      nfev += cvsrch(model,x,f,g,p,1e-2,1e-1,0.4,DBL_EPSILON,1e-20,1e20,60,data);
    
      v_sub(x,oldx,s);
      v_sub(g,oldg,y);

      // update inverse hessian based on julia code
      double sy = in_prod(s,y);

      if (sy==0)
	{
	  printf("cannot find a stepsize that reduces the function along the descent direction\n");
	  exit(1);
	}

      mv_mlt(B,y,u);

      double yBy = in_prod(u,y);
      double c1 = (sy + yBy) / (sy*sy);
      double c2 = 1/sy;

      // not using meschach for this - overcomplicates it.
      for (int i=0;i<n;i++)
	for (int j=0;j<n;j++)
	  B->me[i][j] += c1 * s->ve[i] * s->ve[j] - c2 * ( u->ve[i] * s->ve[j] + u->ve[j] * s->ve[i] );

    }

  (*model)(x,data,g,&f); 
  printf("number function evals: %d, function value: %f\n",nfev,f);
  v_output(x);

  V_FREE(oldx);
  V_FREE(oldg);
  V_FREE(s);
  V_FREE(y);
  V_FREE(g);
  M_FREE(B);
  V_FREE(u);

  return(x);

}

int cstep(

	  double *stx,
	  double *fx,
	  double *dx,
	  double *sty,
	  double *fy,
	  double *dy,
	  double *stp,
	  double fp,
	  double dp,
	  int *brackt,
	  double stpmin,
	  double stpmax

	  )
{

  /*

    Translation of Dianne O'Leary's translation (to matlab) of minpack subroutine cstep.

    Dianne's comments will be prefaced by the matlab comment specifier '%'

    %   Translation of minpack subroutine cstep 
    %   Dianne O'Leary   July 1991
    %     **********
    %
    %     Subroutine cstep
    %
    %     The purpose of cstep is to compute a safeguarded step for
    %     a linesearch and to update an interval of uncertainty for
    %     a minimizer of the function.
    %
    %     The parameter stx contains the step with the least function
    %     value. The parameter stp contains the current step. It is
    %     assumed that the derivative at stx is negative in the
    %     direction of the step. If brackt is set true then a
    %     minimizer has been bracketed in an interval of uncertainty
    %     with endpoints stx and sty.
    %
    %     The subroutine statement is
    %
    %       subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
    %                        stpmin,stpmax,info)
    % 
    %     where
    %
    %       stx, fx, and dx are variables which specify the step,
    %         the function, and the derivative at the best step obtained
    %         so far. The derivative must be negative in the direction
    %         of the step, that is, dx and stp-stx must have opposite 
    %         signs. On output these parameters are updated appropriately.
    %
    %       sty, fy, and dy are variables which specify the step,
    %         the function, and the derivative at the other endpoint of
    %         the interval of uncertainty. On output these parameters are 
    %         updated appropriately.
    %
    %       stp, fp, and dp are variables which specify the step,
    %         the function, and the derivative at the current step.
    %         If brackt is set true then on input stp must be
    %         between stx and sty. On output stp is set to the new step.
    %
    %       brackt is a logical variable which specifies if a minimizer
    %         has been bracketed. If the minimizer has not been bracketed
    %         then on input brackt must be set false. If the minimizer
    %         is bracketed then on output brackt is set true.
    %
    %       stpmin and stpmax are input variables which specify lower 
    %         and upper bounds for the step.
    %
    %       info is an integer output variable set as follows:
    %         If info = 1,2,3,4,5, then the step has been computed
    %         according to one of the five cases below. Otherwise
    %         info = 0, and this indicates improper input parameters.
    %
    %     Subprograms called
    %
    %       FORTRAN-supplied ... abs,max,min,sqrt
    %                        ... dble
    %
    %     Argonne National Laboratory. MINPACK Project. June 1983
    %     Jorge J. More', David J. Thuente
    %
    %     **********
  */

  int info = 0;

  //% Determine if the derivatives have opposite sign.

  double sgnd = dp*((*dx)/fabs((*dx)));

  //% First case. A higher function value.
  //% The minimum is bracketed. If the cubic step is closer to stx than the quadratic step, the cubic step is taken, else the average of the cubic and quadratic steps is taken.

  int bound;
  double theta;
  double ss;
  double gamma;
  double stpf;

  if (fp > (*fx)) 
    {

      info = 1;
      bound = 1;
      theta = 3*((*fx) - fp)/((*stp) - (*stx)) + (*dx) + dp;
      VEC *tmp = v_get(3);
      tmp->ve[0] = theta;
      tmp->ve[1] = (*dx);
      tmp->ve[2] = dp;
      ss = v_norm_inf(tmp);
      V_FREE(tmp);
      gamma = ss*sqrt(pow(theta/ss,2.) - ((*dx)/ss)*(dp/ss));
      if ((*stp) < (*stx)) 
	gamma = -gamma;

      double p = (gamma - (*dx)) + theta;
      double q = ((gamma - (*dx)) + gamma) + dp;
      double r = p/q;
      double stpc = (*stx) + r*((*stp) - (*stx));
      double stpq = (*stx) + (((*dx)/(((*fx)-fp)/((*stp)-(*stx))+(*dx)))/2)*((*stp) - (*stx));
      if (fabs(stpc-(*stx)) < fabs(stpq-(*stx)))
	stpf = stpc;
      else
	stpf = stpc + (stpq - stpc)/2;
         
      *brackt = 1;

    } 

  //% Second case. A lower function value and derivatives of opposite sign. The minimum is bracketed. If the cubic step is closer to stx than the quadratic (secant) step, the cubic step is taken, else the quadratic step is taken.

  else if (sgnd < 0.0) 
    {
      info = 2;
      bound = 0;

      theta = 3*((*fx) - fp)/((*stp) - (*stx)) + (*dx) + dp;
      VEC *tmp = v_get(3);
      tmp->ve[0] = theta;
      tmp->ve[1] = (*dx);
      tmp->ve[2] = dp;
      ss = v_norm_inf(tmp);
      V_FREE(tmp);
      gamma = ss*sqrt(pow(theta/ss,2.) - ((*dx)/ss)*(dp/ss));

      if ((*stp) > (*stx)) 
	gamma = -gamma;
         
      double p = (gamma - dp) + theta;
      double q = ((gamma - dp) + gamma) + (*dx);
      double r = p/q;

      double stpc = (*stp) + r*((*stx) - (*stp));
      double stpq = (*stp) + (dp/(dp-(*dx)))*((*stx) - (*stp));
      if (fabs(stpc-(*stp)) > fabs(stpq-(*stp)))
	stpf = stpc;
      else
	stpf = stpq;
         
      *brackt = 1;
    }

  //% Third case. A lower function value, derivatives of the same sign, and the magnitude of the derivative decreases. The cubic step is only used if the cubic tends to infinity in the direction of the step or if the minimum of the cubic is beyond stp. Otherwise the cubic step is defined to be either stpmin or stpmax. The quadratic (secant) step is also computed and if the minimum is bracketed then the the step closest to stx is taken, else the step farthest away is taken.

  else if (fabs(dp) < fabs(*dx)) 
    {

      info = 3;
      bound = 1;

      theta = 3*((*fx) - fp)/((*stp) - (*stx)) + (*dx) + dp;
      VEC *tmp = v_get(3);
      tmp->ve[0] = theta;
      tmp->ve[1] = (*dx);
      tmp->ve[2] = dp;
      ss = v_norm_inf(tmp);
      V_FREE(tmp);
      // % The case gamma = 0 only arises if the cubic does not tend to infinity in the direction of the step.

      gamma = ss*sqrt(max(0.,pow(theta/ss,2.) - (*dx/ss)*(dp/ss)));

      if ((*stp) > (*stx)) 
	gamma = -gamma;
         
      double p = (gamma - dp) + theta;
      double q = (gamma + ((*dx) - dp)) + gamma;
      double r = p/q;

      double stpc,stpq;

      if (r < 0.0 && gamma != 0.0)
	stpc = (*stp) + r*((*stx) - (*stp));
      else if ((*stp) > (*stx))
	stpc = stpmax;
      else
	stpc = stpmin;

      stpq = (*stp) + (dp/(dp-(*dx)))*((*stx) - (*stp));
      if (*brackt)
	if (fabs((*stp)-stpc) < fabs((*stp)-stpq))
	  stpf = stpc;
	else
	  stpf = stpq;           
      else
	if (fabs((*stp)-stpc) > fabs((*stp)-stpq))
	  stpf = stpc;
	else
	  stpf = stpq;
    }

  //% Fourth case. A lower function value, derivatives of the same sign, and the magnitude of the derivative does not decrease. If the minimum is not bracketed, the step is either stpmin or stpmax, else the cubic step is taken.

  else
    {
      info = 4;
      bound = 0;

      if (*brackt)
	{

	  theta = 3*(fp - (*fy))/((*sty) - *stp) + (*dy) + dp;
	  VEC *tmp = v_get(3);
	  tmp->ve[0] = theta;
	  tmp->ve[1] = (*dy);
	  tmp->ve[2] = dp;
	  ss = v_norm_inf(tmp);
	  V_FREE(tmp);

	  gamma = ss*sqrt(pow(theta/ss,2.) - ((*dy)/ss)*(dp/ss));

	  if (*stp > (*sty))
	    gamma = -gamma;
            
          double p = (gamma - dp) + theta;
          double q = ((gamma - dp) + gamma) + (*dy);
          double r = p/q;
          double stpc = *stp + r*((*sty) - *stp);
          stpf = stpc;
	}
      else if (*stp > (*stx))
	stpf = stpmax;
      else
	stpf = stpmin;
    
    }

  //% Update the interval of uncertainty. This update does not depend on the new step or the case analysis above.

  if (fp > *fx)
    {
      *sty = *stp;
      *fy = fp;
      *dy = dp;
    }
  else
    {
      if (sgnd < 0.0)
	{
	  *sty = *stx;
	  *fy = *fx;
	  *dy = *dx;
	}
          
      *stx = *stp;
      *fx = fp;
      *dx = dp;
    }

  //% Compute the new step and safeguard it.

  stpf = min(stpmax,stpf);
  stpf = max(stpmin,stpf);
  *stp = stpf;
  if (*brackt && bound)
    if (*sty > *stx) 
      *stp = min(*stx+.66*(*sty-*stx),*stp);
    else
      *stp = max(*stx+.66*(*sty-*stx),*stp);
         
  return info;

  //%  last card of subroutine cstep

}

int cvsrch(

	  VEC *(*fcn)(VEC *,struct DATA *,VEC *,double *),
	  VEC *x,
	  double f,
	  VEC *gr,
	  VEC *sd,
	  double stp,
	  double ftol,
	  double gtol,
	  double xtol,
	  double stpmin,
	  double stpmax,
	  int maxfev,
	  struct DATA *d

	  )
{
  /* 

    Translation of Dianne O'Leary's translation (into matlab) of minpack subroutine cvsrch.
    
    Dianne's comments will be prefaced with the matlab comment specifier '%'    
 
    %   Translation of minpack subroutine cvsrch
    %   Dianne O'Leary   July 1991
    %     **********
    %
    %     Subroutine cvsrch
    %
    %     The purpose of cvsrch is to find a step which satisfies 
    %     a sufficient decrease condition and a curvature condition.
    %     The user must provide a subroutine which calculates the
    %     function and the gradient.
    %
    %     At each stage the subroutine updates an interval of
    %     uncertainty with endpoints stx and sty. The interval of
    %     uncertainty is initially chosen so that it contains a 
    %     minimizer of the modified function
    %
    %          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
    %
    %     If a step is obtained for which the modified function 
    %     has a nonpositive function value and nonnegative derivative, 
    %     then the interval of uncertainty is chosen so that it 
    %     contains a minimizer of f(x+stp*s).
    %
    %     The algorithm is designed to find a step which satisfies 
    %     the sufficient decrease condition 
    %
    %           f(x+stp*s) <= f(x) + ftol*stp*(gradf(x)'s),
    %
    %     and the curvature condition
    %
    %           abs(gradf(x+stp*s)'s)) <= gtol*abs(gradf(x)'s).
    %
    %     If ftol is less than gtol and if, for example, the function
    %     is bounded below, then there is always a step which satisfies
    %     both conditions. If no step can be found which satisfies both
    %     conditions, then the algorithm usually stops when rounding
    %     errors prevent further progress. In this case stp only 
    %     satisfies the sufficient decrease condition.
    %
    %     The subroutine statement is
    %
    %        subroutine cvsrch(fcn,n,x,f,g,s,stp,ftol,gtol,xtol,
    %                          stpmin,stpmax,maxfev,info,nfev,wa)
    %     where
    %
    %	fcn is the name of the user-supplied subroutine which
    %         calculates the function and the gradient.  fcn must 
    %      	  be declared in an external statement in the user 
    %         calling program, and should be written as follows.
    %
    %         function [f,g] = fcn(n,x) (Matlab)     (10/2010 change in documentation)
    %	  (derived from Fortran subroutine fcn(n,x,f,g) )
    %         integer n
    %         f
    %         x(n),g(n)
    %	  ----------
    %         Calculate the function at x and
    %         return this value in the variable f.
    %         Calculate the gradient at x and
    %         return this vector in g.
    %	  ----------
    %	  return
    %	  end
    %
    %       n is a positive integer input variable set to the number
    %	  of variables.
    %
    %	x is an array of length n. On input it must contain the
    %	  base point for the line search. On output it contains 
    %         x + stp*s.
    %
    %	f is a variable. On input it must contain the value of f
    %         at x. On output it contains the value of f at x + stp*s.
    %
    %	g is an array of length n. On input it must contain the
    %         gradient of f at x. On output it contains the gradient
    %         of f at x + stp*s.
    %
    %	s is an input array of length n which specifies the
    %         search direction.
    %
    %	stp is a nonnegative variable. On input stp contains an
    %         initial estimate of a satisfactory step. On output
    %         stp contains the final estimate.
    %
    %       ftol and gtol are nonnegative input variables. Termination
    %         occurs when the sufficient decrease condition and the
    %         directional derivative condition are satisfied.
    %
    %	xtol is a nonnegative input variable. Termination occurs
    %         when the relative width of the interval of uncertainty 
    %	  is at most xtol.
    %
    %	stpmin and stpmax are nonnegative input variables which 
    %	  specify lower and upper bounds for the step.
    %
    %	maxfev is a positive integer input variable. Termination
    %         occurs when the number of calls to fcn is at least
    %         maxfev by the end of an iteration.
    %
    %	info is an integer output variable set as follows:
    %	  
    %	  info = 0  Improper input parameters.
    %
    %	  info = 1  The sufficient decrease condition and the
    %                   directional derivative condition hold.
    %
    %	  info = 2  Relative width of the interval of uncertainty
    %		    is at most xtol.
    %
    %	  info = 3  Number of calls to fcn has reached maxfev.
    %
    %	  info = 4  The step is at the lower bound stpmin.
    %
    %	  info = 5  The step is at the upper bound stpmax.
    %
    %	  info = 6  Rounding errors prevent further progress.
    %                   There may not be a step which satisfies the
    %                   sufficient decrease and curvature conditions.
    %                   Tolerances may be too small.
    %
    %       nfev is an integer output variable set to the number of
    %         calls to fcn.
    %
    %	wa is a work array of length n.
    %
    %     Subprograms called
    %
    %	user-supplied......fcn
    %
    %	MINPACK-supplied...cstep
    %
    %	FORTRAN-supplied...abs,max,min
    %	  
    %     Argonne National Laboratory. MINPACK Project. June 1983
    %     Jorge J. More', David J. Thuente
    %
    %     **********

  */
 
  int xtrapf = 4;
  int info = 0;
  int infoc = 1;

  //%     Compute the initial gradient in the search direction
  //%     and check that s is a descent direction.

  double dginit = in_prod(gr,sd); 
  if (dginit >= 0.)
    {
      printf("initial gradient in the search direction must be a descent\n");
      exit(1);
    }

  int brackt = 0;
  int stage1 = 1;
  int nfev = 0;
  double finit = f;
  double dgtest = ftol*dginit;
  double width = stpmax - stpmin;
  double width1 = 2*width;
  VEC *wa = v_get(x->dim);
  v_copy(x,wa);

  //%     The variables stx, fx, dgx contain the values of the step, 
  //%     function, and directional derivative at the best step.
  //%     The variables sty, fy, dgy contain the value of the step,
  //%     function, and derivative at the other endpoint of
  //%     the interval of uncertainty.
  //%     The variables stp, f, dg contain the values of the step,
  //%     function, and derivative at the current step.

  double stx = 0;
  double fx = finit;
  double dgx = dginit;
  double sty = 0;
  double fy = finit;
  double dgy = dginit;

  //%
  //%     Start of iteration.
  //%

  while (1) 
    {

      //%
      //%        Set the minimum and maximum steps to correspond
      //%        to the present interval of uncertainty.
      //%

      double stmin,stmax;

      if (brackt)
	{
	  stmin = min(stx,sty);
	  stmax = max(stx,sty);
	}
      else
	{
	  stmin = stx;
	  stmax = stp + xtrapf*(stp - stx);
	}

      //%
      //%        Force the step to be within the bounds stpmax and stpmin.
      //%

      stp = max(stp,stpmin);
      stp = min(stp,stpmax);

      //%
      //%        If an unusual termination is to occur then let 
      //%        stp be the lowest point obtained so far.
      //%

      if ((brackt && (stp <= stmin || stp >= stmax)) || nfev >= maxfev-1 || infoc == 0 || (brackt && stmax-stmin <= xtol*stmax))
	stp = stx;

      //%
      //%        Evaluate the function and gradient at stp
      //%        and compute the directional derivative.
      //%

      VEC *vtmp = v_get(x->dim);
      vtmp = sv_mlt(stp,sd,vtmp);
      v_add(wa,vtmp,x);
      V_FREE(vtmp);
      gr = (*fcn)(x,d,gr,&f);
      nfev += 1;
      double dg = in_prod(gr,sd);
      double ftest1 = finit + stp*dgtest;

      //%
      //%        Test for convergence.
      //%

      if ((brackt && (stp <= stmin || stp >= stmax)) || infoc == 0)
	info = 6;

      if (stp == stpmax && f <= ftest1 && dg <= dgtest)
	info = 5;

      if (stp == stpmin && (f > ftest1 || dg >= dgtest))
	info = 4;

      if (nfev >= maxfev)
	info = 3;

      if (brackt && stmax - stmin <= xtol*stmax)
	info = 2;

      if (f <= ftest1 && fabs(dg) <= gtol*(-dginit))
	info = 1;

      //%
      //%        Check for termination.
      //%

      if (info != 0){
	if (info != 1) {
	  printf("mthls exit condition: %d\n",info);
	  exit(1);
	}
	V_FREE(wa);
	return nfev;
      }

      //%
      //%        In the first stage we seek a step for which the modified
      //%        function has a nonpositive value and nonnegative derivative.
      //%

      if (stage1 && f <= ftest1 && dg >= min(ftol,gtol)*dginit)
	stage1 = 0;
      //%
      //%        A modified function is used to predict the step only if
      //%        we have not obtained a step for which the modified
      //%        function has a nonpositive function value and nonnegative 
      //%        derivative, and if a lower function value has been  
      //%        obtained but the decrease is not sufficient.
      //%

      if (stage1 && f <= fx && f > ftest1)
	{
	  //%
	  //%           Define the modified function and derivative values.
	  //%

	  double fm = f - stp*dgtest;
	  double fxm = fx - stx*dgtest;
	  double fym = fy - sty*dgtest;
	  double dgm = dg - dgtest;
	  double dgxm = dgx - dgtest;
	  double dgym = dgy - dgtest;
 
	  //% 
	  //%           Call cstep to update the interval of uncertainty 
	  //%           and to compute the new step.
	  //%

	  infoc = cstep(&stx,&fxm,&dgxm,&sty,&fym,&dgym,&stp,fm,dgm,&brackt,stmin,stmax);

	  //%
	  //%           Reset the function and gradient values for f.
	  //%

	  fx = fxm + stx*dgtest;
	  fy = fym + sty*dgtest;
	  dgx = dgxm + dgtest;
	  dgy = dgym + dgtest;

	}
      else
	{

	  //% 
	  //%           Call cstep to update the interval of uncertainty 
	  //%           and to compute the new step.
	  //%

	  infoc = cstep(&stx,&fx,&dgx,&sty,&fy,&dgy,&stp,f,dg,&brackt,stmin,stmax);
	}

      //
      //%        Force a sufficient decrease in the size of the
      //%        interval of uncertainty.
      //%

      if (brackt)
	{

	  if (fabs(sty-stx) >= .66*width1)
	    stp = stx + .5*(sty - stx);
	  width1 = width;
	  width = fabs(sty-stx);
	}

      //%
      //%        End of iteration.
      //%

    }

  //%
  //%     Last card of subroutine cvsrch.
  //%

}