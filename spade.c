// SPADE
// Stock assessment using PArtial Differential Equations
// Alex Campbell 'ghostofsandy' 2015 - 2016

// Testing:

// quick test:
//time ./spade -fn karumba .09 .09 1.3 .07
//number function evals: 42, function value: 790362.592819
//Vector: dim: 4
//   0.100412718    0.121894385     1.10724432    0.079491708
//real	1m1.321s
//user	2m30.212s
//sys	0m0.024s

// longer test:
// ./spade -fn karumba .8 .4 1.3 .1
// should return:
//  number function evals: 352, function value: 790362.592819
//  Vector: dim: 4
//     0.100412718    0.121894385     1.10724432    0.079

#include <fenv.h>
#include <math.h>
#include "spade.h"
#include "common.h"
#include "arg.h"
#include "parameters.h"
#include "machinery/VMGMM.h"
#include "machinery/objfns.h"
#include "optim/optim.h"
#include "plotting/plot.h"
#include "mathprop/mathprop.h"
#include "machinery/alpha/grad_alpha.h"
#include "machinery/beta/grad_beta.h"
#include "machinery/gamma/grad_gamma.h"
#include "machinery/iota/grad_iota.h"
#include "machinery/kappa/grad_kappa.h"
#include "machinery/omega/grad_omega.h"

int feenableexcept(int);

void print_usage() {
    printf(
      "SPADE: Stock assessment using PArtial Differential Equations.\n"
      "\n"
      "Usage:\n"
      "  spade -fn <file> -alpha <a> -beta <b> -gamma <g> -iota <i> -kappa <k>\n"
      "        -omega <w>\n"
      "\n"
      "Options:\n"
      "  -fn <file>\n"
      "      Specifies the common name of the input data files. SPADE will attempt\n"
      "      to read the files <file>-ce.dat and <file>-lf.dat from disk. This\n"
      "      option is required.\n"
      "\n"
      "  -minfish <minfish>      Default: 250\n"
      "\n"
      "  -J <J>                  Default: 400\n"
      "\n"
      "  -timestep <timestep>    Default: 0.025\n"
      "      The model timestep interval represented as a fraction of one year.\n"
      "\n"
      "  -warmup-ratio <ratio>   Default: 1\n"
      "      The SPADE model consists of two stages: a warmup stage followed by\n"
      "      a model stage. This option specifies the number of warmup steps as\n"
      "      a function of the number of model steps. For example, a warmup ratio\n"
      "      of 1 specified the same number of warmup steps as model steps. A warmup\n"
      "      ratio of 0.5 specifies half the number of warmup steps as compared to\n"
      "      model steps. Any numeric value greater than or equal to zero is an\n"
      "      acceptable warmup ratio.\n"
      "\n"
      "Parameters:\n"
      "  To disable a parameter suffix the parameter name with '-disabled'.\n"
      "  Example: 'spade -alpha-disabled <a> ...'\n"
  );
}

int main(int argc, char *argv[])
{
  feenableexcept(FE_DIVBYZERO); 
  feenableexcept(FE_INVALID); 
  feenableexcept(FE_OVERFLOW);

  int N;
  int minfish;
  Real k;
  Real warmup_ratio;

  // Read model-related command line args and set defaults if
  // arguments have not been provided
  if(arg_read_int("minfish", &minfish, argc, argv) == FALSE) {
    minfish = 250;
  }

  if(arg_read_int("J", &J, argc, argv) == FALSE) {
    J = 400;
  }

  if(arg_read_real("timestep", &k, argc, argv) == FALSE) {
    k = 0.025;
  }

  if(arg_read_real("warmup-ratio", &warmup_ratio, argc, argv) == FALSE) {
    // By default we have an equal number of warmup and model steps
    warmup_ratio = 1;
  }

  if(warmup_ratio < 0) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  // Read and parse data files
  char * data_file_name;
  if(arg_read_string("fn", &data_file_name, argc, argv) == FALSE) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  Data data;
  data_read_ce(data_file_name, &data, &N, k);
  data_read_lf(data_file_name, &data, N, k, minfish);

  // The model consists of two stages - a warmup stage followed by the model stage.
  // I (total number of time steps) = warmup_steps + N (number of time steps for model)
  //   * warmup_steps = N (number of time steps for model) * warmup_ratio
  //       For example a warmup_ratio of 1 would imply an equal number of warmup steps and model steps
  //   * warmup_ratio = number of warmup steps as a fraction of number of model steps
  //   * Number of time steps for model is determined by k (time step as fraction of a year) and  data->Y (number of years of data)

  int warmup_steps = floor(N * warmup_ratio);
  data.I = warmup_steps + N;
  data.S = warmup_steps;

  if(!SGNM) {
    data.J = J + data.I;
  } else {
    data.J = J;
  }

  data.k = k;

  // Read optim options
  OptimControl optim;
  optim_control_read("control.optim", &optim);

  // Configure each parameter. This must be updated when a
  // new parameter is created.
  Parameters parameters = {
   .alpha = { .name = "alpha", .grad = &grad_alpha },
   .beta = { .name = "beta", .grad = &grad_beta },
   .gamma = { .name = "gamma", .grad = &grad_gamma },
   .iota = { .name = "iota", .grad = &grad_iota },
   .kappa = { .name = "kappa", .grad = &grad_kappa },
   .omega = { .name = "omega", .grad = &grad_omega }
  };

  // Map all parameters to the parameter array. This must
  // be updated when a new parameter is created.
  parameters.parameter[0] = &parameters.alpha;
  parameters.parameter[1] = &parameters.beta;
  parameters.parameter[2] = &parameters.gamma;
  parameters.parameter[3] = &parameters.iota;
  parameters.parameter[4] = &parameters.kappa;
  parameters.parameter[5] = &parameters.omega;
  parameters.count = PARAMETER_COUNT;

  // Read all parameter values from the command line.
  if(!parameters_read(&parameters, argc, argv)) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  VEC *theta = parameters_to_vec(&parameters);

  h = parameters.omega.value / J;

  char labbuffer[10];
  sprintf(labbuffer,"before");
  //plot(&parameters,&data,labbuffer);
  
  theta = bfgs(VMGMM,theta,&data,&parameters,optim);

  char labbuffer2[10];
  sprintf(labbuffer2,"after");
  //plot(&parameters,&data,labbuffer2);
  
  V_FREE(theta);
  V_FREE(data.cat);
  V_FREE(data.eff);
  free(data.t_id);
  free(data.t_sz);

  for (int i=0;i<data.n;i++)
    free(data.lf[i]);
  free(data.lf);

  return(0);
}


/*

for splines

Li <- array(0,shorterN)
Lin <- array(0,shorterN)

min.lengths <- 20

for (i in 1:nrow(cssf))
    Lin[cssf$idx[i]] <- Lin[cssf$idx[i]] + 1

Li <- Lin > min.lengths

maxnk <- 0
maxnbasis <- 0

szs <- {}
timesidx <- {}

#outxall <- {}
#outyall <- {}


for (i in 1:shorterN) {
    if (Li[i]) {
      cand <- subset(cssf,idx==i)            
      
      cdi <- {}
      for (j in 1:nrow(cand)) 
          cdi <- c(cdi,rep(cand$Length1[j],100/cand$PctCatchSampled[j]))
      szs <- c(szs,length(cdi))
      timesidx <- c(timesidx,i)
      
      binidx <- .bincode(cdi,seq(0,200,by=3))

      m <- table(binidx)
      x <- (as.numeric(unlist(dimnames(m)))-1)*3+1.5
      y <- as.numeric(m) / sum(binidx)
     
      xl <- (min(x)-10):(min(x)-6)
      yl <- rep(0,5)
      xh <- (max(x)+6):(max(x)+10)
      yh <- rep(0,5)
      
      x <- c(xl,min(x)-3,x,max(x)+3,xh)
      y <- c(yl,y[1]/2,y,y[length(y)]/2,yh)            
      
      xsp <- (x - min(x)) / (max(x)-min(x))
      
      spl <- smooth.spline(xsp,y,nknots=length(xsp)-1,spar=.5)
      bspl.basis <- create.bspline.basis(unique(spl$fit$knot))

      #outx <- seq(0,1,length=300)
      #outy <- array(0,300)
      #for (j in 1:300)
      #    outy[j] <- eval.basis(outx[j],bspl.basis) %*% spl$fit$coef
           
      #if (i == 592){
      #    outx <- outx * (max(x)-min(x)) + min(x)
      #    outxall <- outx
      #    outy <- outy*200
      #    outyall <- outy
      #    #plot(outx,outy,type='l',ylim=c(-100,5),xlim=c(50,180))
      #}
      #else if (i>592) {
      #    outx <- outx * (max(x)-min(x)) + min(x)
      #    outxall <- c(outxall,outx)
      #    outy <- (outy*200) - epsilon*(1-exp(-delta*(i-592)/25))
      #    outyall <- c(outyall,outy)
      #    #lines(outx,outy)
      #}
      
      nbasis <- bspl.basis$nbasis
      params <- bspl.basis$params
      rangeval <- bspl.basis$rangeval
      breaks <- c(rangeval[1], params, rangeval[2])
      norder <- nbasis - length(breaks) + 2
      nbreaks <- length(breaks)
      knots <- c(rep(breaks[1], norder -1), breaks, rep(breaks[nbreaks], norder -1))
      nbasis <- nbreaks + norder - 2
      nk <- length(knots)

      if (nbasis>maxnbasis)
          maxnbasis <- nbasis
      if (nk>maxnk)
          maxnk <- nk

    }
}


ntimes <- length(szs)
splcoef <- matrix(0,ntimes,maxnbasis)
knots <- matrix(0,ntimes,maxnk)
nbases <- array(0,ntimes)
nks <- array(0,ntimes)
timesidx <- array(0,ntimes)
szs <- array(0,ntimes)
starts <- array(0,ntimes)
stops <- array(0,ntimes)

jj <- 1

#pdf('lplot.pdf')

for (i in 1:shorterN) {
    if (Li[i]) {
      cand <- subset(cssf,idx==i)
      cdi <- {}
      for (j in 1:nrow(cand)) 
          cdi <- c(cdi,rep(cand$Length1[j],100/cand$PctCatchSampled[j]))

      szs[jj] <- length(cdi)
      timesidx[jj] <- i
      
      binidx <- .bincode(cdi,seq(0,200,by=3))
      
      m <- table(binidx)
      x <- (as.numeric(unlist(dimnames(m)))-1)*3+1.5
      y <- as.numeric(m) / sum(binidx)
                
      xl <- (min(x)-10):(min(x)-6)
      yl <- rep(0,5)
      xh <- (max(x)+6):(max(x)+10)
      yh <- rep(0,5)
      
      x <- c(xl,min(x)-3,x,max(x)+3,xh)
      y <- c(yl,y[1]/2,y,y[length(y)]/2,yh)
      
      starts[jj] <- min(x)
      stops[jj] <- max(x)
      
      xsp <- (x - min(x)) / (max(x)-min(x))
            
      spl <- smooth.spline(xsp,y,nknots=length(xsp)-1,spar=.5)
      bspl.basis <- create.bspline.basis(unique(spl$fit$knot))

            
      #outx <- seq(0,1,length=300)
      #outy <- array(0,300)
      #for (j in 1:300)
      #    outy[j] <- eval.basis(outx[j],bspl.basis) %*% spl$fit$coef
      
      #if (i == 592){
      #    outx <- outx * (max(x)-min(x)) + min(x)      
      #    outy <- outy*200
      #    plot(outx,outy,type='l',xlim=range(outxall),ylim=range(outyall))
      #}
      #else if(i > 592) {
      #    outx <- outx * (max(x)-min(x)) + min(x)
      #    outy <- (outy*200) - epsilon*(1-exp(-delta*(i-592)/25))
      #    lines(outx,outy)
      #}            
      
      
      nbasis <- bspl.basis$nbasis
      params <- bspl.basis$params
      rangeval <- bspl.basis$rangeval
      breaks <- c(rangeval[1], params, rangeval[2])
      norder <- nbasis - length(breaks) + 2
      nbreaks <- length(breaks)
      knotstmp <- c(rep(breaks[1], norder -1), breaks, rep(breaks[nbreaks], norder -1))
      knots[jj,1:length(knotstmp)] <- knotstmp
      nbasis <- nbreaks + norder - 2
      nk <- length(knotstmp)
      splcoeftmp <- spl$fit$coef
      splcoef[jj,1:length(splcoeftmp)] <- splcoeftmp
            
      nbases[jj] <- nbasis
      nks[jj] <- nk
     
      jj <- jj + 1
    }
}

#dev.off()

subst <- c(5,37)
#subst <- c(5,11,26,37)
#subst <- c(7,10,26,40)
#subst <- c(1,2,3,4)
ntimes <- length(subst)

nks <- nks[subst]
nbases <- nbases[subst]
knots <- knots[subst,]
splcoef <- splcoef[subst,]
szs <- szs[subst]
timesidx <- timesidx[subst]
starts <- starts[subst]
stops <- stops[subst]


{
    sink('barra_agg.dat',append=TRUE)
    cat(sprintf("%d\n",ntimes))
    for (i in 1:ntimes)
      cat(sprintf("%d ",nks[i]))
    cat(sprintf("\n"))
    for (i in 1:ntimes)
      cat(sprintf("%d ",nbases[i]))
    cat(sprintf("\n"))
    for (i in 1:ntimes)
      cat(sprintf("%d ",1))
    cat(sprintf("\n"))
    for (i in 1:ntimes) {
        for (j in 1:nks[i])  
            cat(sprintf("%f ",knots[i,j]))
        cat(sprintf("\n"))
    }
    for (i in 1:ntimes) {
        for (j in 1:nbases[i])  
            cat(sprintf("%f ",splcoef[i,j]))
        cat(sprintf("\n"))
    }
    for (i in 1:ntimes)      
      cat(sprintf("%d ",szs[i]))
    cat(sprintf("\n"))
    for (i in 1:ntimes)
      cat(sprintf("%d ",timesidx[i]))
    cat(sprintf("\n"))
    for (i in 1:ntimes)      
      cat(sprintf("%f ",starts[i]))
    cat(sprintf("\n"))
    for (i in 1:ntimes)
      cat(sprintf("%f ",stops[i]))
    cat(sprintf("\n"))
    sink()
}


...

in tpl file:

  init_int T1
  init_ivector Kn1(1,T1)
  init_ivector B1(1,T1)
  init_ivector one1(1,T1)
  init_matrix knots1(1,T1,one1,Kn1)
  init_matrix coef1(1,T1,one1,B1)
  init_vector start1(1,T1)
  init_vector stop1(1,T1)


  lf1=0;
  p1.initialize();
  sn.initialize();
  int q,st,jj;
  dvariable ans2;
  for (int q=1;q<=T1;q++) {
    i = H1(q) + Hw1/2;

    for (j=0;j<=J;j++)
      if (x(i,j) > start1(q))
        break;
    st = j;
    
    for (jj=st;jj<=J;jj++) {

      if (x(i,jj) > stop1(q))
        break;

      ans2 = sple(Kn1(q),knots1(q),(x(i,jj)-start1(q))/(stop1(q)-start1(q)),coef1(q));            
      if (ans2 > 0)
        p1(jj) = ans2;
      else 
        p1(jj) = 1e-3;            

      sn(jj) = s(x(i,jj))*n1(i,jj);

    }

    dvar_vector xtmp = extract_row(x,i);
    dvar_vector xtmp2 = xtmp(st,jj-1);
    dvar_vector p12 = p1(st,jj-1);
    dvar_vector sn2 = sn(st,jj-1);
          
    p12.shift(0);
    sn2.shift(0);
    xtmp2.shift(0);

    p12 *= 1/Qjn(xtmp2,p12,xtmp2.size()-1);    
    sn2 *= 1/Qjn(xtmp2,sn2,xtmp2.size()-1);
         
    lf1 += Qjn(xtmp2,elem_prod(p12,log(elem_div(p12,sn2))),xtmp2.size()-1);

  } 



..

FUNCTION dvariable sple(const int nk, const dvar_vector& knots, const dvariable& xval, const dvar_vector& coef)

  int j,q,r,l,offset;
  dvar_vector val(1,4);
  dvar_vector rdel(1,3);
  dvar_vector ldel(1,3); 

  for (j=1;j<=nk;j++)
    if (knots(j) >= xval)
      break;

  l = j - 4;
  offset = l-1;

  for (q=0;q<=(4-2);q++) {
    rdel(q+1) = knots(j+q) - xval;
    ldel(q+1) = xval - knots(j - (q+1));
  }

  val(1) = 1;
  dvariable saved;
  dvariable term;
  for (q=1;q<=(4-1);q++) {
    saved=0;
    for (r=0;r<=(q-1);r++) {
      term = val(r+1) / (rdel(r+1) + ldel(q-1-r+1));
      val(r+1) = saved + rdel(r+1) * term;
      saved = ldel(q-1-r+1) * term;
    }
    val(q+1) = saved;
  }

  int ncoef = nk - 4;

  dvar_vector design(1,ncoef);
  design.initialize();
  design(offset+1) = val(1);
  design(offset+2) = val(2);
  design(offset+3) = val(3);
  design(offset+4) = val(4);

  return design * coef;


*/
