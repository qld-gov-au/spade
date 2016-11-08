#include <math.h>
#include "../meschach/matrix.h"
#include "../common.h"
#include "../parameters.h"
#include "spade_solve.h"
#include "spade_solve_clean.h"
#include "../model/fishing/selectivity.h"
#include "../model/biology/weight.h"
#include "../model/fishing/effort.h"
#include "../model/fishing/catch.h"
#include "../model/biology/birth.h"
#include "Q.h"
#include "../util/util.h"
#include <time.h>

#include <immintrin.h>

# define ALIGN32_BEG
# define ALIGN32_END __attribute__((aligned(32)))

/* __m128 is ugly to write */
typedef __m256  v8sf; // vector of 8 float (avx)
typedef __m256i v8si; // vector of 8 int   (avx)
typedef __m128i v4si; // vector of 8 int   (avx)

#define _PI32AVX_CONST(Name, Val)                                            \
  static const ALIGN32_BEG int _pi32avx_##Name[4] ALIGN32_END = { Val, Val, Val, Val }

_PI32AVX_CONST(1, 1);
_PI32AVX_CONST(inv1, ~1);
_PI32AVX_CONST(2, 2);
_PI32AVX_CONST(4, 4);


/* declare some AVX constants -- why can't I figure a better way to do that? */
#define _PS256_CONST(Name, Val)                                            \
  static const ALIGN32_BEG float _ps256_##Name[8] ALIGN32_END = { Val, Val, Val, Val, Val, Val, Val, Val }
#define _PI32_CONST256(Name, Val)                                            \
  static const ALIGN32_BEG int _pi32_256_##Name[8] ALIGN32_END = { Val, Val, Val, Val, Val, Val, Val, Val }
#define _PS256_CONST_TYPE(Name, Type, Val)                                 \
  static const ALIGN32_BEG Type _ps256_##Name[8] ALIGN32_END = { Val, Val, Val, Val, Val, Val, Val, Val }

_PS256_CONST(1  , 1.0f);
_PS256_CONST(0p5, 0.5f);
/* the smallest non denormalized float number */
_PS256_CONST_TYPE(min_norm_pos, int, 0x00800000);
_PS256_CONST_TYPE(mant_mask, int, 0x7f800000);
_PS256_CONST_TYPE(inv_mant_mask, int, ~0x7f800000);

_PS256_CONST_TYPE(sign_mask, int, 0x80000000);
_PS256_CONST_TYPE(inv_sign_mask, int, ~0x80000000);

_PI32_CONST256(0, 0);
_PI32_CONST256(1, 1);
_PI32_CONST256(inv1, ~1);
_PI32_CONST256(2, 2);
_PI32_CONST256(4, 4);
_PI32_CONST256(0x7f, 0x7f);

_PS256_CONST(cephes_SQRTHF, 0.707106781186547524);
_PS256_CONST(cephes_log_p0, 7.0376836292E-2);
_PS256_CONST(cephes_log_p1, - 1.1514610310E-1);
_PS256_CONST(cephes_log_p2, 1.1676998740E-1);
_PS256_CONST(cephes_log_p3, - 1.2420140846E-1);
_PS256_CONST(cephes_log_p4, + 1.4249322787E-1);
_PS256_CONST(cephes_log_p5, - 1.6668057665E-1);
_PS256_CONST(cephes_log_p6, + 2.0000714765E-1);
_PS256_CONST(cephes_log_p7, - 2.4999993993E-1);
_PS256_CONST(cephes_log_p8, + 3.3333331174E-1);
_PS256_CONST(cephes_log_q1, -2.12194440e-4);
_PS256_CONST(cephes_log_q2, 0.693359375);

#ifndef __AVX2__

typedef union imm_xmm_union {
  v8si imm;
  v4si xmm[2];
} imm_xmm_union;

#define COPY_IMM_TO_XMM(imm_, xmm0_, xmm1_) {    \
    imm_xmm_union u __attribute__((aligned(32)));  \
    u.imm = imm_;				   \
    xmm0_ = u.xmm[0];                            \
    xmm1_ = u.xmm[1];                            \
}

#define COPY_XMM_TO_IMM(xmm0_, xmm1_, imm_) {                       \
    imm_xmm_union u __attribute__((aligned(32))); \
    u.xmm[0]=xmm0_; u.xmm[1]=xmm1_; imm_ = u.imm; \
  }


#define AVX2_BITOP_USING_SSE2(fn) \
static inline v8si _mm256_##fn(v8si x, int a) \
{ \
  /* use SSE2 instruction to perform the bitop AVX2 */ \
  v4si x1, x2; \
  v8si ret; \
  COPY_IMM_TO_XMM(x, x1, x2); \
  x1 = _mm_##fn(x1,a); \
  x2 = _mm_##fn(x2,a); \
  COPY_XMM_TO_IMM(x1, x2, ret); \
  return(ret); \
}

#warning "Using SSE2 to perform AVX2 bitshift ops"
AVX2_BITOP_USING_SSE2(slli_epi32)
AVX2_BITOP_USING_SSE2(srli_epi32)

#define AVX2_INTOP_USING_SSE2(fn) \
static inline v8si _mm256_##fn(v8si x, v8si y) \
{ \
  /* use SSE2 instructions to perform the AVX2 integer operation */ \
  v4si x1, x2; \
  v4si y1, y2; \
  v8si ret; \
  COPY_IMM_TO_XMM(x, x1, x2); \
  COPY_IMM_TO_XMM(y, y1, y2); \
  x1 = _mm_##fn(x1,y1); \
  x2 = _mm_##fn(x2,y2); \
  COPY_XMM_TO_IMM(x1, x2, ret); \
  return(ret); \
}

#warning "Using SSE2 to perform AVX2 integer ops"
AVX2_INTOP_USING_SSE2(and_si128)
AVX2_INTOP_USING_SSE2(andnot_si128)
AVX2_INTOP_USING_SSE2(cmpeq_epi32)
AVX2_INTOP_USING_SSE2(sub_epi32)
AVX2_INTOP_USING_SSE2(add_epi32)

#endif /* __AVX2__ */



/* natural logarithm computed for 8 simultaneous float 
   return NaN for x <= 0
*/
v8sf log256_ps(v8sf x) {
  v8si imm0;
  v8sf one = *(v8sf*)_ps256_1;

  //v8sf invalid_mask = _mm256_cmple_ps(x, _mm256_setzero_ps());
  v8sf invalid_mask = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LE_OS);

  x = _mm256_max_ps(x, *(v8sf*)_ps256_min_norm_pos);  /* cut off denormalized stuff */

  // can be done with AVX2
  imm0 = _mm256_srli_epi32(_mm256_castps_si256(x), 23);

  /* keep only the fractional part */
  x = _mm256_and_ps(x, *(v8sf*)_ps256_inv_mant_mask);
  x = _mm256_or_ps(x, *(v8sf*)_ps256_0p5);

  // this is again another AVX2 instruction
  imm0 = _mm256_sub_epi32(imm0, *(v8si*)_pi32_256_0x7f);
  v8sf e = _mm256_cvtepi32_ps(imm0);

  e = _mm256_add_ps(e, one);

  /* part2: 
     if( x < SQRTHF ) {
       e -= 1;
       x = x + x - 1.0;
     } else { x = x - 1.0; }
  */
  //v8sf mask = _mm256_cmplt_ps(x, *(v8sf*)_ps256_cephes_SQRTHF);
  v8sf mask = _mm256_cmp_ps(x, *(v8sf*)_ps256_cephes_SQRTHF, _CMP_LT_OS);
  v8sf tmp = _mm256_and_ps(x, mask);
  x = _mm256_sub_ps(x, one);
  e = _mm256_sub_ps(e, _mm256_and_ps(one, mask));
  x = _mm256_add_ps(x, tmp);

  v8sf z = _mm256_mul_ps(x,x);

  v8sf y = *(v8sf*)_ps256_cephes_log_p0;
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p1);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p2);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p3);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p4);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p5);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p6);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p7);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_log_p8);
  y = _mm256_mul_ps(y, x);

  y = _mm256_mul_ps(y, z);
  
  tmp = _mm256_mul_ps(e, *(v8sf*)_ps256_cephes_log_q1);
  y = _mm256_add_ps(y, tmp);


  tmp = _mm256_mul_ps(z, *(v8sf*)_ps256_0p5);
  y = _mm256_sub_ps(y, tmp);

  tmp = _mm256_mul_ps(e, *(v8sf*)_ps256_cephes_log_q2);
  x = _mm256_add_ps(x, y);
  x = _mm256_add_ps(x, tmp);
  x = _mm256_or_ps(x, invalid_mask); // negative arg will be NAN
  return x;
}



_PS256_CONST(exp_hi,	88.3762626647949f);
_PS256_CONST(exp_lo,	-88.3762626647949f);

_PS256_CONST(cephes_LOG2EF, 1.44269504088896341);
_PS256_CONST(cephes_exp_C1, 0.693359375);
_PS256_CONST(cephes_exp_C2, -2.12194440e-4);

_PS256_CONST(cephes_exp_p0, 1.9875691500E-4);
_PS256_CONST(cephes_exp_p1, 1.3981999507E-3);
_PS256_CONST(cephes_exp_p2, 8.3334519073E-3);
_PS256_CONST(cephes_exp_p3, 4.1665795894E-2);
_PS256_CONST(cephes_exp_p4, 1.6666665459E-1);
_PS256_CONST(cephes_exp_p5, 5.0000001201E-1);

v8sf exp256_ps(v8sf x) {
  v8sf tmp = _mm256_setzero_ps(), fx;
  v8si imm0;
  v8sf one = *(v8sf*)_ps256_1;

  x = _mm256_min_ps(x, *(v8sf*)_ps256_exp_hi);
  x = _mm256_max_ps(x, *(v8sf*)_ps256_exp_lo);

  /* express exp(x) as exp(g + n*log(2)) */
  fx = _mm256_mul_ps(x, *(v8sf*)_ps256_cephes_LOG2EF);
  fx = _mm256_add_ps(fx, *(v8sf*)_ps256_0p5);

  /* how to perform a floorf with SSE: just below */
  //imm0 = _mm256_cvttps_epi32(fx);
  //tmp  = _mm256_cvtepi32_ps(imm0);
  
  tmp = _mm256_floor_ps(fx);

  /* if greater, substract 1 */
  //v8sf mask = _mm256_cmpgt_ps(tmp, fx);    
  v8sf mask = _mm256_cmp_ps(tmp, fx, _CMP_GT_OS);    
  mask = _mm256_and_ps(mask, one);
  fx = _mm256_sub_ps(tmp, mask);

  tmp = _mm256_mul_ps(fx, *(v8sf*)_ps256_cephes_exp_C1);
  v8sf z = _mm256_mul_ps(fx, *(v8sf*)_ps256_cephes_exp_C2);
  x = _mm256_sub_ps(x, tmp);
  x = _mm256_sub_ps(x, z);

  z = _mm256_mul_ps(x,x);
  
  v8sf y = *(v8sf*)_ps256_cephes_exp_p0;
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_exp_p1);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_exp_p2);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_exp_p3);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_exp_p4);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, *(v8sf*)_ps256_cephes_exp_p5);
  y = _mm256_mul_ps(y, z);
  y = _mm256_add_ps(y, x);
  y = _mm256_add_ps(y, one);

  /* build 2^n */
  imm0 = _mm256_cvttps_epi32(fx);
  // another two AVX2 instructions
  imm0 = _mm256_add_epi32(imm0, *(v8si*)_pi32_256_0x7f);
  imm0 = _mm256_slli_epi32(imm0, 23);
  v8sf pow2n = _mm256_castsi256_ps(imm0);
  y = _mm256_mul_ps(y, pow2n);
  return y;
}


float hsum_ps_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums        = _mm_add_ss(sums, shuf);
    return        _mm_cvtss_f32(sums);
}

float hsum256_ps_avx(__m256 v) {
    __m128 vlow  = _mm256_castps256_ps128(v);
    __m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
           vlow  = _mm_add_ps(vlow, vhigh);     // add the low 128
    return hsum_ps_sse3(vlow);         // and inline the sse3 version, which is optimal for AVX
    // (no wasted instructions, and all of them are the 4B minimum)
}

void Kveryfast(

	   void *args	 
	  )
{

  Parameters *parameters = (Parameters *)args;
  
  Real a1 = parameters->alpha1.value;
  Real a2 = parameters->alpha2.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;

  int J = d.J+1;
  
  float *x;
  posix_memalign(&x,32,J*sizeof(float));

  __m256 *_x;

  posix_memalign(&_x,32,(J/8)*sizeof(*_x));

  __m256 *_xh;

  posix_memalign(&_xh,32,(J/8+1)*sizeof(*_xh));

  float *xh;
  posix_memalign(&xh,32,J*sizeof(float));

  __m256 *_xn;

  posix_memalign(&_xn,32,(J/8+1)*sizeof(*_xh));

  float *xn;
  posix_memalign(&xn,32,J*sizeof(float));
  
  float * u;
  posix_memalign(&u,32,J*sizeof(float));

  float * r;
  posix_memalign(&r,32,J*sizeof(float));

  float * o;
  posix_memalign(&o,32,J*sizeof(float));
  
  __m256 *_u;

  posix_memalign(&_u,32,(J/8)*sizeof(*_u));

  __m256 *_uh;

  posix_memalign(&_uh,32,(J/8+1)*sizeof(*_uh));

  float *uh;
  posix_memalign(&uh,32,J*sizeof(float));

  __m256 *_un;

  posix_memalign(&_un,32,(J/8+1)*sizeof(*_uh));

  float *un;
  posix_memalign(&un,32,J*sizeof(float));

  __m256 *_r;

  posix_memalign(&_r,32,(J/8)*sizeof(*_r));
  
  __m256i vindex = _mm256_set_epi32(32,28,24,20,16,12,8,4);
  
       
  /* 
     initialize x
  */  
  for (int j=0;j<J;j++) {
    x[j] = ww*(float)j/J;
  }

  float *px = x;

  int it=0;
  for (int j=0;j<J;j+=8)
    {
      _x[it] = _mm256_load_ps(px);
      px += 8;
      it++;
    }

  
  /*
     ok, now initialize u
  */
  
  // prelims
  Real zeta = sqrt( 81*kk*kk*ww*ww*pow(a1+2*a2*ww,2.) - 12*kk*pow(a1*ww+kk,3.) );
  Real eta = 9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta;
  Real Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*kk*(a1*ww+kk) / pow(eta,1./3);

  Real ubar = (Z - bb - kk) / gg; 
  Real vbar = (kk*ww*ubar) / (bb+gg*ubar+kk);
  Real wbar = (2*kk*ww*vbar) / (bb+gg*ubar+2*kk);  
  
  // set
  for (int j=0;j<J;j++)
    u[j] = (a1*vbar+a2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));
  
  float *pu = u;

  it=0;
  for (int j=0;j<J;j+=8)
    {
      _u[it] = _mm256_load_ps(pu);
      pu += 8;
      it++;
    }

  __m256 const_ww = _mm256_set1_ps((float)ww);
  __m256 const_k2kk = _mm256_set1_ps((float)(k/2)*kk);
  __m256 const_phi_iota1 = _mm256_set1_ps((float)(phi*iota1));
  __m256 const_2iota2 = _mm256_set1_ps((float)(2*iota2));
  __m256 const_phisq = _mm256_set1_ps((float)(phi*phi));
  __m256 const_neg1 = _mm256_set1_ps(-1.0f);
  __m256 const_nk2 = _mm256_set1_ps((float)(-(k/2)));
  
  clock_t begin = clock();
  
  int i=1;
  float t = k*(i-1);
  it=0;
  for (int j=0;j<J;j+=8)
    {
    
      __m256 xtmp1 = _mm256_sub_ps(const_ww,_x[it]);
      
      xtmp1 = _mm256_mul_ps(const_k2kk,xtmp1);
      _xh[it] = _mm256_add_ps(_x[it],xtmp1);      
      
      __m256 tmp = _mm256_sub_ps(_x[it],const_phi_iota1);

      tmp = _mm256_mul_ps(tmp,tmp);
      tmp = _mm256_mul_ps(const_neg1,tmp);
      
      __m256 tmp2 = _mm256_mul_ps(const_2iota2,const_phisq);
      tmp = _mm256_div_ps(tmp,tmp2);

      tmp = exp256_ps(tmp);

      __m256 const_iie = _mm256_set1_ps((float)(ii*_e(d.eff,t)));

      tmp = _mm256_mul_ps(const_iie,tmp);
      
      __m256 const_bbpggUmkk = _mm256_set1_ps((float)(bb+gg*ubar-kk));

      tmp = _mm256_add_ps(const_bbpggUmkk,tmp);
      
      tmp = _mm256_mul_ps(const_nk2,tmp);

      tmp = exp256_ps(tmp);

      _uh[it] = _mm256_mul_ps(_u[it],tmp);

      _mm256_store_ps(&(xh[j]),_xh[it]);
      _mm256_store_ps(&(uh[j]),_uh[it]);
      
      it++;
      
    }
  
  float uh_0 = xh[0] * _b(a1,a2,xh[0])*uh[0];

  _xh[J/8] = _mm256_set1_ps(ww);  // zero past the end of the array
  _uh[J/8] = _mm256_setzero_ps();  // zero past the end of the array
  
  __m256 const_a1 = _mm256_set1_ps((float)a1);
  __m256 const_a2 = _mm256_set1_ps((float)a2);
    
  it=0;
  for (int j=0;j<J;j+=8)
    {
      
      __m256 xh_shift = _mm256_i32gather_ps(&(_xh[it]),vindex,1);
            
      __m256 x2_1 = _mm256_mul_ps(_xh[it],_xh[it]);
      __m256 tmpb_1 = _mm256_mul_ps(const_a2,x2_1);

      __m256 tmpb_i = _mm256_mul_ps(const_a1,_xh[it]);
      
      tmpb_1 = _mm256_add_ps(tmpb_i,tmpb_1);

      __m256 pt1 = _mm256_mul_ps(tmpb_1,_uh[it]);

      __m256 x2_2 = _mm256_mul_ps(xh_shift,xh_shift);
      __m256 tmpb_2 = _mm256_mul_ps(const_a2,x2_2);

      tmpb_i = _mm256_mul_ps(const_a1,xh_shift);
      
      tmpb_2 = _mm256_add_ps(tmpb_i,tmpb_2);

      __m256 uh_shift = _mm256_i32gather_ps(&(_uh[it]),vindex,1);      
      
      __m256 pt2 = _mm256_mul_ps(tmpb_2,uh_shift);

      __m256 sm = _mm256_add_ps(pt1,pt2);

      __m256 subx = _mm256_sub_ps(xh_shift,_xh[it]);
      
      __m256 ans = _mm256_mul_ps(sm,subx);

      //      _mm256_store_ps(&(uh[j]),ans);
      
      uh_0 += hsum256_ps_avx(ans);
      
      it++;
    }
  			 
  uh_0 /= (2*kk*ww - xh[0]*_b(a1,a2,0));

  float Uh = .5 * xh[0] * (uh_0 + uh[0]);
  for (int j=0;j<J-1;j++)
    Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);

  float th = k*(i-.5);

  __m256 const_kkk = _mm256_set1_ps((float)k*kk);
  __m256 const_nk = _mm256_set1_ps((float)-k);
  __m256 const_iieth = _mm256_set1_ps((float)(ii*_e(d.eff,th)));
  __m256 const_bbpggUhmkk = _mm256_set1_ps((float)(bb+gg*Uh-kk));
  
  it=0;
  for (int j=0;j<J;j+=8)
    {
      
      __m256 xtmp1 = _mm256_sub_ps(const_ww,_xh[it]);
      
      xtmp1 = _mm256_mul_ps(const_kkk,xtmp1);
      _xn[it] = _mm256_add_ps(_x[it],xtmp1);      
      
      __m256 tmp = _mm256_sub_ps(_xh[it],const_phi_iota1);

      tmp = _mm256_mul_ps(tmp,tmp);
      tmp = _mm256_mul_ps(const_neg1,tmp);
      
      __m256 tmp2 = _mm256_mul_ps(const_2iota2,const_phisq);
      tmp = _mm256_div_ps(tmp,tmp2);

      tmp = exp256_ps(tmp);

      tmp = _mm256_mul_ps(const_iieth,tmp);
      
      tmp = _mm256_add_ps(const_bbpggUhmkk,tmp);
      
      tmp = _mm256_mul_ps(const_nk,tmp);

      tmp = exp256_ps(tmp);

      _un[it] = _mm256_mul_ps(_u[it],tmp);

      _mm256_store_ps(&(xn[j]),_xn[it]);
      _mm256_store_ps(&(un[j]),_un[it]);
      
      it++;
      
    }

  float un_0 = xn[0] * _b(a1,a2,xn[0])*un[0]; //      un[0] = x[1] * _b(a1,a2,xn[1])*un[1];

  _xn[J/8] = _mm256_set1_ps(ww);  // zero past the end of the array
  _un[J/8] = _mm256_setzero_ps();  // zero past the end of the array
        
  it=0;
  for (int j=0;j<J;j+=8)
    {

      __m256 xn_shift = _mm256_i32gather_ps(&(_xn[it]),vindex,1);
            
      __m256 x2_1 = _mm256_mul_ps(_xn[it],_xn[it]);
      __m256 tmpb_1 = _mm256_mul_ps(const_a2,x2_1);

      __m256 tmpb_i = _mm256_mul_ps(const_a1,_xn[it]);
      
      tmpb_1 = _mm256_add_ps(tmpb_i,tmpb_1);

      __m256 pt1 = _mm256_mul_ps(tmpb_1,_un[it]);

      __m256 x2_2 = _mm256_mul_ps(xn_shift,xn_shift);
      __m256 tmpb_2 = _mm256_mul_ps(const_a2,x2_2);

      tmpb_i = _mm256_mul_ps(const_a1,xn_shift);
      
      tmpb_2 = _mm256_add_ps(tmpb_i,tmpb_2);

      __m256 un_shift = _mm256_i32gather_ps(&(_un[it]),vindex,1);      
      
      __m256 pt2 = _mm256_mul_ps(tmpb_2,un_shift);

      __m256 sm = _mm256_add_ps(pt1,pt2);

      __m256 subx = _mm256_sub_ps(xn_shift,_xn[it]);
      
      __m256 ans = _mm256_mul_ps(sm,subx);
      
      un_0 += hsum256_ps_avx(ans);
      
      it++;
    }

  un_0 /= (2*kk*ww - xn[0]*_b(a1,a2,0)); //    un[0] /= (2*kk*ww - xn[1]*_b(a1,a2,0));

  __m256 const_1en3eta1 = _mm256_set1_ps((float)(1e-3*eta1));
  __m256 const_eta2 = _mm256_set1_ps((float)eta2);
  __m256 const_iiet1 = _mm256_set1_ps((float)(ii*_e(d.eff,k*i)));

  it=0;
  for (int j=0;j<J;j+=8)
    {

      __m256 tmpw = _mm256_mul_ps(const_1en3eta1,exp256_ps(_mm256_mul_ps(const_eta2, log256_ps(_xn[it]))));
      
      __m256 tmps = _mm256_sub_ps(_xn[it],const_phi_iota1);

      tmps = _mm256_mul_ps(tmps,tmps);
      tmps = _mm256_mul_ps(const_neg1,tmps);
      
      __m256 tmp2 = _mm256_mul_ps(const_2iota2,const_phisq);
      tmps = _mm256_div_ps(tmps,tmp2);

      tmps = exp256_ps(tmps);

      __m256 tmp = _mm256_mul_ps(tmpw,tmps);
      tmp = _mm256_mul_ps(tmp,_un[it]);
      
      _r[it] = _mm256_mul_ps(tmp,const_iiet1);
      
      _mm256_store_ps(&(r[j]),_r[it]);

      it++;
    }

  float U = .5 * xn[0] * (un[0] + un_0);
  float C = 0;
  for (int j=0;j<J-1;j++)
    {
      U += .5 * (xn[j+1] - xn[j]) * (un[j+1] + un[j]);
      C += .5 * (xn[j+1] - xn[j]) * (r[j+1] + r[j]);
    }

  __m256 const_ct = _mm256_set1_ps((float)(_c(d.cat,k*i)));
  __m256 const_1mqp = _mm256_set1_ps((float)(1-d.Qp[i]));
  __m256 const_C = _mm256_set1_ps(C);

  float *pp = d.p[i];
  
  it=0;
  for (int j=0;j<J;j+=8)
    {

      __m256 _p = _mm256_load_ps(pp);

      __m256 tmp = _mm256_mul_ps(const_1mqp,_mm256_div_ps(_r[it],const_C));

      tmp = _mm256_add_ps(_p,tmp);

      tmp = _mm256_mul_ps(const_ct,tmp);

      __m256 ans = _mm256_sub_ps(tmp,_r[it]);

      ans = _mm256_mul_ps(ans,ans);

      _mm256_store_ps(&(o[j]),ans);

      pp+=8;
      it++;

      }
      
  
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("time: %lf\n",time_spent);
  
  printf("\n%lf\n",C);

  //float *xhtest;
  //posix_memalign(&xhtest,32,J*sizeof(float));

  //float *uhtest;
  //posix_memalign(&uhtest,32,J*sizeof(float));


  // recompute because corrupted memory??
  
  for (int j=0;j<J;j++) {
    x[j] = ww*(float)j/J;
  }

  for (int j=0;j<J;j++)
    u[j] = (a1*vbar+a2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));
  
  float *xhtest = (float *)calloc(J,sizeof(float));
  float *uhtest = (float *)calloc(J,sizeof(float));  
  
  begin = clock();
  
  for (int j=0;j<J;j++)
    {
      xhtest[j] = x[j] + k/2 * kk*(ww - x[j]);    
      uhtest[j] = u[j] * exp( -k/2 * (bb + gg*ubar + s(x[j])* ii * _e(d.eff,t) - kk) );
    }

  float uh_0test = xhtest[0] * _b(a1,a2,xhtest[0])*uhtest[0];
    
  for (int j=0;j<J-1;j++)
    uh_0test += (_b(a1,a2,xhtest[j])*uhtest[j] + _b(a1,a2,xhtest[j+1])*uhtest[j+1]) * (xhtest[j+1]-xhtest[j]);

  uh_0test /= (2*kk*ww - xhtest[0]*_b(a1,a2,0));

  float Uhtest = .5 * xhtest[0] * (uh_0test + uhtest[0]);
  
  for (int j=0;j<J-1;j++)
    Uhtest += .5 * (xhtest[j+1] - xhtest[j]) * (uhtest[j+1] + uhtest[j]);

  float *xntest = (float *)calloc(J+1,sizeof(float));
  float *untest = (float *)calloc(J+1,sizeof(float));  

  for (int j=J+1;j>0;j--)
    {
      xntest[j] = x[j-1] + k * kk*(ww-xhtest[j-1]);
      untest[j] = u[j-1] * exp ( -k * (bb + gg*Uh + s(xhtest[j-1]) * ii * _e(d.eff,th) - kk));
    }
  
  xntest[0] = 0;
  untest[0] = x[1] * _b(a1,a2,xntest[1])*untest[1];
      
  for (int j=1;j<=J;j++)
    untest[0] += (_b(a1,a2,xntest[j])*untest[j] + _b(a1,a2,xntest[j+1])*untest[j+1]) * (xntest[j+1]-xntest[j]);

  untest[0] /= (2*kk*ww - xntest[1]*_b(a1,a2,0));

  float *rtest = (float *)calloc(J+1,sizeof(float));  

  for (int j=0;j<=J;j++)
    rtest[j] = 1e-3*w(xntest[j])*s(xntest[j])*untest[j]*ii*_e(d.eff,k*i);
            
  float Utest = 0;
  float Ctest = 0;
  for (int j=0;j<J;j++)
    {
      Utest += .5 * (xntest[j+1] - xntest[j]) * (untest[j+1] + untest[j]);
      Ctest += .5 * (xntest[j+1] - xntest[j]) * (rtest[j+1] + rtest[j]);
    }

  float *otest = (float *)calloc(J+1,sizeof(float));  
  
  for (int j=1;j<=J;j++)
    otest[j] = pow(_c(d.cat,k*i) * (d.p[i][j-1] + (1-d.Qp[i]) * rtest[j]/Ctest) - rtest[j],2.0);
  
  end = clock();  

  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("time: %lf\n",time_spent);
  
  printf("\n%lf\n",Ctest);  
    
  exit(1);
  
  /*for (int j=0;j<J;j++)
    {
      printf("%f ",ret[j]);      
      printf("%f\n",uhtest[j]);      
    }

  end = clock();
  
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("time: %lf\n",time_spent);
  
  /*
  it=0;
  for (int j=0;j<J;j+=8)
    {
      _mm256_store_ps(&ret[j],_xh[it]);
      it+=1;
    }

  printf("\n");
  for (int j=0;j<J;j++)
    printf("%lf\n",ret[j]);  

  it = 0;  
  for (int j=0;j<J;j+=8)
    {
      _mm256_store_ps(&ret[j],_xh_shift[it]);
      it+=1;
    }

  printf("\n");
  for (int j=0;j<J;j++)
    printf("%lf\n",ret[j]);  
  */
  
  //  exit(1);

  /*
  __m256 *_xh_shift;

  posix_memalign(&_xh_shift,32,(J/8)*sizeof(*_xh_shift));

  __m256 *_uh_shift;

  posix_memalign(&_uh_shift,32,(J/8)*sizeof(*_uh_shift));

  __m256 *_ans;

  posix_memalign(&_ans,32,(J/8)*sizeof(*_ans));
  */

  //printf("%lf\n",uh_0);
  //exit(1);


  


  /*
  __m256i pi = _mm256_set_epi32(,7,6,5,4,3,2,1);

  
  it=0;
  for (int j=0;j<J;j+=8)
    {
      _xh_shift[it] = _mm256_permutexvar_ps(pi,_xh[it]);
      _uh_shift[it] = _mm256_permutexvar_ps(pi,_uh[it]);
      it+=1;
    }

  it=0;
  for (int j=0;j<(J-8);j+=8)
    {
      _xh_shift[it] = _mm256_blend_ps(_xh_shift[it],_xh_shift[it+1],0x80);      
      _uh_shift[it] = _mm256_blend_ps(_uh_shift[it],_uh_shift[it+1],0x80);      
      it+=1;
    }
 
  _mm256_store_ps(&xh[0],_xh[0]);
  _mm256_store_ps(&uh[0],_uh[0]);
  
  */


  
    //  for (int j=0;j<J;j++)
  // uh_0 += ret[j];
  
  /*
  it=0;
  int ait=0;
      
  for (int j=0;j<J;j+=8)
    {
      _ans[ait] = _mm256_hadd_ps(_ans[it],_ans[it+1]);
      it += 2;
      ait++;
    }



  it=0;
  ait=0;
      
  for (int j=0;j<J;j+=16)
    {
      _ans[ait] = _mm256_hadd_ps(_ans[it],_ans[it+1]);
      it += 2;
      ait++;
    }




  it=0;
  ait=0;
      
  for (int j=0;j<J;j+=32)
    {
      _ans[ait] = _mm256_hadd_ps(_ans[it],_ans[it+1]);
      it += 2;
      ait++;
    }
  



  it=0;
  ait=0;
      
  for (int j=0;j<J;j+=64)
    {
      _ans[ait] = _mm256_hadd_ps(_ans[it],_ans[it+1]);
      it += 2;
      ait++;
    }

  it=0;
  ait=0;
      
  for (int j=0;j<J;j+=128)
    {
      _ans[ait] = _mm256_hadd_ps(_ans[it],_ans[it+1]);
      it += 2;
      ait++;
    }
  */


  
      /*
  
  for (int jj=3;jj<8;jj++)
    {
    int jplus = pow(2,jj);

      int it=0;
      int ait=0;
      
      for (int j=0;j<J;j+=jplus)
	{
	  _ans[ait] = _mm256_hadd_ps(_ans[it],_ans[it+1]);
	  it += 2;
	  ait++;
	}
  
	}*/



  
  /*
  
      Real uh_0 = xh[0] * _b(a1,a2,xh[0])*uh[0];

      for (int j=0;j<d.J;j++)
	uh_0 += (_b(a1,a2,xh[j])*uh[j] + _b(a1,a2,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);
      
      uh_0 /= (2*kk*ww - xh[0]*_b(a1,a2,0));

  

      _mm256_store_ps(&ret[j],_xh[it]);

  */


  /*
  _xh_shift[0] = _mm256_permutexvar_ps(pi,_xh[0]);

  _xh_shift[1] = _mm256_permutexvar_ps(pi,_xh[1]);  
  
  _xh_shift[0] = _mm256_blend_ps(_xh_shift[0],_xh_shift[1],0x80);
  
  _mm256_store_ps(&ret[0],_xh[0]);
  _mm256_store_ps(&ret[8],_xh[1]);

  for (int j=0;j<10;j++)
    printf("%lf\n",ret[j]);

  _mm256_store_ps(&ret[0],_xh_shift[0]);

  printf("\n");
  for (int j=0;j<8;j++)
    printf("%lf\n",ret[j]);  
  
  exit(1);
  */



  /*

  __m256 tmpx1 = _mm256_i32gather_ps(&(_xh[0]),vindex1,1);
    
  _mm256_store_ps(&(ret[0]),tmpx1);

  for (int j=0;j<8;j++)
    printf("%lf\n",ret[j]);

  _mm256_store_ps(&(ret[8]),_xh[0]);
  _mm256_store_ps(&(ret[16]),_xh[1]);

  printf("\n");
  for (int j=8;j<18;j++)
    printf("%lf\n",ret[j]);

  
  exit(1);

  __m256 mask = _mm256_setzero_ps();
  // __m256 dummy1;      
  //  __m256 tmpx1 = _mm256_mask_i32gather_ps(dummy1,&(_xh[0]),vindex,mask,1);
  */

  
  
}

  

void Kfast(

	   void *args	 
	  )
{

  Parameters *parameters = (Parameters *)args;
  
  Real a1 = parameters->alpha1.value;
  Real a2 = parameters->alpha2.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;
  
  Real * restrict x = (Real *)calloc(d.J+1,sizeof(Real));
  Real * restrict u = (Real *)calloc(d.J+1,sizeof(Real));

  Real * restrict r = (Real *)calloc(d.J+2,sizeof(Real));
  Real * restrict o = (Real *)calloc(d.J+2,sizeof(Real));
  
  /* 
     initialize x
  */  
  for (int j=0;j<=d.J;j++) {
    x[j] = h*j;
  }

  /*
     ok, now initialize u
  */
  
  // prelims
  Real zeta = sqrt( 81*kk*kk*ww*ww*pow(a1+2*a2*ww,2.) - 12*kk*pow(a1*ww+kk,3.) );
  Real eta = 9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta;
  Real Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*kk*(a1*ww+kk) / pow(eta,1./3);

  Real ubar = (Z - bb - kk) / gg; 
  Real vbar = (kk*ww*ubar) / (bb+gg*ubar+kk);
  Real wbar = (2*kk*ww*vbar) / (bb+gg*ubar+2*kk);  
  
  // set
  for (int j=0;j<=d.J;j++)
    {
      u[j] = (a1*vbar+a2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));
      r[j] = 1e-3*w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,0);
    }

  Real ff = 0;
  
  Real U = 0;
  Real C = 0;
  for (int j=0;j<d.J;j++)
    {
      U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
    } 
  
  for (int j=0;j<=d.J;j++) 
    o[j] = pow(_c(d.cat,0) * (d.p[0][j] + (1-d.Qp[0]) * r[j]/C) - r[j],2.0);                      

  for (int j=0;j<d.J;j++)   
    ff += k * .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);
  
  Real * restrict xh = (Real *) calloc(d.J+1,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(d.J+1,sizeof(Real));  

  Real * restrict xn = (Real *) calloc(d.J+2,sizeof(Real));  
  Real * restrict un = (Real *) calloc(d.J+2,sizeof(Real));  

  for (int i=1;i<=d.I;i++)
    {
  
      Real t = k*(i-1);
      Real th = k*(i-.5);
	 
      for (int j=0;j<=d.J;j++)
	{	  
	  xh[j] = x[j] + k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -k/2 * (bb + gg*U + s(x[j])* ii * _e(d.eff,t) - kk) );
	}
      
      // only valid for models where fish are size zero at birth
      Real uh_0 = xh[0] * _b(a1,a2,xh[0])*uh[0];

      for (int j=0;j<d.J;j++)
	uh_0 += (_b(a1,a2,xh[j])*uh[j] + _b(a1,a2,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);

      uh_0 /= (2*kk*ww - xh[0]*_b(a1,a2,0));
            
      Real Uh = .5 * xh[0] * (uh_0 + uh[0]);
      for (int j=0;j<d.J;j++)
	Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);
      
      for (int j=d.J+1;j>0;j--) 
	{
	  xn[j] = x[j-1] + k * kk*(ww-xh[j-1]);      
	  un[j] = u[j-1] * exp ( -k * (bb + gg*Uh + s(xh[j-1]) * ii * _e(d.eff,th) - kk));
	}      

      xn[0] = 0;
      un[0] = x[1] * _b(a1,a2,xn[1])*un[1];
      
      for (int j=1;j<=d.J;j++)
	un[0] += (_b(a1,a2,xn[j])*un[j] + _b(a1,a2,xn[j+1])*un[j+1]) * (xn[j+1]-xn[j]);

      un[0] /= (2*kk*ww - xn[1]*_b(a1,a2,0));
      
      for (int j=0;j<=d.J+1;j++)
	r[j] = 1e-3*w(xn[j])*s(xn[j])*un[j]*ii*_e(d.eff,k*i);
            
      U = 0;
      C = 0;
      for (int j=0;j<=d.J;j++)
	{
	  U += .5 * (xn[j+1] - xn[j]) * (un[j+1] + un[j]);
	  C += .5 * (xn[j+1] - xn[j]) * (r[j+1] + r[j]);
	}
      
      for (int j=0;j<=d.J+1;j++)
	o[j] = pow(_c(d.cat,k*i) * (d.p[i][j] + (1-d.Qp[i]) * r[j]/C) - r[j],2.0);
      
      for (int j=0;j<=d.J;j++)
	ff += k * .5 * (xn[j+1] - xn[j]) * (o[j+1] + o[j]);

      // remove
      for (int j=0;j<idx[i-1];j++)
	x[j] = xn[j];
      for (int j=idx[i-1];j<=d.J;j++)
	x[j] = xn[j+1];      

      for (int j=0;j<idx[i-1];j++)
	u[j] = un[j];
      for (int j=idx[i-1];j<=d.J;j++)
	u[j] = un[j+1];      
      
    }
  
  free(o);
  free(r);
  free(x);
  free(u);
  free(xh);
  free(uh);
  free(xn);
  free(un);

  parameters->ff = ff;
  
}





Real K_no(

	  Parameters *parameters

	  )
{

  Real a1 = parameters->alpha1.value;
  Real a2 = parameters->alpha2.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;
  
  Real * restrict x = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict u = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict r = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict o = (Real *)calloc(d.J,sizeof(Real));
  
  /* 
     initialize x
  */
  int J = (d.J+1) - (d.I+1);
  
  // 'active': first J values in first 'row'
  for (int j=0;j<J;j++) {
    x[j] = h*j;
  }

  // 'neutral': all the rest (J+1.. I+J)
  for (int j=J;j<=d.J;j++)
    x[j] = ww;

  /*
     ok, now initialize u
  */
  
  // prelims
  Real zeta = sqrt( 81*kk*kk*ww*ww*pow(a1+2*a2*ww,2.) - 12*kk*pow(a1*ww+kk,3.) );
  Real eta = 9*a1*kk*kk*ww + 18*a2*kk*kk*ww*ww + kk*zeta;
  Real Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*kk*(a1*ww+kk) / pow(eta,1./3);

  Real ubar = (Z - bb - kk) / gg; 
  Real vbar = (kk*ww*ubar) / (bb+gg*ubar+kk);
  Real wbar = (2*kk*ww*vbar) / (bb+gg*ubar+2*kk);  
  
  // set
  for (int j=0;j<=d.J;j++)
    {
      u[j] = (a1*vbar+a2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));
      r[j] = 1e-3*w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,0);
    }

  Real ff = 0;
  
  Real U = 0;
  Real C = 0;
  for (int j=0;j<d.J;j++)
    {
      U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
    } 
  
  for (int j=0;j<=d.J;j++) 
    o[j] = pow(_c(d.cat,0) * (d.p[0][j] + (1-d.Qp[0]) * r[j]/C) - r[j],2.0);                      
  for (int j=0;j<=d.J;j++)   
    ff += k * .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);
  
  Real * restrict xh = (Real *) calloc(d.J,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(d.J,sizeof(Real));  

  for (int i=1;i<=d.I;i++)
    {
  
      Real t = k*(i-1);
      Real th = k*(i-.5);
	 
      for (int j=0;j<=d.J;j++)
	{	  
	  xh[j] = x[j] + k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -k/2 * (bb + gg*U + s(x[j])* ii * _e(d.eff,t) - kk) );
	}
      
      // only valid for models where fish are size zero at birth
      Real uh_0 = xh[0] * _b(a1,a2,xh[0])*uh[0];

      for (int j=0;j<=d.J;j++)
	uh_0 += (_b(a1,a2,xh[j])*uh[j] + _b(a1,a2,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);

      uh_0 /= (2*kk*ww - xh[0]*_b(a1,a2,0));
            
      Real Uh = .5 * xh[0] * (uh_0 + uh[0]);
      for (int j=0;j<d.J;j++)
	Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);
      
      for (int j=d.J;j>0;j--) 
	{
	  x[j] = x[j-1] + k * kk*(ww-xh[j-1]);      
	  u[j] = u[j-1] * exp ( -k * (bb + gg*Uh + s(xh[j-1]) * ii * _e(d.eff,th) - kk));
	}      

      x[0] = 0;
      u[0] = x[1] * _b(a1,a2,x[1])*u[1];
      
      for (int j=1;j<d.J;j++)
	u[0] += (_b(a1,a2,x[j])*u[j] + _b(a1,a2,x[j+1])*u[j+1]) * (x[j+1]-x[j]);

      u[0] /= (2*kk*ww - x[1]*_b(a1,a2,0));

      for (int j=0;j<=d.J;j++)
	r[j] = 1e-3*w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,k*i);
            
      U = 0;
      C = 0;
      for (int j=0;j<d.J;j++)
	{
	  U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
	  C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
	}
      
      for (int j=0;j<=d.J;j++)
	o[j] = pow(_c(d.cat,k*i) * (d.p[i][j] + (1-d.Qp[i]) * r[j]/C) - r[j],2.0);
      
      for (int j=0;j<=d.J;j++)
	ff += k * .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);
      
    }
  
  free(o);
  free(r);
  free(x);
  free(u);
  free(xh);
  free(uh);

  return ff;
	  
}


/*
Real K_no2(

	  Parameters *parameters

	  )
{

  Real aa = 1; //parameters->alpha.value;
  Real bb = parameters->beta.value;
  Real gg = parameters->gamma.value*1e-7;
  Real kk = parameters->kappa.value;
  Real ww = parameters->omega.value;
  Real ii = parameters->iota.value*1e-3;
  
  Real * restrict x = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict u = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict r = (Real *)calloc(d.J,sizeof(Real));
  Real * restrict o = (Real *)calloc(d.J,sizeof(Real));
  
   
  //   initialize x
  
  int J = (d.J+1) - (d.I+1);
  
  // 'active': first J values in first 'row'
  for (int j=0;j<J;j++) {
    x[j] = h*j;
  }

  // 'neutral': all the rest (J+1.. I+J)
  for (int j=J;j<=d.J;j++)
    x[j] = ww;

  
  //   ok, now initialize u
  
  
  // prelims
  Real zeta = sqrt( 81*kk*kk*ww*ww*pow(aa*A1+2*aa*A2*ww,2.) - 12*kk*pow(aa*A1*ww+kk,3.) );
  Real eta = 9*aa*A1*kk*kk*ww + 18*aa*A2*kk*kk*ww*ww + kk*zeta;
  Real Z = pow(eta,1./3) / (3*pow(2./3,1./3)) + pow(2./3,1./3)*kk*(aa*A1*ww+kk) / pow(eta,1./3);

  Real ubar = (Z - bb - kk) / gg; 
  Real vbar = (kk*ww*ubar) / (bb+gg*ubar+kk);
  Real wbar = (2*kk*ww*vbar) / (bb+gg*ubar+2*kk);  
  
  // set
  for (int j=0;j<=d.J;j++)
    {
      u[j] = (aa*A1*vbar+aa*A2*wbar)*pow(ww-x[j],(bb+gg*ubar)/kk-1) / (kk*pow(ww,(bb+gg*ubar)/kk));
      r[j] = w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,d.k,d.k*(0-d.N));
    }

  Real ff = 0;
  
  Real U = 0;
  Real C = 0;
  for (int j=0;j<d.J;j++)
    {
      U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
      C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
    } 
  
  for (int j=0;j<=d.J;j++) 
    o[j] = pow(_c(d.cat,d.k,d.k*(0-d.N)) * (d.p[0][j] + (1-d.Qp[0]) * r[j]/C) - r[j],2.0);                      
  for (int j=0;j<=d.J;j++)   
    ff += .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);

  //  printf("\n%lf\n",ff);
  
  Real * restrict xh = (Real *) calloc(d.J,sizeof(Real));  
  Real * restrict uh = (Real *) calloc(d.J,sizeof(Real));  

  for (int i=1;i<=d.I;i++)
    {
  
      Real t = d.k*(i-d.N-1);
      Real th = d.k*(i-d.N-.5);
	 
      for (int j=0;j<d.J;j++)
	{	  
	  xh[j] = x[j] + d.k/2 * kk*(ww - x[j]);
	  uh[j] = u[j] * exp( -d.k/2 * (bb + gg*U + s(x[j])* ii * _e(d.eff,d.k,t) - kk) );
	}
      
      // only valid for models where fish are size zero at birth
      Real uh_0 = xh[0] * b(aa,xh[0])*uh[0];

      for (int j=0;j<d.J-1;j++)
	uh_0 += (b(aa,xh[j])*uh[j] + b(aa,xh[j+1])*uh[j+1]) * (xh[j+1]-xh[j]);

      uh_0 /= (2*kk*ww - xh[0]*b(aa,0));
            
      Real Uh = .5 * xh[0] * (uh_0 + uh[0]);
      for (int j=0;j<d.J-1;j++)
	Uh += .5 * (xh[j+1] - xh[j]) * (uh[j+1] + uh[j]);
      
      for (int j=d.J;j>0;j--) 
	{
	  x[j] = x[j-1] + d.k * kk*(ww-xh[j-1]);      
	  u[j] = u[j-1] * exp ( -d.k * (bb + gg*Uh + s(xh[j-1]) * ii * _e(d.eff,d.k,th) - kk));
	}      

      x[0] = 0;
      u[0] = x[1] * b(aa,x[1])*u[1];
      
      for (int j=1;j<d.J;j++)
	u[0] += (b(aa,x[j])*u[j] + b(aa,x[j+1])*u[j+1]) * (x[j+1]-x[j]);

      u[0] /= (2*kk*ww - x[1]*b(aa,0));

      for (int j=0;j<=d.J;j++)
	r[j] = w(x[j])*s(x[j])*u[j]*ii*_e(d.eff,d.k,d.k*(i-d.N));
            
      U = 0;
      C = 0;
      for (int j=0;j<d.J;j++)
	{
	  U += .5 * (x[j+1] - x[j]) * (u[j+1] + u[j]);
	  C += .5 * (x[j+1] - x[j]) * (r[j+1] + r[j]);
	}
      
      for (int j=0;j<=d.J;j++)
	o[j] = pow(_c(d.cat,d.k,d.k*(i-d.N)) * (d.p[i][j] + (1-d.Qp[i]) * r[j]/C) - r[j],2.0);
      
      for (int j=0;j<=d.J;j++)
	ff += .5 * (x[j+1] - x[j]) * (o[j+1] + o[j]);

      //      printf("%lf\n",ff);
      
    }

  // exit(1);
  
  free(o);
  free(r);
  free(x);
  free(u);
  free(xh);
  free(uh);

  return ff;
	  
}
*/

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
  
  solve(parameters,d->eff,d->k,d->S,d->Y,&core_args);

  MeMAT *x = core_args.x;
  MeMAT *u = core_args.u;

  Real iota = parameters->iota.value;

  iota *= 1e-3;
  int S = d->S;
  Real k = d->k;

  int lfi=0;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      MeVEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      MeVEC *v = v_get(x->n);

      for (int j=0;j<x->n;j++)
  v->ve[j] = iota*e(d->eff,k,k*(i-S),d->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < d->n && d->t_id[lfi]==i) 
  {
      
    MeVEC *dt = v_get(d->t_sz[lfi]);

    for (int j=0;j<dt->dim;j++)
      dt->ve[j] = d->lf[lfi][j];

    Real bw = get_bw(dt);

    MeVEC *l = v_get(xt->dim);

    for (int j=0;j<xt->dim;j++)
      for (int jj=0;jj<dt->dim;jj++)
        l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

    Real al = 1e3*c(d->cat,k,k*(i - S)) / Q(xt,l); 

    MeVEC *ld = v_get(x->n);

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


  if (MESCHACH)
    solve(parameters,d->eff,d->k,d->S,d->Y,core_args);
  else
    solve_clean(parameters,d->eff,d->k,d->S,d->Y,core_args);
  
  MeMAT *x = core_args->x;
  MeMAT *u = core_args->u;

  Real iota = parameters->iota.value;

  iota *= 1e-3;
  int S = d->S;
  Real k = d->k;

  int lfi=0;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  int J; 
  int bigJ;
  if (!SGNM){
    bigJ = x->n - 1;
    J = x->n - x->m;
  } else {
    J = x->n - 1;   
  }
  
  for (int i=S;i<x->m;i++)
    {

      int terminator;
      if(!SGNM)
	terminator = J+i;
      else
	terminator = J;

      MeVEC *xt = v_get(terminator+1);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(terminator+1);
      get_row(u,i,ut);
      
      if (!SGNM)
	{
	  xt = v_resize(xt,terminator+1);
	  ut = v_resize(ut,terminator+1);
	}
      
      MeVEC *v = v_get(terminator+1);

      for (int j=0;j<=terminator;j++)
	v->ve[j] = iota*e(d->eff,k,k*(i-S),d->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < d->n && d->t_id[lfi]==i) 
	{
      
 	  MeVEC *dt = v_get(d->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = d->lf[lfi][j];

	  Real bw = get_bw(dt);

	  MeVEC *l = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(d->cat,k,k*(i - S)) / Q(xt,l); 

	  MeVEC *ld = v_get(terminator+1);

	  //printf("\n");
	  for (int j=0;j<=terminator;j++){
	    
	    ld->ve[j] = pow(v->ve[j] - al*l->ve[j],2.);
		    //ld->ve[j] = fabs(v->ve[j] - al*l->ve[j]);
	    //printf("%Lf %Lf\n",xt->ve[j],ld->ve[j]);
	  }
	  //exit(1);
	  
	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	} 
      else 
	{
	  ht->ve[i] = pow(Q(xt,v)-1e3*c(d->cat,k,k*(i - S)),2.);
		  //ht->ve[i] = fabs(Q(xt,v)-1e3*c(d->cat,k,k*(i - S)));
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

	 MeMAT *p,
	 MeMAT *x,
	 MeMAT *u,
	 Data *data,
	 Real iota
	
	 )
{

  iota *=1e-3;
  int lfi=0;

  int S = data->S;
  Real k = data->k;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  int J; 
  int bigJ;
  if (!SGNM){
    bigJ = x->n - 1;
    J = x->n - x->m;

  } else {
    J = x->n - 1;
   
  }

  
  for (int i=S;i<x->m;i++)
    {

      int terminator;
      if(!SGNM)
	terminator = J+i;
      else
	terminator = J;

      MeVEC *xt = v_get(terminator+1);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(terminator+1);
      get_row(u,i,ut);      
      MeVEC *pt = v_get(terminator+1);
      get_row(p,i,pt); 

      if (!SGNM)
	{
	  xt = v_resize(xt,terminator+1);
	  ut = v_resize(ut,terminator+1);
	  pt = v_resize(pt,terminator+1);
	}
      
      MeVEC *v = v_get(terminator+1);
      MeVEC *pv = v_get(terminator+1);

      for (int j=0;j<=terminator;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S),data->Y)*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j] + 1e-3*e(data->eff,k,k*(i-S),data->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S),data->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      //      if(data->t_id[lfi]==i) 
      if(lfi < data->n && data->t_id[lfi]==i) 
	{
      
	  MeVEC *dt = v_get(data->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = data->lf[lfi][j];

	  Real bw = get_bw(dt);

	  MeVEC *l = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(data->cat,k,k*(i - S)) / Q(xt,l);

	  MeVEC *ld = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++) {
	    ld->ve[j] = 2*(v->ve[j] - al*l->ve[j])*pv->ve[j];

	    /*
	    if (v->ve[j] < al*l->ve[j])
	      ld->ve[j] = -pv->ve[j];
	    else
	      ld->ve[j] = pv->ve[j];
	    */
	  }
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

	    MeMAT *p,
	    MeMAT *x,
	    MeMAT *u,
	    Data *data,
	    Real iota

	    )
{

  int S = data->S;
  Real k = data->k;

  int lfi=0;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  for (int i=S;i<x->m;i++)
    {

      MeVEC *xt = v_get(x->n);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(x->n);
      get_row(u,i,ut);      
      MeVEC *pt = v_get(x->n);
      get_row(p,i,pt); 

      MeVEC *v = v_get(x->n);
      MeVEC *pv = v_get(x->n);

      for (int j=0;j<x->n;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S))*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      if(data->t_id[lfi]==i) 
	{
      
	  MeVEC *l = v_get(xt->dim);
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

	    MeMAT *p,
	    MeMAT *x,
	    MeMAT *u,
	    Data *data,
	    Real iota

	    )
{

  iota *=1e-3;
  int S = data->S;
  Real k = data->k;

  int lfi=0;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  int J; 
  int bigJ;
  
  if (!SGNM){
    bigJ = x->n - 1;
    J = x->n - x->m;

  } else {
    J = x->n - 1;
   
  }

  for (int i=S;i<x->m;i++)
    {

      int terminator;
      if(!SGNM)
	terminator = J+i;
      else
	terminator = J;

      MeVEC *xt = v_get(terminator+1);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(terminator+1);
      get_row(u,i,ut);      
      MeVEC *pt = v_get(terminator+1);
      get_row(p,i,pt); 

      if (!SGNM)
	{
	  xt = v_resize(xt,terminator+1);
	  ut = v_resize(ut,terminator+1);
	  pt = v_resize(pt,terminator+1);
	}
      
      MeVEC *v = v_get(terminator+1);
      MeVEC *pv = v_get(terminator+1);
      
      for (int j=0;j<=terminator;j++) 
	{
	  pv->ve[j] = iota*e(data->eff,k,k*(i-S),data->Y)*s(xt->ve[j])*w(xt->ve[j])*pt->ve[j];
	  v->ve[j] = iota*e(data->eff,k,k*(i-S),data->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];
	}

      if(lfi < data->n && data->t_id[lfi]==i) 
	{
      
	  MeVEC *dt = v_get(data->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = data->lf[lfi][j];

	  Real bw = get_bw(dt);

	  MeVEC *l = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(data->cat,k,k*(i - S)) / Q(xt,l);

	  MeVEC *ld = v_get(terminator+1);

	  //printf("\n");
	  for (int j=0;j<=terminator;j++)
	    {
	    
	      ld->ve[j] = 2*(v->ve[j] - al*l->ve[j])*pv->ve[j];
	    //printf("%Lf %Lf\n",xt->ve[j],pv->ve[j]);
	  
	  //exit(1);
	  /*
	      //	      for (int j=0;j<xt->dim;j++)
	      if (v->ve[j] < al*l->ve[j])
		ld->ve[j] = -pv->ve[j];
	      else
	      ld->ve[j] = pv->ve[j];*/
	    }

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
	  ht->ve[i] = Q(xt,pv);*/
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

/*
Real newK(
	  
     Parameters *parameters,
		 Data *d,
		 Solve_Core_Args *core_args
		 
	
		 )
{


  solve(parameters,d->eff,d->k,d->S,d->Y,core_args);

  MeMAT *x = core_args->x;
  MeMAT *u = core_args->u;

  Real iota = parameters->iota.value;

  iota *= 1e-3;
  int S = d->S;
  Real k = d->k;

  int lfi=0;

  MeVEC *ht = v_get(x->m);
  MeVEC *tt = v_get(x->m);

  int J; 
  int bigJ;
  if (BIGMeMATRICES){
    bigJ = x->n - 1;
    J = x->n - x->m;

  } else {
    J = x->n - 1;
   
  }
  
  for (int i=S;i<x->m;i++)
    {

      int terminator;
      if(BIGMeMATRICES)
	terminator = J+i;
      else
	terminator = J;

      MeVEC *xt = v_get(terminator+1);
      get_row(x,i,xt);      
      MeVEC *ut = v_get(terminator+1);
      get_row(u,i,ut);
      
      if (BIGMeMATRICES)
	{
	  xt = v_resize(xt,terminator+1);
	  ut = v_resize(ut,terminator+1);
	}
      
      MeVEC *v = v_get(terminator+1);

      for (int j=0;j<=terminator;j++)
	v->ve[j] = iota*e(d->eff,k,k*(i-S),d->Y)*s(xt->ve[j])*w(xt->ve[j])*ut->ve[j];

      if(lfi < d->n && d->t_id[lfi]==i) 
	{
      
 	  MeVEC *dt = v_get(d->t_sz[lfi]);

	  for (int j=0;j<dt->dim;j++)
	    dt->ve[j] = d->lf[lfi][j];

	  Real bw = get_bw(dt);

	  MeVEC *l = v_get(terminator+1);

	  for (int j=0;j<=terminator;j++)
	    for (int jj=0;jj<dt->dim;jj++)
	      l->ve[j] += exp( -pow((xt->ve[j] - dt->ve[jj])/bw,2.) );

	  Real al = 1e3*c(d->cat,k,k*(i - S)) / Q(xt,l); 

	  MeVEC *ld = v_get(terminator+1);

	  //printf("\n");
	  for (int j=0;j<=terminator;j++){
	    
	    ld->ve[j] = pow(v->ve[j] - al*l->ve[j],2.);
		    //ld->ve[j] = fabs(v->ve[j] - al*l->ve[j]);
	    //printf("%Lf %Lf\n",xt->ve[j],ld->ve[j]);
	  }
	  //exit(1);
	  
	  ht->ve[i] = Q(xt,ld);

	  lfi += 1;

	  V_FREE(dt);
	  V_FREE(l);
	  V_FREE(ld);

	}
*/
