
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include <R.h>
#define MATHLIB_STANDALONE 1
//#define xmalloc(nbytes) malloc_or_exit(nbytes, __FILE__, __LINE__)
#define len 4 * 300 * 4
#define offset2(x, y, xlim, ylim) ((unsigned long long) ((x) * (ylim) + (y)))
#define offset3(x, y, z, xlim, ylim, zlim) ((unsigned long long)((x) * (ylim) * (zlim) + offset2(y, z, ylim, zlim)))
#define offset4(w, x, y, z, wlim, xlim, ylim, zlim) ((unsigned long long)((w) * (xlim) * (ylim) * (zlim) + offset3(x, y, z, xlim, ylim, zlim)))

/*
void *malloc_or_exit(size_t nbytes, const char *file, int line)
{
  void *p=malloc(nbytes);
  if (p==NULL){
    fprintf(stderr, "%s, line %d:" "unable to allocate %lu bytes, calling exit()\n", file, line, (unsigned long long)nbytes);
  }
  return p;
}
*/

double *dvector(int n)
{return (double *)R_alloc(n, sizeof(double));}

unsigned long long perm(unsigned long long n, unsigned long long m){
    if (m == 0) return 1;
    return n * perm(n - 1, m - 1);
}


double _trChat_Cov0(double *Yvec, int N, int P, int T, int s1, int h1, int s2, int h2){
    int i, j, p;
    double trCov = 0;
    for (i = 0; i < N; ++i){
        for (j = 0; j < N; ++j){
            if (i == j) continue;
            double sum1 = 0, sum2 = 0;
            for (p = 0; p < P; ++p){
                sum1 += Yvec[offset3(s1 - 1, i, p, T, N, P)] * Yvec[offset3(s2 - 1, j, p, T, N, P)];
                sum2 += Yvec[offset3(h1 - 1, i, p, T, N, P)] * Yvec[offset3(h2 - 1, j, p, T, N, P)];
            }
            trCov += sum1 * sum2;
        }
    }
    return trCov / perm(N, 2);
}

double _trChat_Cov1(double *Yvec, int N, int P, int T, int s1, int h1, int s2, int h2){
    int i, j, k, p;
    double q1 = 0, q2 = 0, q3 = 0;
    double trCov = 0;
    for (i = 0; i < N; ++i){
        double out_sum1 = 0, out_sum2 = 0;
        for (j = 0; j < N; ++j){
            double sum1 = 0, sum2 = 0;
            for (p = 0; p < P; ++p){
                sum1 += Yvec[offset3(s1 - 1, i, p, T, N, P)] * Yvec[offset3(s2 - 1, j, p, T, N, P)];
                sum2 += Yvec[offset3(h1 - 1, i, p, T, N, P)] * Yvec[offset3(h2 - 1, j, p, T, N, P)];
            }
            q2 += sum1 * sum2;
            out_sum1 += sum1;
            out_sum2 += sum2;
        }
        q1 += out_sum1 * out_sum2;
    }
    for (j = 0; j < N; ++j){
        for (k = 0; k < N; ++k){
            if (k == j) continue;
            double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
            for (p = 0; p < P; ++p){
                sum1 += Yvec[offset3(s1 - 1, j, p, T, N, P)] * Yvec[offset3(s2 - 1, j, p, T, N, P)];
                sum2 += Yvec[offset3(h1 - 1, j, p, T, N, P)] * Yvec[offset3(h2 - 1, k, p, T, N, P)];
                sum3 += Yvec[offset3(s1 - 1, k, p, T, N, P)] * Yvec[offset3(s2 - 1, j, p, T, N, P)];
                sum4 += Yvec[offset3(h1 - 1, k, p, T, N, P)] * Yvec[offset3(h2 - 1, k, p, T, N, P)];
            }
            q3 += sum1 * sum2 + sum3 * sum4;
        }
    }
    trCov = q1 - q2 - q3;
    return trCov / perm(N, 3);
}

double _trChat_Cov2(double *Yvec, int N, int P, int T, int s1, int h1, int s2, int h2){
    int i, j, p;
    double q1 = 0, q2 = 0, q3 = 0, q2_1 = 0, q2_2 = 0, q2_3 = 0, q2_4 = 0, q3_1 = 0, q3_2 = 0;
    double trCov = 0;
    double sum1 = 0, sum2 = 0;
    for (i = 0; i < N; ++i){
        for (j = 0; j < N; ++j){
            if (i == j) continue;
            for (p = 0; p < P; ++p){
                sum1 += Yvec[offset3(s1 - 1, i, p, T, N, P)] * Yvec[offset3(s2 - 1, j, p, T, N, P)];
                sum2 += Yvec[offset3(h1 - 1, i, p, T, N, P)] * Yvec[offset3(h2 - 1, j, p, T, N, P)];
            }
        }
    }
    q1 = sum1 * sum2;
    q2_1 = _trChat_Cov1(Yvec, N, P, T, s1, h1, s2, h2);
    q2_2 = _trChat_Cov1(Yvec, N, P, T, s2, h2, s1, h1);
    q2_3 = _trChat_Cov1(Yvec, N, P, T, s1, h2, s2, h1);
    q2_4 = _trChat_Cov1(Yvec, N, P, T, s2, h1, s1, h2);
    q2 = (q2_1 + q2_2 + q2_3 + q2_4) * perm(N, 3);
    q3_1 = _trChat_Cov0(Yvec, N, P, T, s1, h1, s2, h2);
    q3_2 = _trChat_Cov0(Yvec, N, P, T, s1, h2, s2, h1);
    q3 = (q3_1 + q3_2) * perm(N, 2);
    trCov = q1 - q2 - q3;
    return trCov / perm(N, 4);
}

double trChat(double *Yvec, int N, int P, int T, int s1, int h1, int s2, int h2){
    double trCov0 = 0, trCov1_1 = 0, trCov1_2 = 0, trCov2 = 0;
    // trCov0
    trCov0 = _trChat_Cov0(Yvec, N, P, T, s1, h1, s2, h2);
    // trCov1
    trCov1_1 = _trChat_Cov1(Yvec, N, P, T, s2, h2, s1, h1);
    trCov1_2 = _trChat_Cov1(Yvec, N, P, T, s1, h1, s2, h2);
    // trCov2
    trCov2 = _trChat_Cov2(Yvec, N, P, T, s1, h1, s2, h2);
    // final
    return trCov0 - trCov1_1 - trCov1_2 + trCov2;
}

double _F(double *trCovvec, int T, int s1, int s2, int h1, int h2){
    double tr1111, tr1112, tr1122, tr1211, tr1212, tr1221, tr1222, tr2211, tr2212, tr2222;
    tr1111 = trCovvec[offset4(s1 - 1, h1 - 1, s1 - 1, h1 - 1, T, T, T, T)];
    tr1112 = trCovvec[offset4(s1 - 1, h1 - 1, s1 - 1, h2 - 1, T, T, T, T)];
    tr1122 = trCovvec[offset4(s1 - 1, h2 - 1, s1 - 1, h2 - 1, T, T, T, T)];
    tr1211 = trCovvec[offset4(s1 - 1, h1 - 1, s2 - 1, h1 - 1, T, T, T, T)];
    tr1212 = trCovvec[offset4(s1 - 1, h1 - 1, s2 - 1, h2 - 1, T, T, T, T)];
    tr1221 = trCovvec[offset4(s1 - 1, h2 - 1, s2 - 1, h1 - 1, T, T, T, T)];
    tr1222 = trCovvec[offset4(s1 - 1, h2 - 1, s2 - 1, h2 - 1, T, T, T, T)];
    tr2211 = trCovvec[offset4(s2 - 1, h1 - 1, s2 - 1, h1 - 1, T, T, T, T)];
    tr2212 = trCovvec[offset4(s2 - 1, h1 - 1, s2 - 1, h2 - 1, T, T, T, T)];
    tr2222 = trCovvec[offset4(s2 - 1, h2 - 1, s2 - 1, h2 - 1, T, T, T, T)];
    double ans = (tr1111*tr1111 + tr1122*tr1122 + tr2211*tr2211 + tr2222*tr2222
            - 2*tr1112*tr1112 - 2*tr1211*tr1211 + 2*tr1212*tr1212 + 2*tr1221*tr1221 - 2*tr1222*tr1222 - 2*tr2212*tr2212);
    return ans;
}

double _R(double *trCovvec, int T, int t, int q){
    double ans = 0;
    int s, h;
    for (s = t + 1; s <= T; ++s)
        for (h = q + 1; h <= T; ++h)
            ans += _F(trCovvec, T, t, s, q, h);
    for (s = t + 1; s <= T; ++s)
        for (h = 1; h <= q - 1; ++h)
            ans -= _F(trCovvec, T, t, s, h, q);
    for (s = 1; s <= t - 1; ++s)
        for (h = q + 1; h <= T; ++h)
            ans -= _F(trCovvec, T, s, t, q, h);
    for (s = 1; s <= t - 1; ++s)
        for (h = 1; h <= q - 1; ++h)
            ans += _F(trCovvec, T, s, t, h, q);
    return ans;
}

double _W(double *trCovvec, int T, int t, int q){
    int s, h1, h2;
    double ans = 0;
    for (h1 = 1; h1 <= q; ++h1)
        for (h2 = q + 1; h2 <= T; ++h2){
            for (s = t+1; s <= T; ++s)
                ans += _F(trCovvec, T, t, s, h1, h2);
            for (s = 1; s < t; ++s)
                ans -= _F(trCovvec, T, s, t, h1, h2);
        }
    return ans;
}

double* varTnk12_matrix(double *trCovvec, int N, int T){
    // Calculate W
    double* W = dvector((T - 1) * (T - 1));
    int t, q, s, h;
    for (t = 1; t < T; ++t){
        double tmp = 0;
        for (s = t+1; s <= T; ++s)
            for (h = 2; h <= T; ++h)
                tmp += _F(trCovvec, T, t, s, 1, h);
        for (s = 1; s < t; ++s)
            for (h = 2; h <= T; ++h)
                tmp -= _F(trCovvec, T, s, t, 1, h);
        W[offset2(t - 1, 0, T - 1, T - 1)] = tmp;
    }
    for (q = 2; q < T; ++q)
        for (t = 1; t < T; ++t)
            W[offset2(t - 1, q - 1, T - 1, T - 1)] = W[offset2(t - 1, q - 2, T - 1, T - 1)] + _R(trCovvec, T, t, q);

    // Calculate V
    double* V = dvector((T - 1) * (T - 1));
    /// V(1,1)
    double tmp = 0;
    for (s = 2; s <= T; ++s)
        for (h = 2; h <= T; ++h)
            tmp += _F(trCovvec, T, 1, s, 1, h);
    V[offset2(0, 0, T - 1, T - 1)] = tmp;
    /// V(1, q)
    for (q = 2; q < T; ++q){
        tmp = V[offset2(0, q-2, T-1, T-1)];
        for (s = 2; s <= T; ++s){
            for (h = 1; h < q; ++h)
                tmp -= _F(trCovvec, T, 1, s, h, q);
            for (h = q+1; h <= T; ++h)
                tmp += _F(trCovvec, T, 1, s, q, h);
        }
        V[offset2(0, q-1, T-1, T-1)] = tmp;
    }
    for (t = 2; t < T; ++t){
        for (q = 1; q < t; ++q)
            V[offset2(t - 1, q - 1, T - 1, T - 1)] = V[offset2(q - 1, t - 1, T - 1, T - 1)];
        for (q = t; q < T; ++q)
            V[offset2(t - 1, q - 1, T - 1, T - 1)] = V[offset2(t - 2, q - 1, T - 1, T - 1)] + W[offset2(t - 1, q - 1, T - 1, T - 1)];
    }

    for (t = 1; t < T; ++t)
        for (q = 1; q < T; ++q)
            V[offset2(t - 1, q - 1, T - 1, T - 1)] *= 4.0 / ((long double)N * N * t * q * (T - t) * (T - q));

    return V;
}

int Dkhat(double *Yvec,int n, int p, int T, double *Dk, int *khat)
{
  int k,s,k0,h,h1;
  double *trSigk = dvector(T),Ck,Ck1,Ck2,mk,mTk,largest;

  for (k=0; k<T; ++k)
    trSigk[k]=trChat(Yvec, n, p, T, k + 1, k + 1, k + 1, k + 1);

    // initialize when k = 0;
  Ck=0;
    for (s=0; s<T-1; ++s)
        Ck += trChat(Yvec, n, p, T, 1, 1, s + 2, s + 2);

    mTk=0.0;
    for (h1=0; h1<T-1; ++h1)
        mTk += trSigk[h1+1] / (T-1);

    Dk[0]= trSigk[0] + mTk - 2.0 * Ck / (T - 1);

    // when k is not 0
  for (k=1;k<T-1;k++) {
        Ck1 = Ck2 = 0;
    k0 = k+1;

        mk=0.0;
    for (h=0; h<k; h++){
      mk += trSigk[h];
            Ck1 += trChat(Yvec, n, p, T, k+1, k+1, h+1, h+1);
    }

    mTk=0.0;
    for (h1=0; h1<T-k0; h1++){
      mTk += trSigk[h1+k0];
            Ck2 += trChat(Yvec, n, p, T, h1 + k0 + 1, h1 + k0 + 1, k+1, k+1);
    }

        Dk[k]=(k*(T-k)*Dk[k-1] + (T-2*k0+1)*trSigk[k] - mk + mTk + 2*Ck1 - 2*Ck2)/(k0*(T-k0));
    }

  khat[0] = 1;
  largest = Dk[0];
    for (k=1; k < T - 1; ++k)
        if (Dk[k]>largest){
            khat[0]=k+1;
            largest=Dk[k];
        }

  //R_Free(trSigk);
  return 0;
}

int testandCP(double *Yvec, int n, int p, int T, int *khat, double *stdDk, double *maxstdDk, double *CorrMat, int thread_count)
{
   int T0=T-1, k1, k2, s, s1, h, h1, idx;
   double *tr = dvector((unsigned long long)T*T*T*T), *Dk=dvector(T-1), *Varkvec=dvector((T-1) * (T-1)), *CovMat=dvector(T0*(T0+1)/2), trvecval;

   Dkhat(Yvec,n, p, T, Dk, khat);

   /* Execute the below four for loops in parallel.
      Private variables are declared so there is no conflict as multiple threads
      execute.
   */

   omp_set_num_threads(thread_count);


  #pragma omp parallel
    {
    #pragma omp for schedule(dynamic) private(s1, h, h1, trvecval) collapse(4)
      for (s=1;s<(T+1);s++)
      {
          for (s1=1; s1<(T+1); s1++)
          {
              for (h=1; h<(T+1); h++)
              {
                  for (h1=1; h1<(T+1); h1++)
                  {
                      trvecval=trChat(Yvec, n, p, T, s, s1, h, h1);
                      tr[(h1-1) + T*(h-1) + (unsigned long long)T*T*(s1-1) + (unsigned long long)T*T*T*(s-1)]=trvecval;
                  }
              }
          }
      }
    } // pragma parallel


   Varkvec = varTnk12_matrix(tr, n , T);

    for (idx = 0, s = 1; s < T; ++s)
        for (h = s; h < T; ++h, ++idx)
            CovMat[idx]=Varkvec[offset2(s-1,h-1,T-1,T-1)];

   maxstdDk[0]=Dk[0]/sqrt(CovMat[0]);

    for (k1=0; k1<T0; k1++){
        stdDk[k1]=Dk[k1]/sqrt(CovMat[(2*T0-k1+1)*k1/2]);
        if (stdDk[k1] > maxstdDk[0])
            maxstdDk[0] = stdDk[k1];
    }

    for (k1=0;k1<T0;k1++)
        for (k2=0;k2<T0-k1;k2++)
          CorrMat[(2*T0-k1+1)*k1/2+k2] /= sqrt(CovMat[(2*T0-k1+1)*k1/2]*CovMat[(2*T0-k1-k2+1)*(k1+k2)/2]);

   //R_Free(Dk);
   //R_Free(tr);
   //R_Free(CovMat);
   return 0;
}


void testandCP4r(double *Yvec, int *n, int *p, int *T, int *khat, double *stdDk, double *maxstdDk, double *CorrMat, int *thread_count)
{
  testandCP(Yvec, n[0], p[0], T[0], khat, stdDk, maxstdDk, CorrMat, thread_count[0]);
}

/*
 *
 *  To compile the code run the following lines:
 *
 *   export PKG_CFLAGS="-fopenmp"
 *   export PKG_LIBS="-lgomp"
 *   R CMD SHLIB test-and-cpi-new2.c
 *
*/

