#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

#define QEPS 1e-6
#define MPI 3.1415926
#define MPI1 0.1591549   //1.0/(2*MPI)
#define MEPS 1e-10
#define EPS 1.0e-8
#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))

//-------------- Cox subgroup test ---------------------------------------------
double** vec2matrix(double *vec, int n, int p) {
    // vec is a array with n*p elements, by row saved as a matrix
    double **M;
    M = (double**)malloc(n*sizeof(double*));
    for ( int i = 0; i < n; i++ ) M[i] = vec + i*p;
    return M;
}

int** ivec2matrix(int *vec, int n, int p) {
    // vec is a array with n*p elements, by row saved as a matrix
    int **M;
    M = (int**)malloc(n*sizeof(int*));
    for ( int i = 0; i < n; i++ ) M[i] = vec + i*p;
    return M;
}

// for benchmarks --------------------------------------------------------------
void calculate_v4score(double *exb, int *status, int *yrmin, int *yrmax, double *eta, double *vec, int n) {
    // Y is ascending
    // j >= yrmin[i] equal to Y[j] >= Y[i]
    // j <= yrmax[i] equal to Y[j] <= Y[i]
    int i, j;
    double tmp;

    tmp = 0.0;
    for ( j = i = n-1; i >= 0; i--) {
        for ( ; j >= yrmin[i]; j--) tmp += exb[j];
        eta[i] = tmp;
    }

    tmp = 0.0;
    for ( j = i = 0; i < n; i++) {
        for ( ; j <= yrmax[i]; j++) tmp += status[j] / eta[j];
        vec[i] = status[i] - tmp*exb[i];
    }
}

void calculate_score(
    double **X, double *exb, int *status, int *yrmin, int *yrmax, double *eta, 
    double *score, int n, int p) {
    // Y is ascending
    // j >= yrmin[i] equal to Y[j] >= Y[i]
    // j <= yrmax[i] equal to Y[j] <= Y[i]
    int i, j, k;
    double tmp;
    double vec;

    for ( i = 0; i < p; i++) score[i] = 0.0;

    tmp = 0.0;
    for ( j = i = n-1; i >= 0; i--) {
        for ( ; j >= yrmin[i]; j--) tmp += exb[j];
        eta[i] = tmp;
    }

    tmp = 0.0;
    for ( j = i = 0; i < n; i++) {
        for ( ; j <= yrmax[i]; j++) tmp += status[j] / eta[j];
        vec = status[i] - tmp*exb[i];
        for ( k = 0; k < p; k++) score[k] += vec * X[i][k];
    }
}

void cholesky_decomposition_(double **M, int n)
{
    // M is a n by n symmetry matrix, upper triangle matrix saved is enough
    int i, j, k;
    double tmp;

    for ( i = 0; i < n; i++ ) {
        tmp = M[i][i];
        for ( j = 0; j < i; j++ ) tmp -= M[j][i] * M[j][i];
        tmp = sqrt(fabs(tmp));
        for ( j = 0; j < i; j++ ) {
            for ( k = i; k < n; k++ ) M[i][k] -= M[j][k] * M[j][i];
        }
        if ( tmp > QEPS ) {
            for ( k = i; k < n; k++ ) M[i][k] /= tmp;
        }
    }
}

int MatrixInvSymmetric(double **a,int n){
	int i,j,k,m;
    double w,g,*b;
    b = (double*)malloc(n*sizeof(double));

    for (k=0; k<=n-1; k++){
        w=a[0][0];
        if (fabs(w)+1.0==1.0){
            free(b); return(-1);
        }
        m=n-k-1;
        for (i=1; i<=n-1; i++){
            g=a[i][0]; b[i]=g/w;
            if (i<=m) b[i]=-b[i];
            for (j=1; j<=i; j++)
                a[i-1][j-1]=a[i][j]+g*b[j];
        }
        a[n-1][n-1]=1.0/w;
        for (i=1; i<=n-1; i++)
            a[n-1][i-1]=b[i];
    }
    for (i=0; i<=n-2; i++)
        for (j=i+1; j<=n-1; j++)
            a[i][j]=a[j][i];
    free(b);
    return(0);
}


void calculate_Hessian_inv(
    double **X, double *exb, int *status, double *eta, int *yrmin, int *yrmax, 
    double *ejXj, double **Hessian, int n, int p
    ) {
    // Y is ascending
    // j >= yrmin[i] equal to Y[j] >= Y[i]
    // j <= yrmax[i] equal to Y[j] <= Y[i]
    int i, j, i1, i2;
    double tmp;

    for ( i1 = 0; i1 < p; i1++ ) {
        ejXj[i1] = 0.0;
        Hessian[i1][i1] = 0.0;
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i1][i2] = 0.0;
        }
    }

    double dtmp1 = 0.0;
    for ( j = i = 0; i < n; i++ ) {
        for ( ; j <= yrmax[i]; j++) dtmp1 += status[j] / eta[j];

        tmp = exb[i] * dtmp1;
        for ( i1 = 0; i1 < p; i1++ ) {
            Hessian[i1][i1] += tmp * X[i][i1] * X[i][i1];
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i1][i2] += tmp * X[i][i1] * X[i][i2];
            }
        }
    }

    for ( j = i = n-1; i >= 0; i-- ) {
        for ( ; j >= yrmin[i]; j--) {
            for ( i1 = 0; i1 < p; i1++ ) {
                ejXj[i1] += X[j][i1] * exb[j];
            }
        }

        if (status[i] == 0) continue;

        tmp = 1.0 / ( eta[i] * eta[i] );
        for ( i1 = 0; i1 < p; i1++ ) {
            Hessian[i1][i1] -= tmp * ejXj[i1] * ejXj[i1];
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i1][i2] -= tmp * ejXj[i1] * ejXj[i2];
            }
        }
    }

    for ( i1 = 0; i1 < p; i1++ ) {
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i2][i1] = Hessian[i1][i2];
        }
    }

    // for (i=0; i<p; i++){
    //     for (j = 0; j < p; j++ ) {
    //         Rprintf("%f   ",Hessian[i][j]);
    //     }
    //     Rprintf("\n");
    // }
    // Rprintf("\n");

    MatrixInvSymmetric(Hessian, p);
}


void calculate_D(double *exb, int *status, double *eta, int *yrmax, double **D, int n) {
    // Y is ascending
    // j >= yrmin[i] equal to Y[j] >= Y[i]
    // j <= yrmax[i] equal to Y[j] <= Y[i]
    int i, j, k;

    // for ( i = 0; i < n; i++) {
    //     D[i][i] = 0.0;
    //     for ( j = i+1; j < n; j++ ) {
    //         D[j][i] = D[i][j] = 0.0;
    //     }
    // }
    
    double dtmp1 = 0.0;
    double dtmp2 = 0.0;
    double dtmp3;
    for ( j = i = 0; i < n; i++) {
        for ( ; j <= yrmax[i]; j++) {
            if (status[j]) {
                // dtmp3 = status[j] / eta[j];
                dtmp3 = 1.0 / eta[j];
                dtmp1 += dtmp3;
                dtmp2 += dtmp3*dtmp3;
            }
        }

        // update D
        dtmp3 = -exb[i]*dtmp2;
        D[i][i] = exb[i]*(dtmp1 + dtmp3);
        for ( k = i+1; k < n; k++ ) {
            D[k][i] = D[i][k] = exb[k]*dtmp3;
        }
    }
}

// =============================================================================
// for EstCoxphNewton ----------------------------------------------------------
// =============================================================================
void calculate_exb(double *CoxCoef, double **x, int n, int px, double *exb) {
    int i, j;
    for ( i = 0; i < n; i++ ) {
        exb[i] = 0.0;
        for ( j = 0; j < px; j++ ) {
            exb[i] += x[i][j] * CoxCoef[j];
        }
        exb[i] = exp(exb[i]);
    }
}

void EstCoxphNewton(
    double *beta0, double **x, int *status, int *yrmin, int *yrmax, int n, int p, double tol
){
    int i, j;
    double *score, *exb, *ejXj, *eta, *Hessianvec;

    score = (double*)malloc(sizeof(double)*p);
    exb   = (double*)malloc(sizeof(double)*n);
    ejXj  = (double*)malloc(sizeof(double)*p);
    eta   = (double*)malloc(sizeof(double)*n);
    Hessianvec = (double*)malloc(sizeof(double)*p*p);
    double **Hessian = vec2matrix(Hessianvec, p, p);

    double tmp;
    double increment;

    for (j = 0; j < p; j++) beta0[j] = 0.0;

    increment = tol + 1.0;
    int k = 0;
    int dodc = 0;
    while (increment > tol) {
        calculate_exb(beta0, x, n, p, exb);
        calculate_score(x, exb, status, yrmin, yrmax, eta, score, n, p);
        calculate_Hessian_inv(x, exb, status, eta, yrmin, yrmax, ejXj, Hessian, n, p);
    
        increment = 0.0;
        for ( i = 0; i < p; i++ ) {
            tmp = 0.0;
            for ( j = 0; j < p; j++ ) tmp += Hessian[i][j] * score[j];
            if (fabs(tmp) > increment) increment = fabs(tmp);
            beta0[i] += tmp;
        }

        k++;
        if (k > 20) {
            dodc = 1;
            break;
        }
    }

    if (dodc == 1) {
        Rprintf("Newton method not convergence in 20 steps, try DC method instand.\n");
        increment = tol + 1.0;
        for (j = 0; j < p; j++) beta0[j] = 0.0;
        while (increment > tol) {
            calculate_exb(beta0, x, n, p, exb);
            calculate_score(x, exb, status, yrmin, yrmax, eta, score, n, p);
        
            increment = 0.0;
            for ( i = 0; i < p; i++ ) {
                tmp = score[i] / n * 0.1;
                if (fabs(tmp) > increment) increment = fabs(tmp);
                beta0[i] += tmp;
            }
        }
    }

    free(Hessianvec);
    free(Hessian);
    free(eta);
    free(score);
    free(exb);
    free(ejXj);
}

SEXP _EstCoxphNewton(
    SEXP X_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N, SEXP P, SEXP TOL) {
    int n  = INTEGER(N)[0];
    int p  = INTEGER(P)[0];
    double tol   = REAL(TOL)[0];

    SEXP Beta, list, list_names;
    PROTECT(Beta        = allocVector(REALSXP,  p));
    PROTECT(list_names  = allocVector(STRSXP,   1));
    PROTECT(list        = allocVector(VECSXP,   1));
    double **X = vec2matrix(REAL(X_), n, p);

    EstCoxphNewton(
        REAL(Beta), X, INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), n, p, tol
    );

    SET_STRING_ELT(list_names,  0,  mkChar("Beta"));
    SET_VECTOR_ELT(list,        0,  Beta);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(X);
    UNPROTECT(3);
    return list;
}

// =============================================================================
// for SUP ---------------------------------------------------------------------
// =============================================================================
void calculate_xDx_xD_vzM_vX(
    double *Vec, double **D, double **x, double **z, double **xDx, double **xD, 
    double **vzM, double *vx, int n, int px, int pz
) {
    int i, j, k;

    // calculate xD
    for ( k = 0; k < n; k++ ) {
        for ( j = 0; j < px; j++ ) {
            xD[k][j] = 0.0;
            for ( i = 0; i < n; i++ ) {
                xD[k][j] += D[k][i] * x[i][j];
            }
        }
    }

    // calculate xDx
    for ( i = 0; i < px; i++ ) {
        for ( j = i; j < px; j++ ) {
            xDx[i][j] = 0.0;
            for ( k = 0; k < n; k++ ) {
                xDx[i][j] += x[k][i] * xD[k][j];
            }
        }
    }
    for ( i = 0; i < px; i++ ) {
        for ( j = i; j < px; j++ ) {
            xDx[j][i] = xDx[i][j];
        }
    }

    // calculate vx
    for ( j = 0; j < px; j++ ) {
        vx[j] = 0.0;
        for ( k = 0; k < n; k++ ) {
            vx[j] += Vec[k] * x[k][j];
        }
    }
    
    // calculate vzM
    for ( j = 0; j < pz; j++ ) {
        for ( k = 0; k < n; k++ ) {
            vzM[k][j] = Vec[k] * z[k][j];
        }
    }
}

void calculate_SUP_given_id(
    double **xDx, double **xD, double **D, double **z, double **vzM, 
    int *id, int n, int p, int px, int pz, double *Diz, 
    double *Score, double *Hessian, double **HessianM, double *test
) {
    int i, j, k, h;

    for ( i = px; i < p; i++ ) Score[i] = 0.0;
    for ( i = 0; i < p*p; i++ ) Hessian[i] = 0.0;

    for ( k = 0; k < n; k++ ) {
        if (id[k]) {
            for ( j = 0; j < pz; j++ ) Score[j+px] += vzM[k][j];
            for ( i = 0; i < pz; i++ ) {
                for ( j = 0; j < px; j++ ) {
                    HessianM[j][i+px] += xD[k][j] * z[k][i];
                }
            }

            for ( j = 0; j < pz; j++ ) Diz[j] = 0.0;
            for ( h = 0; h < n; h++ ) {
                if (id[h]) {
                    for ( j = 0; j < pz; j++ ) Diz[j] += D[k][h] * z[h][j];
                }
            }
            for ( i = 0; i < pz; i++ ) {
                for ( j = 0; j <= i; j++ ) {
                    HessianM[j+px][i+px] += Diz[j] * z[k][i];
                }
            }
        }
    }

    for ( i = 0; i < px; i++ ) {
        HessianM[i][i] = xDx[i][i];
        for ( j = i+1; j < px; j++ ) {
            HessianM[i][j] = xDx[i][j];
        }
    }

    for ( i = 0; i < p; i++ ) {
        for ( j = i+1; j < p; j++ ) {
            HessianM[j][i] = HessianM[i][j];
        }
    }

    // printArrayDouble(Score, p, p);
    // Rprintf("Hessian.\n");
    // printArrayDouble(Hessian, p*p, p);

    MatrixInvSymmetric(HessianM, p);

    // calculate_test
    test[0] = 0.0;
    for ( i = 0; i < p; i++ ) {
        test[0] += HessianM[i][i]*Score[i]*Score[i];
        for ( j = i+1; j < p; j++ ) {
            test[0] += 2.0*HessianM[i][j]*Score[i]*Score[j];
        }
    }
}

void calculate_SUP00(
    double **xDx, double **xD, double **D, double **z, double **vzM, double *vx, 
    int **idM, int **sample, int n, int p, int px, int pz, int K, int B, 
    double *test, double *testboot, double *pvalue
) {
    int i, j, bb;
    double testi;

    double *Score, *Hessian, *Diz;
    Score = (double*)malloc(sizeof(double)*p);
    Hessian = (double*)malloc(sizeof(double)*p*p);
    Diz = (double*)malloc(sizeof(double)*pz);
    double **HessianM = vec2matrix(Hessian, p, p);
    for ( i = 0; i < px; i++ ) Score[i] = vx[i];

    // SUP test staitsitic
    test[0] = 0.0;
    for ( i = 0; i < K; i++ ) {
        calculate_SUP_given_id(xDx, xD, D, z, vzM, idM[i], n, p, px, pz, Diz, Score, Hessian, HessianM, &testi);
        if (test[0] < testi) test[0] = testi;
    }

    // bootstrap
    int *idBoot;
    idBoot = (int*)malloc(sizeof(int)*n);
    pvalue[0] = 0.0;
    for ( bb = 0; bb < B; bb++ ) {
        testboot[bb] = 0.0;
        for ( i = 0; i < K; i++ ) {
            for ( j = 0; j < n; j++ ) idBoot[j] = idM[i][sample[j][bb]];
            calculate_SUP_given_id(xDx, xD, D, z, vzM, idBoot, n, p, px, pz, Diz, Score, Hessian, HessianM, &testi);
            if (testboot[bb] < testi) testboot[bb] = testi;
        }
        // pvalue
        if (testboot[bb] > test[0]) pvalue[0] += 1.0;
    }
    pvalue[0] /= B;
    
    free(idBoot);
    free(Score);
    free(Diz);
    free(Hessian);
    free(HessianM);
}

SEXP _calculate_SUP00(
    SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP X_, SEXP Z_,  
    SEXP ID_, SEXP SAMPLE_, SEXP N, SEXP PZ, SEXP PX, SEXP B, SEXP K, SEXP TOL_) {
    int n  = INTEGER(N)[0];
    int pz = INTEGER(PZ)[0];
    int px = INTEGER(PX)[0];
    int b  = INTEGER(B)[0];
    int k  = INTEGER(K)[0];
    int p = pz + px;
    double tol = REAL(TOL_)[0];

    SEXP TestR, TestB, Pval, list, list_names;
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(TestB       = allocVector(REALSXP,  b));
    PROTECT(Pval        = allocVector(REALSXP,  1));
    PROTECT(list_names  = allocVector(STRSXP,   3));
    PROTECT(list        = allocVector(VECSXP,   3));

    double *eta, *Dvec, *Vec, *exb, *Beta;
    double **D;
    eta  = (double*)malloc(sizeof(double)*n);
    Vec = (double*)malloc(sizeof(double)*n);
    Dvec  = (double*)malloc(sizeof(double)*n*n);
    D = vec2matrix(Dvec, n, n);
    exb  = (double*)malloc(sizeof(double)*n);
    Beta  = (double*)malloc(sizeof(double)*px);
    double **x = vec2matrix(REAL(X_), n, px);
    EstCoxphNewton(Beta, x, INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), n, px, tol);
    calculate_exb(Beta, x, n, px, exb);
    calculate_v4score(exb, INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), eta, Vec, n);
    calculate_D(exb, INTEGER(Status_), eta, INTEGER(YRMAX_), D, n);
    free(eta);
    free(exb);

    double **z = vec2matrix(REAL(Z_), n, pz);
    double *xDxvec, *xDvec, *vzMvec, *vx;
    double **xDx, **xD, **vzM;
    vx      = (double*)malloc(sizeof(double)*px);
    xDxvec  = (double*)malloc(sizeof(double)*px*px);
    xDvec   = (double*)malloc(sizeof(double)*n*px);
    vzMvec  = (double*)malloc(sizeof(double)*n*pz);
    xDx     = vec2matrix(Dvec,   px, px);
    xD      = vec2matrix(xDvec,  n, px);
    vzM     = vec2matrix(vzMvec, n, pz);
    calculate_xDx_xD_vzM_vX(Vec, D, x, z, xDx, xD, vzM, vx, n, px, pz);
    free(Vec);
    free(x);

    int **sample = ivec2matrix(INTEGER(SAMPLE_), n, b);
    int **idM = ivec2matrix(INTEGER(ID_), k, n);
    calculate_SUP00(xDx, xD, D, z, vzM, vx, idM, sample, n, p, px, pz, k, b, REAL(TestR), REAL(TestB), REAL(Pval));

    SET_STRING_ELT(list_names,  0,  mkChar("TestR"));
    SET_STRING_ELT(list_names,  1,  mkChar("TestB"));
    SET_STRING_ELT(list_names,  2,  mkChar("Pvalue"));
    SET_VECTOR_ELT(list,        0,  TestR);
    SET_VECTOR_ELT(list,        1,  TestB);
    SET_VECTOR_ELT(list,        2,  Pval);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(Beta);
    free(Dvec);
    free(D);
    free(z);
    free(xDxvec);
    free(xDx);
    free(xDvec);
    free(xD);
    free(vzMvec);
    free(vzM);
    free(vx);
    free(sample);
    free(idM);
    UNPROTECT(5);
    return list;
}


void calculate_exb_sub(double *CoxCoef, double **x, double **z, int *id, int n, int px, int pz, double *exb) {
    int i, j;
    double *CoxCoefZ = CoxCoef + px;
    for ( i = 0; i < n; i++ ) {
        exb[i] = 0.0;
        for ( j = 0; j < px; j++ ) {
            exb[i] += x[i][j] * CoxCoef[j];
        }
        if (id[i]) {
            for ( j = 0; j < pz; j++ ) {
                exb[i] += z[i][j] * CoxCoefZ[j];
            }
        }
        exb[i] = exp(exb[i]);
    }
}

void calculate_score_sub(
    double **X, double **Z, int *id, double *exb, int *status, int *yrmin, int *yrmax, double *eta, 
    double *score, int n, int px, int pz) {
    // Y is ascending
    // j >= yrmin[i] equal to Y[j] >= Y[i]
    // j <= yrmax[i] equal to Y[j] <= Y[i]
    int i, j, k;
    double tmp;
    double vec;
    int p = px + pz;

    double *scoreZ = score + px;
    for ( i = 0; i < p; i++) score[i] = 0.0;

    tmp = 0.0;
    for ( j = i = n-1; i >= 0; i--) {
        for ( ; j >= yrmin[i]; j--) tmp += exb[j];
        eta[i] = tmp;
    }

    tmp = 0.0;
    for ( j = i = 0; i < n; i++) {
        for ( ; j <= yrmax[i]; j++) tmp += status[j] / eta[j];
        vec = status[i] - tmp*exb[i];
        for ( k = 0; k < px; k++) score[k] += vec * X[i][k];
        if (id[i]) for ( k = 0; k < pz; k++) scoreZ[k] += vec * Z[i][k];
    }
}

void calculate_Hessian_inv_sub(
    double **x, double **z, int *id, double *exb, int *status, double *eta, int *yrmin, int *yrmax, 
    double *ejXj, double **Hessian, int n, int px, int pz
    ) {
    int i, j, i1, i2;
    double tmp;
    int p = px + pz;

    for ( i1 = 0; i1 < p; i1++ ) {
        ejXj[i1] = 0.0;
        Hessian[i1][i1] = 0.0;
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i1][i2] = 0.0;
        }
    }

    double dtmp1 = 0.0;
    for ( j = i = 0; i < n; i++ ) {
        for ( ; j <= yrmax[i]; j++) dtmp1 += status[j] / eta[j];

        tmp = exb[i] * dtmp1;
        for ( i1 = 0; i1 < px; i1++ ) {
            Hessian[i1][i1] += tmp * x[i][i1] * x[i][i1];
            for ( i2 = i1+1; i2 < px; i2++ ) {
                Hessian[i1][i2] += tmp * x[i][i1] * x[i][i2];
            }
        }
        if (id[i]) {
            for ( i1 = 0; i1 < pz; i1++ ) {
                Hessian[i1+px][i1+px] += tmp * z[i][i1] * z[i][i1];
                for ( i2 = i1+1; i2 < pz; i2++ ) {
                    Hessian[i1+px][i2+px] += tmp * z[i][i1] * z[i][i2];
                }
                for ( i2 = 0; i2 < px; i2++ ) {
                    Hessian[i2][i1+px] += tmp * x[i][i2] * z[i][i1];
                }
            }
        }
    }
    
    for ( j = i = n-1; i >= 0; i-- ) {
        for ( ; j >= yrmin[i]; j--) {
            for ( i1 = 0; i1 < px; i1++ ) {
                ejXj[i1] += x[j][i1] * exb[j];
            }
            if (id[j]) {
                for ( i1 = 0; i1 < pz; i1++ ) {
                    ejXj[i1+px] += z[j][i1] * exb[j];
                }
            }
        }

        if (status[i] == 0) continue;

        tmp = 1.0 / ( eta[i] * eta[i] );
        for ( i1 = 0; i1 < p; i1++ ) {
            Hessian[i1][i1] -= tmp * ejXj[i1] * ejXj[i1];
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i1][i2] -= tmp * ejXj[i1] * ejXj[i2];
            }
        }
    }

    for ( i1 = 0; i1 < p; i1++ ) {
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i2][i1] = Hessian[i1][i2];
        }
    }
    MatrixInvSymmetric(Hessian, p);
}

// =============================================================================
// for LRT ---------------------------------------------------------------------
// =============================================================================
void calculate_LRT_given_id(
    double *beta0, double *score, double *exb, double *ejXj, double *eta, 
    double **x, double **z, int *id, int *status, int *yrmin, 
    int *yrmax, int n, int px, int pz, int b, double tol, double **HessianH0, double **HessianH1minusH0, 
    double **rbyDsqrt, double **rbyDsqrtX, double **vzM,  
    double *ScoreH0, double *test, double *testboot, int maxs, int loops
) {
    // EstCoxphNewtonSubgroup
    int i, j, bb;
    int p = px + pz;

    double tmp;
    double increment;

    for (j = 0; j < p; j++) beta0[j] = 0.0;

    increment = tol + 1.0;
    int k = 0;
    int dodc = 0;
    while (increment > tol) {
        calculate_exb_sub(beta0, x, z, id, n, px, pz, exb);
        calculate_score_sub(x, z, id, exb, status, yrmin, yrmax, eta, score, n, px, pz);
        calculate_Hessian_inv_sub(x, z, id, exb, status, eta, yrmin, yrmax, ejXj, HessianH1minusH0, n, px, pz);

        increment = 0.0;
        for ( i = 0; i < p; i++ ) {
            tmp = 0.0;
            for ( j = 0; j < p; j++ ) tmp += HessianH1minusH0[i][j] * score[j];
            if (fabs(tmp) > increment) increment = fabs(tmp);
            beta0[i] += tmp;
        }

        k++;
        if (k > maxs) {
            dodc = 1;
            break;
        }
    }

    if (dodc == 1) {
        Rprintf("Newton method not convergence in %d steps for %d-th gamma, try DC method instand.\n", maxs, loops+1);
        increment = tol + 1.0;
        for (j = 0; j < p; j++) beta0[j] = 0.0;
        while (increment > tol) {
            calculate_exb_sub(beta0, x, z, id, n, px, pz, exb);
            calculate_score_sub(x, z, id, exb, status, yrmin, yrmax, eta, score, n, px, pz);
        
            increment = 0.0;
            for ( i = 0; i < p; i++ ) {
                tmp = score[i] / n * 0.1;
                if (fabs(tmp) > increment) increment = fabs(tmp);
                beta0[i] += tmp;
            }
        }
    }
    calculate_exb_sub(beta0, x, z, id, n, px, pz, exb);
    calculate_score_sub(x, z, id, exb, status, yrmin, yrmax, eta, score, n, px, pz);
    calculate_Hessian_inv_sub(x, z, id, exb, status, eta, yrmin, yrmax, ejXj, HessianH1minusH0, n, px, pz);
    // HessianH1minusH0
    for ( i = 0; i < px; i++ ) {
        HessianH1minusH0[i][i] -= HessianH0[i][i];
        for (j = i+1; j < px; j++) {
            HessianH1minusH0[i][j] -= HessianH0[i][j];
            HessianH1minusH0[j][i] = HessianH1minusH0[i][j];
        }
    }

    // calculate ScoreH0
    for ( j = px; j < p; j++ ) ScoreH0[j] = 0.0;
    for ( i = 0; i < n; i++ ) {
        if (id[i]) {
            for ( j = 0; j < pz; j++ ) ScoreH0[j+px] += vzM[i][j];
        }
    }

    // calculate_test
    test[0] = 0.0;
    for ( i = 0; i < p; i++ ) {
        test[0] += HessianH1minusH0[i][i]*ScoreH0[i]*ScoreH0[i];
        for ( j = i+1; j < p; j++ ) {
            test[0] += 2.0*HessianH1minusH0[i][j]*ScoreH0[i]*ScoreH0[j];
        }
    }

    // bootstrap
    double *scoreZ = score + px;
    for ( bb = 0; bb < b; bb++) {
        for ( i = 0; i < px; i++ ) score[i] = rbyDsqrtX[bb][i];
        for ( i = px; i < p; i++ ) score[i] = 0.0;
        for ( j = 0; j < n; j++ ) {
            if (id[j]) {
                for ( i = 0; i < pz; i++ ) {
                    scoreZ[i] += rbyDsqrt[bb][j] * z[j][i];
                }
            }
        }

        // calculate_test
        testboot[bb] = 0.0;
        for ( i = 0; i < p; i++ ) {
            testboot[bb] += HessianH1minusH0[i][i]*score[i]*score[i];
            for ( j = i+1; j < p; j++ ) {
                testboot[bb] += 2.0*HessianH1minusH0[i][j]*score[i]*score[j];
            }
        }
    }

    // printArrayDouble(testboot, b, b);
    // printArrayDouble(test, 1, 1);
    // printArrayDouble2(HessianH1minusH0, p, p);
}

void calculate_LRT00(
    double **x, double **z, int **idM, int *status, int *yrmin, int *yrmax, int n, 
    int px, int pz, int b, double tol, double **HessianH0, int K, double *vx, 
    double **rbyDsqrt, double **rbyDsqrtX, double **vzM, double *test, 
    double *testboot, double *pval, int maxs, 
    double *testvec, double *bootvec, int saveall
) {
    int i, j;
    int p = px + pz;

    double *beta0, *score, *exb, *ejXj, *eta, *Hessianvec;
    beta0 = (double*)malloc(sizeof(double)*p);
    score = (double*)malloc(sizeof(double)*p);
    exb   = (double*)malloc(sizeof(double)*n);
    ejXj  = (double*)malloc(sizeof(double)*p);
    eta   = (double*)malloc(sizeof(double)*n);
    Hessianvec = (double*)malloc(sizeof(double)*p*p);
    double **HessianH1minusH0 = vec2matrix(Hessianvec, p, p);

    double *ScoreH0;
    ScoreH0 = (double*)malloc(sizeof(double)*p);
    for ( i = 0; i < px; i++ ) ScoreH0[i] = vx[i];

    double testi;
    double *testbooti;
    testbooti = (double*)malloc(sizeof(double)*b);
    
    test[0] = 0.0;
    for ( j = 0; j < b; j++ ) testboot[j] = 0.0;
    double *bootvecK = bootvec;
    
    for ( i = 0; i < K; i++ ) {
        calculate_LRT_given_id(
            beta0, score, exb, ejXj, eta, x, z, idM[i], status, yrmin, yrmax, 
            n, px, pz, b, tol, HessianH0, HessianH1minusH0, rbyDsqrt, rbyDsqrtX, 
            vzM, ScoreH0, &testi, testbooti, maxs, i
        );
        if (saveall == 1) {
            testvec[i] = testi;
            for ( j = 0; j < b; j++ ) bootvecK[j] = testbooti[j];
            bootvecK += b;
        }
        if (test[0] < testi) test[0] = testi;
        for ( j = 0; j < b; j++ ) {
            if (testboot[j] < testbooti[j]) testboot[j] = testbooti[j];
        }
    }

    // pval
    pval[0] = 0.0;
    for ( j = 0; j < b; j++ ) {
        if (testboot[j] > test[0]) pval[0] += 1.0;
    }
    pval[0] /= b;
    
    free(Hessianvec);
    free(HessianH1minusH0);
    free(beta0);
    free(eta);
    free(score);
    free(exb);
    free(ejXj);
}

void calculate_LRT_Related(
    double **x, double **z, double *exb, int *status, double *eta, int *yrmin, int *yrmax, 
    int n, int px, int pz, int b, double *vec, double **rvmatrix,
    double **HessianH0, double *vx, double **rbyDsqrt, double **rbyDsqrtX, double **vzM
    ) {
    int i, j, k;

    // Hessian_inv
    double *ejXj;
    ejXj  = (double*)malloc(sizeof(double)*px);
    calculate_Hessian_inv(x, exb, status, eta, yrmin, yrmax, ejXj, HessianH0, n, px);
    free(ejXj);

    // calculate vx
    for ( j = 0; j < px; j++ ) {
        vx[j] = 0.0;
        for ( k = 0; k < n; k++ ) {
            vx[j] += vec[k] * x[k][j];
        }
    }

    // calculate vzM
    for ( j = 0; j < pz; j++ ) {
        for ( k = 0; k < n; k++ ) {
            vzM[k][j] = vec[k] * z[k][j];
        }
    }

    double *Dvec;
    Dvec = (double*)malloc(sizeof(double)*n*n);
    for ( i = 0; i < n*n; i++) Dvec[i] = 0.0;
    double **Dsqrt = vec2matrix(Dvec, n, n);
    
    // calculate_Dsqrt
    double dtmp1 = 0.0;
    double dtmp2 = 0.0;
    double dtmp3;
    for ( j = i = 0; i < n; i++) {
        for ( ; j <= yrmax[i]; j++) {
            if (status[j]) {
                // dtmp3 = status[j] / eta[j];
                dtmp3 = 1.0 / eta[j];
                dtmp1 += dtmp3;
                dtmp2 += dtmp3*dtmp3;
            }
        }

        // update D
        dtmp3 = -exb[i]*dtmp2;
        Dsqrt[i][i] = exb[i]*(dtmp1 + dtmp3);
        for ( k = i+1; k < n; k++ ) {
            Dsqrt[i][k] = exb[k]*dtmp3;
        }
    }
    cholesky_decomposition_(Dsqrt, n);

    // rbyDsqrt = rvmatrix %*% Dsqrt, B by n matrix
    for ( i = 0; i < b; i++) {
        for ( j = 0; j < n; j++ ) {
            rbyDsqrt[i][j] = 0.0;
            for ( k = 0; k < n; k++ ) {
                rbyDsqrt[i][j] += rvmatrix[i][k] * Dsqrt[k][j];
            }
        }
    }

    // rbyDsqrtX = rbyDsqrt %*% x, B by px matrix
    for ( i = 0; i < b; i++) {
        for ( j = 0; j < px; j++ ) {
            rbyDsqrtX[i][j] = 0.0;
            for ( k = 0; k < n; k++ ) {
                rbyDsqrtX[i][j] += rbyDsqrt[i][k] * x[k][j];
            }
        }
    }
    
    free(Dvec);
    free(Dsqrt);
}

SEXP _calculate_LRT00(
    SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP X_, SEXP Z_,  SEXP ID_, 
    SEXP Rvmatrix_, SEXP N, SEXP PZ, SEXP PX, SEXP B, SEXP K, SEXP TOL, SEXP MAXSTEP, SEXP SAVEALL
    ) {
    int n  = INTEGER(N)[0];
    int pz = INTEGER(PZ)[0];
    int px = INTEGER(PX)[0];
    int b  = INTEGER(B)[0];
    int k  = INTEGER(K)[0];
    int maxs  = INTEGER(MAXSTEP)[0];
    int saveall  = INTEGER(SAVEALL)[0];
    // int p = pz + px;
    double tol = REAL(TOL)[0];
    double **x = vec2matrix(REAL(X_), n, px);
    double **z = vec2matrix(REAL(Z_), n, pz);
    int *yrmin = INTEGER(YRMIN_);
    int *yrmax = INTEGER(YRMAX_);

    SEXP TestR, TestB, Pval, list, list_names, TestRVEC, TestBVEC;
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(TestB       = allocVector(REALSXP,  b));
    PROTECT(Pval        = allocVector(REALSXP,  1));
    PROTECT(list_names  = allocVector(STRSXP,   5));
    PROTECT(list        = allocVector(VECSXP,   5));

    if (saveall) {
        PROTECT(TestRVEC = allocVector(REALSXP,  k));
        PROTECT(TestBVEC = allocVector(REALSXP,  k*b));
    } else {
        PROTECT(TestRVEC = allocVector(REALSXP,  1));
        PROTECT(TestBVEC = allocVector(REALSXP,  1));
    }

    // calculate H0
    double *BetaH0, *exb, *eta, *Vec;
    BetaH0 = (double*)malloc(sizeof(double)*px);
    exb    = (double*)malloc(sizeof(double)*n);
    eta    = (double*)malloc(sizeof(double)*n);
    Vec    = (double*)malloc(sizeof(double)*n);
    EstCoxphNewton(BetaH0, x, INTEGER(Status_), yrmin, yrmax, n, px, tol);
    calculate_exb(BetaH0, x, n, px, exb);
    free(BetaH0);
    calculate_v4score(exb, INTEGER(Status_), yrmin, yrmax, eta, Vec, n);
    double *HessianH0Vec, *vx, *rbyDsqrtVec, *rbyDsqrtXVec, *vzMvec;
    HessianH0Vec = (double*)malloc(sizeof(double)*px*px);
    vx           = (double*)malloc(sizeof(double)*px);
    rbyDsqrtVec  = (double*)malloc(sizeof(double)*b*n);
    rbyDsqrtXVec = (double*)malloc(sizeof(double)*b*px);
    vzMvec       = (double*)malloc(sizeof(double)*n*pz);
    double **HessianH0 = vec2matrix(HessianH0Vec, px, px);
    double **rbyDsqrt  = vec2matrix(rbyDsqrtVec, b, n);
    double **rbyDsqrtX = vec2matrix(rbyDsqrtXVec, b, px);
    double **vzM       = vec2matrix(vzMvec, n, pz);

    double **rvmatrix = vec2matrix(REAL(Rvmatrix_), b, n);
    calculate_LRT_Related(
        x, z, exb, INTEGER(Status_), eta, yrmin, yrmax, n, px, pz, b, Vec, 
        rvmatrix, HessianH0, vx, rbyDsqrt, rbyDsqrtX, vzM
    );
    free(exb);
    free(eta);
    free(Vec);

    int **idM = ivec2matrix(INTEGER(ID_), k, n);
    calculate_LRT00(
        x, z, idM, INTEGER(Status_), yrmin, yrmax, n, px, pz, b, tol, HessianH0, 
        k, vx, rbyDsqrt, rbyDsqrtX, vzM, REAL(TestR), REAL(TestB), REAL(Pval), maxs,
        REAL(TestRVEC), REAL(TestBVEC), saveall
    );

    SET_STRING_ELT(list_names,  0,  mkChar("TestR"));
    SET_STRING_ELT(list_names,  1,  mkChar("TestB"));
    SET_STRING_ELT(list_names,  2,  mkChar("Pvalue"));
    SET_STRING_ELT(list_names,  3,  mkChar("TestRVEC"));
    SET_STRING_ELT(list_names,  4,  mkChar("TestBVEC"));
    SET_VECTOR_ELT(list,        0,  TestR);
    SET_VECTOR_ELT(list,        1,  TestB);
    SET_VECTOR_ELT(list,        2,  Pval);
    SET_VECTOR_ELT(list,        3,  TestRVEC);
    SET_VECTOR_ELT(list,        4,  TestBVEC);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(x);
    free(z);
    free(HessianH0Vec);
    free(HessianH0);
    free(vx);
    free(rbyDsqrtVec);
    free(rbyDsqrt);
    free(rbyDsqrtXVec);
    free(rbyDsqrtX);
    free(vzMvec);
    free(vzM);
    free(rvmatrix);
    free(idM);
    UNPROTECT(7);
    return list;
}

// =============================================================================
// for ST ---------------------------------------------------------------------
// =============================================================================
void calculate_ST_given_id(
    double **z, double **vzM, double **rvmatrix, 
    int *id, int n, int pz, int b, 
    double *Score, double *Hessian, double **HessianM, double *test, double *testboot
) {
    int i, j, k, bb;

    for ( i = 0; i < pz; i++ ) Score[i] = 0.0;
    for ( i = 0; i < pz*pz; i++ ) Hessian[i] = 0.0;

    for ( k = 0; k < n; k++ ) {
        if (id[k]) {
            for ( j = 0; j < pz; j++ ) Score[j] += vzM[k][j];
            for ( i = 0; i < pz; i++ ) {
                for ( j = i; j < pz; j++ ) {
                    HessianM[i][j] += vzM[k][i] * vzM[k][j];
                }
            }
        }
    }

    for ( i = 0; i < pz; i++ ) {
        for ( j = i+1; j < pz; j++ ) {
            HessianM[j][i] = HessianM[i][j];
        }
    }

    MatrixInvSymmetric(HessianM, pz);

    // calculate_test
    test[0] = 0.0;
    for ( i = 0; i < pz; i++ ) {
        test[0] += HessianM[i][i]*Score[i]*Score[i];
        for ( j = i+1; j < pz; j++ ) {
            test[0] += 2.0*HessianM[i][j]*Score[i]*Score[j];
        }
    }

    // bootstrap
    for ( bb = 0; bb < b; bb++ ) {
        for ( i = 0; i < pz; i++ ) Score[i] = 0.0;
        for ( i = 0; i < n; i++ ) {
            if (id[i]) {
                for ( j = 0; j < pz; j++ ) Score[j] += rvmatrix[bb][i] * vzM[i][j];
            }
        }
        
        // calculate_test
        testboot[bb] = 0.0;
        for ( i = 0; i < pz; i++ ) {
            testboot[bb] += HessianM[i][i]*Score[i]*Score[i];
            for ( j = i+1; j < pz; j++ ) {
                testboot[bb] += 2.0*HessianM[i][j]*Score[i]*Score[j];
            }
        }
    }
}

void calculate_ST00(
    double **z, int **idM, int n, int pz, int b, int K, double **rvmatrix, 
    double **vzM, double *test, double *testboot, double *pval, 
    double *testvec, double *bootvec, int saveall
) {
    int i, j;

    double *score, *Hessianvec;
    score = (double*)malloc(sizeof(double)*pz);
    Hessianvec = (double*)malloc(sizeof(double)*pz*pz);
    double **Hessian = vec2matrix(Hessianvec, pz, pz);

    double testi;
    double *testbooti;
    testbooti = (double*)malloc(sizeof(double)*b);
    
    test[0] = 0.0;
    for ( j = 0; j < b; j++ ) testboot[j] = 0.0;
    double *bootvecK = bootvec;
    
    for ( i = 0; i < K; i++ ) {
        calculate_ST_given_id(
            z, vzM, rvmatrix, idM[i], n, pz, b, score, Hessianvec, Hessian, &testi, testbooti
        );
        if (saveall == 1) {
            testvec[i] = testi;
            for ( j = 0; j < b; j++ ) bootvecK[j] = testbooti[j];
            bootvecK += b;
        }
        if (test[0] < testi) test[0] = testi;
        for ( j = 0; j < b; j++ ) {
            if (testboot[j] < testbooti[j]) testboot[j] = testbooti[j];
        }
    }

    // pval
    pval[0] = 0.0;
    for ( j = 0; j < b; j++ ) {
        if (testboot[j] > test[0]) pval[0] += 1.0;
    }
    pval[0] /= b;
    
    free(Hessianvec);
    free(Hessian);
    free(score);
}

SEXP _calculate_ST00(
    SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP X_, SEXP Z_,  
    SEXP ID_, SEXP Rvmatrix_, SEXP N, SEXP PZ, SEXP PX, SEXP B, SEXP K, SEXP TOL, SEXP SAVEALL
    ) {
    int n  = INTEGER(N)[0];
    int pz = INTEGER(PZ)[0];
    int px = INTEGER(PX)[0];
    int b  = INTEGER(B)[0];
    int k  = INTEGER(K)[0];
    int saveall  = INTEGER(SAVEALL)[0];
    // int p = pz + px;
    double tol = REAL(TOL)[0];
    double **x = vec2matrix(REAL(X_), n, px);
    double **z = vec2matrix(REAL(Z_), n, pz);
    int *yrmin = INTEGER(YRMIN_);
    int *yrmax = INTEGER(YRMAX_);

    SEXP TestR, TestB, Pval, list, list_names, TestRVEC, TestBVEC;
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(TestB       = allocVector(REALSXP,  b));
    PROTECT(Pval        = allocVector(REALSXP,  1));
    PROTECT(list_names  = allocVector(STRSXP,   5));
    PROTECT(list        = allocVector(VECSXP,   5));

    if (saveall) {
        PROTECT(TestRVEC = allocVector(REALSXP,  k));
        PROTECT(TestBVEC = allocVector(REALSXP,  k*b));
    } else {
        PROTECT(TestRVEC = allocVector(REALSXP,  1));
        PROTECT(TestBVEC = allocVector(REALSXP,  1));
    }

    // calculate H0
    double *BetaH0, *exb, *eta, *Vec, *vzMvec;
    BetaH0       = (double*)malloc(sizeof(double)*px);
    exb          = (double*)malloc(sizeof(double)*n);
    eta          = (double*)malloc(sizeof(double)*n);
    Vec          = (double*)malloc(sizeof(double)*n);
    vzMvec       = (double*)malloc(sizeof(double)*n*pz);
    double **vzM = vec2matrix(vzMvec, n, pz);
    EstCoxphNewton(BetaH0, x, INTEGER(Status_), yrmin, yrmax, n, px, tol);
    calculate_exb(BetaH0, x, n, px, exb);
    free(BetaH0);
    calculate_v4score(exb, INTEGER(Status_), yrmin, yrmax, eta, Vec, n);
    // calculate_ST_Related
    for ( int j = 0; j < pz; j++ ) {
        for ( int k = 0; k < n; k++ ) {
            vzM[k][j] = Vec[k] * z[k][j];
        }
    }
    free(exb);
    free(eta);
    free(Vec);

    double **rvmatrix = vec2matrix(REAL(Rvmatrix_), b, n);
    int **idM = ivec2matrix(INTEGER(ID_), k, n);
    calculate_ST00(
        z, idM, n, pz, b, k, rvmatrix, vzM, REAL(TestR), REAL(TestB), REAL(Pval),
        REAL(TestRVEC), REAL(TestBVEC), saveall
    );

    SET_STRING_ELT(list_names,  0,  mkChar("TestR"));
    SET_STRING_ELT(list_names,  1,  mkChar("TestB"));
    SET_STRING_ELT(list_names,  2,  mkChar("Pvalue"));
    SET_STRING_ELT(list_names,  3,  mkChar("TestRVEC"));
    SET_STRING_ELT(list_names,  4,  mkChar("TestBVEC"));
    SET_VECTOR_ELT(list,        0,  TestR);
    SET_VECTOR_ELT(list,        1,  TestB);
    SET_VECTOR_ELT(list,        2,  Pval);
    SET_VECTOR_ELT(list,        3,  TestRVEC);
    SET_VECTOR_ELT(list,        4,  TestBVEC);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(x);
    free(z);
    free(vzMvec);
    free(vzM);
    free(rvmatrix);
    free(idM);
    UNPROTECT(7);
    return list;
}


