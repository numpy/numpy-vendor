/*
 *             Automatically Tuned Linear Algebra Software v3.8.3
 *                    (C) Copyright 1999 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef ATLAS_LAPACK_H
   #define ATLAS_LAPACK_H

#include "atlas_misc.h"
#include "cblas.h"

#ifdef PATL

#include "atlas_cblastypealias.h"
/*
 * predefined type macro names
 */
#define ATL_getriR   Mjoin(PATL,getriR)
#define ATL_getriC   Mjoin(PATL,getriC)
#define ATL_getri    Mjoin(PATL,getri)
#define ATL_lauumRL  Mjoin(PATL,lauumRL)
#define ATL_lauumRU  Mjoin(PATL,lauumRU)
#define ATL_lauumCL  Mjoin(PATL,lauumCL)
#define ATL_lauumCU  Mjoin(PATL,lauumCU)
#define ATL_lauum    Mjoin(PATL,lauum)
#define ATL_trtriRL  Mjoin(PATL,trtriRL)
#define ATL_trtriRU  Mjoin(PATL,trtriRU)
#define ATL_trtriCL  Mjoin(PATL,trtriCL)
#define ATL_trtriCU  Mjoin(PATL,trtriCU)
#define ATL_trtri    Mjoin(PATL,trtri)
#define ATL_potrfU   Mjoin(PATL,potrfU)
#define ATL_potrfL   Mjoin(PATL,potrfL)
#define ATL_potrs    Mjoin(PATL,potrs)
#define ATL_potrf    Mjoin(PATL,potrf)
#define ATL_getrfR   Mjoin(PATL,getrfR)
#define ATL_getrfC   Mjoin(PATL,getrfC)
#define ATL_getrs    Mjoin(PATL,getrs)
#define ATL_getrf    Mjoin(PATL,getrf)
#define ATL_laswp    Mjoin(PATL,laswp)

#endif

int ATL_sgetri(const enum CBLAS_ORDER Order, const int N, TYPE *A, const int lda,
               const int *ipiv, TYPE *wrk, int *lwrk);
int ATL_sgetriR(const int N, TYPE *A, const int lda, const int *ipiv,
                TYPE *wrk, const int lwrk);
int ATL_sgetriC(const int N, TYPE *A, const int lda, const int *ipiv,
                TYPE *wrk, const int lwrk);
void ATL_slauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, float *A, const int lda);
int ATL_spotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const int N, float *A, const int lda);
void ATL_spotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const int NRHS, const float *A, const int lda,
                float *B, const int ldb);
int ATL_sgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
               float *A, const int lda, int *ipiv);
void ATL_sgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
                const int N, const int NRHS, const float *A, const int lda,
                const int *ipiv, float *B, const int ldb);
void ATL_slaswp(const int N, float *A, const int lda0, const int K1,
               const int K2, const int *ipiv, const int inci);
int ATL_sgetrfC(const int M, const int N, float *A, const int lda,
                int *ipiv);
int ATL_sgetrfR(const int M, const int N, float *A, const int lda,
                int *ipiv);
void ATL_slauumRU(const int N, float *A, const int lda);
void ATL_slauumRL(const int N, float *A, const int lda);
void ATL_slauumCU(const int N, float *A, const int lda);
void ATL_slauumCL(const int N, float *A, const int lda);
int ATL_spotrfU(const int N, float *A, const int lda);
int ATL_spotrfL(const int N, float *A, const int lda);
int ATL_strtri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const enum CBLAS_DIAG Diag, const int N,
               float *A, const int lda);
int ATL_strtriRU(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_strtriRL(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_strtriCU(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_strtriCL(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);

int ATL_dgetri(const enum CBLAS_ORDER Order, const int N, TYPE *A, const int lda,
               const int *ipiv, TYPE *wrk, int *lwrk);
int ATL_dgetriR(const int N, TYPE *A, const int lda, const int *ipiv,
                TYPE *wrk, const int lwrk);
int ATL_dgetriC(const int N, TYPE *A, const int lda, const int *ipiv,
                TYPE *wrk, const int lwrk);
void ATL_dlauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, double *A, const int lda);
int ATL_dpotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const int N, double *A, const int lda);
void ATL_dpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const int NRHS, const double *A, const int lda,
                double *B, const int ldb);
int ATL_dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
               double *A, const int lda, int *ipiv);
void ATL_dgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
                const int N, const int NRHS, const double *A, const int lda,
                const int *ipiv, double *B, const int ldb);
void ATL_dlaswp(const int N, double *A, const int lda0, const int K1,
               const int K2, const int *ipiv, const int inci);
int ATL_dgetrfC(const int M, const int N, double *A, const int lda,
                int *ipiv);
int ATL_dgetrfR(const int M, const int N, double *A, const int lda,
                int *ipiv);
void ATL_dlauumRU(const int N, double *A, const int lda);
void ATL_dlauumRL(const int N, double *A, const int lda);
void ATL_dlauumCU(const int N, double *A, const int lda);
void ATL_dlauumCL(const int N, double *A, const int lda);
int ATL_dpotrfU(const int N, double *A, const int lda);
int ATL_dpotrfL(const int N, double *A, const int lda);
int ATL_dtrtri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const enum CBLAS_DIAG Diag, const int N,
               double *A, const int lda);
int ATL_dtrtriRU(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_dtrtriRL(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_dtrtriCU(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_dtrtriCL(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);

int ATL_cgetri(const enum CBLAS_ORDER Order, const int N, TYPE *A, const int lda,
               const int *ipiv, TYPE *wrk, int *lwrk);
int ATL_cgetriR(const int N, TYPE *A, const int lda, const int *ipiv,
                TYPE *wrk, const int lwrk);
int ATL_cgetriC(const int N, TYPE *A, const int lda, const int *ipiv,
                TYPE *wrk, const int lwrk);
void ATL_clauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, float *A, const int lda);
int ATL_cpotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const int N, float *A, const int lda);
void ATL_cpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const int NRHS, const float *A, const int lda,
                float *B, const int ldb);
int ATL_cgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
               float *A, const int lda, int *ipiv);
void ATL_cgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
                const int N, const int NRHS, const float *A, const int lda,
                const int *ipiv, float *B, const int ldb);
void ATL_claswp(const int N, float *A, const int lda0, const int K1,
               const int K2, const int *ipiv, const int inci);
int ATL_cgetrfC(const int M, const int N, float *A, const int lda,
                int *ipiv);
int ATL_cgetrfR(const int M, const int N, float *A, const int lda,
                int *ipiv);
void ATL_clauumRU(const int N, float *A, const int lda);
void ATL_clauumRL(const int N, float *A, const int lda);
void ATL_clauumCU(const int N, float *A, const int lda);
void ATL_clauumCL(const int N, float *A, const int lda);
int ATL_cpotrfRU(const int N, float *A, const int lda);
int ATL_cpotrfRL(const int N, float *A, const int lda);
int ATL_cpotrfU(const int N, float *A, const int lda);
int ATL_cpotrfL(const int N, float *A, const int lda);
int ATL_ctrtri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const enum CBLAS_DIAG Diag, const int N,
               float *A, const int lda);
int ATL_ctrtriRU(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_ctrtriRL(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_ctrtriCU(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_ctrtriCL(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);

int ATL_zgetri(const enum CBLAS_ORDER Order, const int N, TYPE *A, const int lda,
               const int *ipiv, TYPE *wrk, int *lwrk);
int ATL_zgetriR(const int N, TYPE *A, const int lda, const int *ipiv,
                TYPE *wrk, const int lwrk);
int ATL_zgetriC(const int N, TYPE *A, const int lda, const int *ipiv,
                TYPE *wrk, const int lwrk);
void ATL_zlauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, double *A, const int lda);
int ATL_zpotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const int N, double *A, const int lda);
void ATL_zpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const int NRHS, const double *A, const int lda,
                double *B, const int ldb);
int ATL_zgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
               double *A, const int lda, int *ipiv);
void ATL_zgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
                const int N, const int NRHS, const double *A, const int lda,
                const int *ipiv, double *B, const int ldb);
void ATL_zlaswp(const int N, double *A, const int lda0, const int K1,
               const int K2, const int *ipiv, const int inci);
int ATL_zgetrfC(const int M, const int N, double *A, const int lda,
                int *ipiv);
int ATL_zgetrfR(const int M, const int N, double *A, const int lda,
                int *ipiv);
void ATL_zlauumRU(const int N, double *A, const int lda);
void ATL_zlauumRL(const int N, double *A, const int lda);
void ATL_zlauumCU(const int N, double *A, const int lda);
void ATL_zlauumCL(const int N, double *A, const int lda);
int ATL_zpotrfRU(const int N, double *A, const int lda);
int ATL_zpotrfRL(const int N, double *A, const int lda);
int ATL_zpotrfU(const int N, double *A, const int lda);
int ATL_zpotrfL(const int N, double *A, const int lda);
int ATL_ztrtri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const enum CBLAS_DIAG Diag, const int N,
               double *A, const int lda);
int ATL_ztrtriRU(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_ztrtriRL(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_ztrtriCU(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_ztrtriCL(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);

#endif
