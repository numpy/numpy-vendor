#ifndef ATLAS_PTALIAS3_H
#define ATLAS_PTALIAS3_H
/*
 * Real BLAS
 */
   #define ATL_strsm   ATL_spttrsm
   #define ATL_strmm   ATL_spttrmm
   #define ATL_ssyr2k  ATL_sptsyr2k
   #define ATL_ssyrk   ATL_sptsyrk
   #define ATL_ssymm   ATL_sptsymm
   #define ATL_sgemm   ATL_sptgemm

   #define ATL_dtrsm   ATL_dpttrsm
   #define ATL_dtrmm   ATL_dpttrmm
   #define ATL_dsyr2k  ATL_dptsyr2k
   #define ATL_dsyrk   ATL_dptsyrk
   #define ATL_dsymm   ATL_dptsymm
   #define ATL_dgemm   ATL_dptgemm

/*
 * Complex BLAS
 */
   #define ATL_ctrsm     ATL_cpttrsm
   #define ATL_ctrmm     ATL_cpttrmm
   #define ATL_cher2k    ATL_cpther2k
   #define ATL_csyr2k    ATL_cptsyr2k
   #define ATL_cherk     ATL_cptherk
   #define ATL_csyrk     ATL_cptsyrk
   #define ATL_chemm     ATL_cpthemm
   #define ATL_csymm     ATL_cptsymm
   #define ATL_cgemm     ATL_cptgemm

   #define ATL_ztrsm     ATL_zpttrsm
   #define ATL_ztrmm     ATL_zpttrmm
   #define ATL_zher2k    ATL_zpther2k
   #define ATL_zsyr2k    ATL_zptsyr2k
   #define ATL_zherk     ATL_zptherk
   #define ATL_zsyrk     ATL_zptsyrk
   #define ATL_zhemm     ATL_zpthemm
   #define ATL_zsymm     ATL_zptsymm
   #define ATL_zgemm     ATL_zptgemm

#endif
