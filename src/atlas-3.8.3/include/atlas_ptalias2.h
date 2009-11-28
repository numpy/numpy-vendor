#define ATLAS_PTALIAS2_H /* no threaded routs for Level 1 and 2 yet */
#ifndef ATLAS_PTALIAS2_H
#define ATLAS_PTALIAS2_H
/*
 * Real BLAS
 */
   #define ATL_sspr2   ATL_sptspr2
   #define ATL_ssyr2   ATL_sptsyr2
   #define ATL_sspr    ATL_sptspr
   #define ATL_ssyr    ATL_sptsyr
   #define ATL_sger    ATL_sptger
   #define ATL_stpsv   ATL_spttpsv
   #define ATL_stbsv   ATL_spttbsv
   #define ATL_strsv   ATL_spttrsv
   #define ATL_stpmv   ATL_spttpmv
   #define ATL_stbmv   ATL_spttbmv
   #define ATL_strmv   ATL_spttrmv
   #define ATL_sspmv   ATL_sptspmv
   #define ATL_ssbmv   ATL_sptsbmv
   #define ATL_ssymv   ATL_sptsymv
   #define ATL_sgbmv   ATL_sptgbmv
   #define ATL_sgemv   ATL_sptgemv

   #define ATL_dspr2   ATL_dptspr2
   #define ATL_dsyr2   ATL_dptsyr2
   #define ATL_dspr    ATL_dptspr
   #define ATL_dsyr    ATL_dptsyr
   #define ATL_dger    ATL_dptger
   #define ATL_dtpsv   ATL_dpttpsv
   #define ATL_dtbsv   ATL_dpttbsv
   #define ATL_dtrsv   ATL_dpttrsv
   #define ATL_dtpmv   ATL_dpttpmv
   #define ATL_dtbmv   ATL_dpttbmv
   #define ATL_dtrmv   ATL_dpttrmv
   #define ATL_dspmv   ATL_dptspmv
   #define ATL_dsbmv   ATL_dptsbmv
   #define ATL_dsymv   ATL_dptsymv
   #define ATL_dgbmv   ATL_dptgbmv
   #define ATL_dgemv   ATL_dptgemv

/*
 * Complex BLAS
 */
   #define ATL_chpr2     ATL_cpthpr2
   #define ATL_cher2     ATL_cpther2
   #define ATL_chpr      ATL_cpthpr
   #define ATL_cher      ATL_cpther
   #define ATL_cgerc     ATL_cptgerc
   #define ATL_cgeru     ATL_cptgeru
   #define ATL_ctpsv     ATL_cpttpsv
   #define ATL_ctbsv     ATL_cpttbsv
   #define ATL_ctrsv     ATL_cpttrsv
   #define ATL_ctpmv     ATL_cpttpmv
   #define ATL_ctbmv     ATL_cpttbmv
   #define ATL_ctrmv     ATL_cpttrmv
   #define ATL_chpmv     ATL_cpthpmv
   #define ATL_chbmv     ATL_cpthbmv
   #define ATL_chemv     ATL_cpthemv
   #define ATL_cgbmv     ATL_cptgbmv
   #define ATL_cgemv     ATL_cptgemv

   #define ATL_zhpr2     ATL_zpthpr2
   #define ATL_zher2     ATL_zpther2
   #define ATL_zhpr      ATL_zpthpr
   #define ATL_zher      ATL_zpther
   #define ATL_zgerc     ATL_zptgerc
   #define ATL_zgeru     ATL_zptgeru
   #define ATL_ztpsv     ATL_zpttpsv
   #define ATL_ztbsv     ATL_zpttbsv
   #define ATL_ztrsv     ATL_zpttrsv
   #define ATL_ztpmv     ATL_zpttpmv
   #define ATL_ztbmv     ATL_zpttbmv
   #define ATL_ztrmv     ATL_zpttrmv
   #define ATL_zhpmv     ATL_zpthpmv
   #define ATL_zhbmv     ATL_zpthbmv
   #define ATL_zhemv     ATL_zpthemv
   #define ATL_zgbmv     ATL_zptgbmv
   #define ATL_zgemv     ATL_zptgemv

#endif
