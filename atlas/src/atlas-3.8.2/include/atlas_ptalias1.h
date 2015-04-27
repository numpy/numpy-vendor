#define ATLAS_PTALIAS1_H /* no threaded routs for Level 1 and 2 yet */
#ifndef ATLAS_PTALIAS1_H
#define ATLAS_PTALIAS1_H
/*
 * Real BLAS
 */
   #define ATL_dsdot   ATL_dsptdot
   #define ATL_sdsdot  ATL_sdsptdot
   #define ATL_sasum   ATL_sptasum
   #define ATL_snrm2   ATL_sptnrm2
   #define ATL_sdot    ATL_sptdot
   #define ATL_saxpy   ATL_sptaxpy
   #define ATL_scopy   ATL_sptcopy
   #define ATL_sscal   ATL_sptscal
   #define ATL_sswap   ATL_sptswap
   #define ATL_srotm   ATL_sptrotm
   #define ATL_srot    ATL_sptrot
   #define ATL_srotmg  ATL_sptrotmg
   #define ATL_srotg   ATL_sptrotg
   #define ATL_isamax  ATL_isptamax

   #define ATL_dasum   ATL_dptasum
   #define ATL_dnrm2   ATL_dptnrm2
   #define ATL_ddot    ATL_dptdot
   #define ATL_daxpy   ATL_dptaxpy
   #define ATL_dcopy   ATL_dptcopy
   #define ATL_dscal   ATL_dptscal
   #define ATL_dswap   ATL_dptswap
   #define ATL_drotm   ATL_dptrotm
   #define ATL_drot    ATL_dptrot
   #define ATL_drotmg  ATL_dptrotmg
   #define ATL_drotg   ATL_dptrotg
   #define ATL_idamax  ATL_idptamax

/*
 * Complex BLAS
 */
   #define ATL_cdotc_sub ATL_cptdotc_sub
   #define ATL_cdotu_sub ATL_cptdotu_sub
   #define ATL_caxpy     ATL_cptaxpy
   #define ATL_ccopy     ATL_cptcopy
   #define ATL_cscal     ATL_cptscal
   #define ATL_cswap     ATL_cptswap
   #define ATL_icamax    ATL_icptamax
   #define ATL_csscal    ATL_csptscal
   #define ATL_scnrm2    ATL_scptnrm2
   #define ATL_scasum    ATL_scptasum

   #define ATL_zdotc_sub ATL_zptdotc_sub
   #define ATL_zdotu_sub ATL_zptdotu_sub
   #define ATL_zaxpy     ATL_zptaxpy
   #define ATL_zcopy     ATL_zptcopy
   #define ATL_zscal     ATL_zptscal
   #define ATL_zswap     ATL_zptswap
   #define ATL_izamax    ATL_izptamax
   #define ATL_zdscal    ATL_zdptscal
   #define ATL_dznrm2    ATL_dzptnrm2
   #define ATL_dzasum    ATL_dzptasum

#endif
