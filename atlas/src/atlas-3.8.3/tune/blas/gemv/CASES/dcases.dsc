10
  1  8  0  0 ATL_gemvN_mm.c     "R. Clint Whaley"
  2  0  1  1 ATL_gemvN_1x1_1.c  "R. Clint Whaley"
  3 16 32  1 ATL_gemvN_1x1_1a.c "R. Clint Whaley"
  4  0  4  2 ATL_gemvN_4x2_0.c  "R. Clint Whaley"
  5  0  4  4 ATL_gemvN_4x4_1.c  "R. Clint Whaley"
  6  0  8  4 ATL_gemvN_8x4_1.c  "R. Clint Whaley"
  7  0 16  2 ATL_gemvN_16x2_1.c "R. Clint Whaley"
  8  0 16  4 ATL_gemvN_16x4_1.c "R. Clint Whaley"
  9 16 32  4 ATL_gemvN_32x4_1.c "R. Clint Whaley"
 10 16 32  2 ATL_gemvN_SSE.c    "Camm Maguire" \
    gcc
    -fomit-frame-pointer -O2
 8
101 8  0  0 ATL_gemvT_mm.c      "R. Clint Whaley"
102 0  2  8 ATL_gemvT_2x8_0.c   "R. Clint Whaley"
103 0  4  8 ATL_gemvT_4x8_1.c   "R. Clint Whaley"
104 0  4 16 ATL_gemvT_4x16_1.c  "R. Clint Whaley"
105 0  2 16 ATL_gemvT_2x16_1.c  "R. Clint Whaley"
106 0  1  1 ATL_gemvT_1x1_1.c   "R. Clint Whaley"
107 0  2 32 ATL_gemvT_SSE.c     "Camm Maguire" \
    gcc
    -fomit-frame-pointer -O2
108 0  1 128 ATL_gemv_SSE.c     "Camm Maguire" \
    gcc
    -fomit-frame-pointer -O2

<flag> <Yunroll> <Xunroll> <filename> <author>
<flag> :  add to get all combinations of flags
          0: normal, ATLAS supplied
          8: GEMM-based primitive
         16: Outer-product or AXPY-based primitive (only valid for NoTrans)
         32: Do not block primitive for GEMV call
