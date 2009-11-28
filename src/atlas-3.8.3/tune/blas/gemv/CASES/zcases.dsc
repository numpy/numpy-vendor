6
  1  8  0 0 ATL_cgemvN_mm.c      "R. Clint Whaley"
  2  0  1 1 ATL_cgemvN_1x1_1.c   "R. Clint Whaley"
  3 16 32 1 ATL_cgemvN_1x1_1a.c  "R. Clint Whaley"
  4  0  4 2 ATL_cgemvN_2x2_0.c   "R. Clint Whaley"
  5  0  4 4 ATL_cgemvN_4x2_1.c   "R. Clint Whaley"
  6 16 16 2 ATL_gemvN_SSE.c      "Camm Maguire" \
    gcc
    -fomit-frame-pointer -O2
6
101  8 0  0 ATL_cgemvT_mm.c      "R. Clint Whaley"
102  0 2  8 ATL_cgemvT_2x2_0.c   "R. Clint Whaley"
103  0 4  8 ATL_cgemvT_2x4_1.c   "R. Clint Whaley"
104  0 1  1 ATL_cgemvT_1x1_1.c   "R. Clint Whaley"
105  0 2 16 ATL_gemvT_SSE.c      "Camm Maguire" \
    gcc
    -fomit-frame-pointer -O2
106  0 1 128 ATL_gemv_SSE.c      "Camm Maguire" \
    gcc
    -fomit-frame-pointer -O2
