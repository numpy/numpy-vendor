
########################################################################
#  
#  -- Automatically Tuned Linear Algebra Software (ATLAS)
#     (C) Copyright 2000 All Rights Reserved
#
########################################################################
#
#  -- Copyright notice and Licensing terms:
#
#  Redistribution  and  use in  source and binary forms, with or without
#  modification, are  permitted provided  that the following  conditions
#  are met:
#
#  1. Redistributions  of  source  code  must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce  the above copyright
#     notice,  this list of conditions, and the  following disclaimer in
#     the documentation and/or other materials provided with the distri-
#     bution.
#  3. The name of the University,  the ATLAS group,  or the names of its
#     contributors  may not be used to endorse or promote products deri-
#     ved from this software without specific written permission.
#
#  -- Disclaimer:
#
#  THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
#  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
#  CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#  TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
#  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
#  RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
#  CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
#  THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
########################################################################
#
include Make.inc
mySRCdir = $(SRCdir)/src/pthreads/blas/level3
#
########################################################################
#
INC3dep   = \
   $(INCAdir)/atlas_pthreads.h     $(INCSdir)/atlas_aux.h      \
   $(INCSdir)/atlas_enum.h      $(INCSdir)/atlas_misc.h     \
   $(INCSdir)/atlas_kern3.h     $(INCSdir)/atlas_rblas3.h   \
   $(INCSdir)/atlas_lvl3.h      $(INCSdir)/atlas_level3.h   \
   $(INCSdir)/atlas_ptmisc.h    $(INCSdir)/atlas_ptlvl3.h   \
   $(INCSdir)/atlas_ptlevel3.h
#
## Object files ########################################################
#
ATL_alvl3obj = \
   ATL_Sgemm_pt.o ATL_Ssymm_pt.o ATL_Ssyr2k_pt.o ATL_Ssyrk_pt.o \
   ATL_Strmm_pt.o ATL_Strsm_pt.o
#
ATL_slvl3obj = \
   ATL_sptgemm.o ATL_sptl3settype.o ATL_sptsymm.o ATL_sptsyr2k.o \
   ATL_sptsyrk.o ATL_spttrmm.o ATL_spttrsm.o
#
ATL_dlvl3obj = \
   ATL_dptgemm.o ATL_dptl3settype.o ATL_dptsymm.o ATL_dptsyr2k.o \
   ATL_dptsyrk.o ATL_dpttrmm.o ATL_dpttrsm.o
#
ATL_clvl3obj = \
   ATL_cptgemm.o ATL_cpthemm.o ATL_cpther2k.o ATL_cptherk.o \
   ATL_cptl3settype.o ATL_cptsymm.o ATL_cptsyr2k.o ATL_cptsyrk.o \
   ATL_cpttrmm.o ATL_cpttrsm.o
#
ATL_zlvl3obj = \
   ATL_zptgemm.o ATL_zpthemm.o ATL_zpther2k.o ATL_zptherk.o \
   ATL_zptl3settype.o ATL_zptsymm.o ATL_zptsyr2k.o ATL_zptsyrk.o \
   ATL_zpttrmm.o ATL_zpttrsm.o
#
## Targets #############################################################
#
all   : lib
#
lvl3  : lib
#
lib   : slib dlib clib zlib
#
slvl3 : slib
dlvl3 : dlib
clvl3 : clib
zlvl3 : zlib
#
slib  : slib.grd
dlib  : dlib.grd
clib  : clib.grd
zlib  : zlib.grd
#
slib.grd  : $(ATL_slvl3obj) $(ATL_alvl3obj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(ATL_slvl3obj) $(ATL_alvl3obj)
	$(RANLIB) $(ATLASlib)
	touch slib.grd
#
dlib.grd  : $(ATL_dlvl3obj) $(ATL_alvl3obj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(ATL_dlvl3obj) $(ATL_alvl3obj)
	$(RANLIB) $(ATLASlib)
	touch dlib.grd
#
clib.grd  : $(ATL_clvl3obj) $(ATL_alvl3obj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(ATL_clvl3obj) $(ATL_alvl3obj)
	$(RANLIB) $(ATLASlib)
	touch clib.grd
#
zlib.grd  : $(ATL_zlvl3obj) $(ATL_alvl3obj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(ATL_zlvl3obj) $(ATL_alvl3obj)
	$(RANLIB) $(ATLASlib)
	touch zlib.grd
#
########################################################################
#
ATL_Strsm_pt.o               : $(mySRCdir)/ATL_Strsm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS $(mySRCdir)/ATL_Strsm.c
ATL_Strmm_pt.o               : $(mySRCdir)/ATL_Strmm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS $(mySRCdir)/ATL_Strmm.c
ATL_Ssyr2k_pt.o              : $(mySRCdir)/ATL_Ssyr2k.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS $(mySRCdir)/ATL_Ssyr2k.c
ATL_Ssyrk_pt.o               : $(mySRCdir)/ATL_Ssyrk.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS $(mySRCdir)/ATL_Ssyrk.c
ATL_Ssymm_pt.o               : $(mySRCdir)/ATL_Ssymm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS $(mySRCdir)/ATL_Ssymm.c
ATL_Sgemm_pt.o               : $(mySRCdir)/ATL_Sgemm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS $(mySRCdir)/ATL_Sgemm.c
#
ATL_sptl3settype.o           : $(mySRCdir)/ATL_ptl3settype.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSREAL $(mySRCdir)/ATL_ptl3settype.c
ATL_sptgemm.o                : $(mySRCdir)/ATL_ptgemm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSREAL $(mySRCdir)/ATL_ptgemm.c
ATL_sptsymm.o                : $(mySRCdir)/ATL_ptsymm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSREAL $(mySRCdir)/ATL_ptsymm.c
ATL_sptsyrk.o                : $(mySRCdir)/ATL_ptsyrk.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSREAL $(mySRCdir)/ATL_ptsyrk.c
ATL_sptsyr2k.o               : $(mySRCdir)/ATL_ptsyr2k.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSREAL $(mySRCdir)/ATL_ptsyr2k.c
ATL_spttrmm.o                : $(mySRCdir)/ATL_pttrmm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSREAL $(mySRCdir)/ATL_pttrmm.c
ATL_spttrsm.o                : $(mySRCdir)/ATL_pttrsm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSREAL $(mySRCdir)/ATL_pttrsm.c
#
ATL_dptl3settype.o           : $(mySRCdir)/ATL_ptl3settype.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDREAL $(mySRCdir)/ATL_ptl3settype.c
ATL_dptgemm.o                : $(mySRCdir)/ATL_ptgemm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDREAL $(mySRCdir)/ATL_ptgemm.c
ATL_dptsymm.o                : $(mySRCdir)/ATL_ptsymm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDREAL $(mySRCdir)/ATL_ptsymm.c
ATL_dptsyrk.o                : $(mySRCdir)/ATL_ptsyrk.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDREAL $(mySRCdir)/ATL_ptsyrk.c
ATL_dptsyr2k.o               : $(mySRCdir)/ATL_ptsyr2k.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDREAL $(mySRCdir)/ATL_ptsyr2k.c
ATL_dpttrmm.o                : $(mySRCdir)/ATL_pttrmm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDREAL $(mySRCdir)/ATL_pttrmm.c
ATL_dpttrsm.o                : $(mySRCdir)/ATL_pttrsm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDREAL $(mySRCdir)/ATL_pttrsm.c
#
ATL_cptl3settype.o           : $(mySRCdir)/ATL_ptl3settype.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSCPLX $(mySRCdir)/ATL_ptl3settype.c
ATL_cptgemm.o                : $(mySRCdir)/ATL_ptgemm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSCPLX $(mySRCdir)/ATL_ptgemm.c
ATL_cpthemm.o                : $(mySRCdir)/ATL_pthemm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSCPLX $(mySRCdir)/ATL_pthemm.c
ATL_cptherk.o                : $(mySRCdir)/ATL_ptherk.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSCPLX $(mySRCdir)/ATL_ptherk.c
ATL_cpther2k.o               : $(mySRCdir)/ATL_pther2k.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSCPLX $(mySRCdir)/ATL_pther2k.c
ATL_cptsymm.o                : $(mySRCdir)/ATL_ptsymm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSCPLX $(mySRCdir)/ATL_ptsymm.c
ATL_cptsyrk.o                : $(mySRCdir)/ATL_ptsyrk.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSCPLX $(mySRCdir)/ATL_ptsyrk.c
ATL_cptsyr2k.o               : $(mySRCdir)/ATL_ptsyr2k.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSCPLX $(mySRCdir)/ATL_ptsyr2k.c
ATL_cpttrmm.o                : $(mySRCdir)/ATL_pttrmm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSCPLX $(mySRCdir)/ATL_pttrmm.c
ATL_cpttrsm.o                : $(mySRCdir)/ATL_pttrsm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DSCPLX $(mySRCdir)/ATL_pttrsm.c
#
ATL_zptl3settype.o           : $(mySRCdir)/ATL_ptl3settype.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDCPLX $(mySRCdir)/ATL_ptl3settype.c
ATL_zptgemm.o                : $(mySRCdir)/ATL_ptgemm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDCPLX $(mySRCdir)/ATL_ptgemm.c
ATL_zpthemm.o                : $(mySRCdir)/ATL_pthemm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDCPLX $(mySRCdir)/ATL_pthemm.c
ATL_zptherk.o                : $(mySRCdir)/ATL_ptherk.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDCPLX $(mySRCdir)/ATL_ptherk.c
ATL_zpther2k.o               : $(mySRCdir)/ATL_pther2k.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDCPLX $(mySRCdir)/ATL_pther2k.c
ATL_zptsymm.o                : $(mySRCdir)/ATL_ptsymm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDCPLX $(mySRCdir)/ATL_ptsymm.c
ATL_zptsyrk.o                : $(mySRCdir)/ATL_ptsyrk.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDCPLX $(mySRCdir)/ATL_ptsyrk.c
ATL_zptsyr2k.o               : $(mySRCdir)/ATL_ptsyr2k.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDCPLX $(mySRCdir)/ATL_ptsyr2k.c
ATL_zpttrmm.o                : $(mySRCdir)/ATL_pttrmm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDCPLX $(mySRCdir)/ATL_pttrmm.c
ATL_zpttrsm.o                : $(mySRCdir)/ATL_pttrsm.c $(INC3dep)
	$(ICC) -o $@ -c $(ICCFLAGS) -DATL_USEPTHREADS -DDCPLX $(mySRCdir)/ATL_pttrsm.c
#
########################################################################
#
sclean           :
	rm -f $(ATL_slvl3obj) $(ATL_alvl3obj) s*.grd a*.grd
#
dclean           :
	rm -f $(ATL_dlvl3obj) $(ATL_alvl3obj) d*.grd a*.grd
#
cclean           :
	rm -f $(ATL_clvl3obj) $(ATL_alvl3obj) c*.grd a*.grd
#
zclean           :
	rm -f $(ATL_zlvl3obj) $(ATL_alvl3obj) z*.grd a*.grd
#
#
clean            :
	rm -f *.o *.grd
#
########################################################################
