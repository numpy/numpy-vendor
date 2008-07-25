include Make.inc
mySRCdir = $(SRCdir)/src/blas/gemv

sobj = ATL_sgemv.o

dobj = ATL_dgemv.o

cobj = ATL_cgemv.o

zobj = ATL_zgemv.o

lib : all
all : slib dlib clib zlib
lvl2 : all
dlvl2 : dlib
slvl2 : slib
clvl2 : clib
zlvl2 : zlib
clean : sclean dclean cclean zclean

dlib : 
	$(PMAKE) -f Make_dmv dlib
	$(PMAKE) dlib.grd
dlib.grd : $(dobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(dobj)
	$(RANLIB) $(ATLASlib)
	touch dlib.grd
dclean : 
	rm -f $(dobj) d*.grd
slib : 
	$(PMAKE) -f Make_smv slib
	$(PMAKE) slib.grd
slib.grd : $(sobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(sobj)
	$(RANLIB) $(ATLASlib)
	touch slib.grd
sclean : 
	rm -f $(sobj) s*.grd
clib : 
	$(PMAKE) -f Make_cmv clib
	$(PMAKE) clib.grd
clib.grd : $(cobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(cobj)
	$(RANLIB) $(ATLASlib)
	touch clib.grd
cclean : 
	rm -f $(cobj) c*.grd
zlib : 
	$(PMAKE) -f Make_zmv zlib
	$(PMAKE) zlib.grd
zlib.grd : $(zobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(zobj)
	$(RANLIB) $(ATLASlib)
	touch zlib.grd
zclean : 
	rm -f $(zobj) z*.grd

ATL_sgemv.o : $(mySRCdir)/ATL_gemv.c $(INCAdir)/atlas_smv.h \
                   $(INCAdir)/atlas_smvN.h $(INCAdir)/atlas_smvT.h
	$(SKC) -o $@ -c $(sKCFLAGS) -DSREAL $(mySRCdir)/ATL_gemv.c
ATL_dgemv.o : $(mySRCdir)/ATL_gemv.c $(INCAdir)/atlas_dmv.h \
                   $(INCAdir)/atlas_dmvN.h $(INCAdir)/atlas_dmvT.h
	$(DKC) -o $@ -c $(dKCFLAGS) -DDREAL $(mySRCdir)/ATL_gemv.c

ATL_cgemv.o : $(mySRCdir)/ATL_gemv.c $(INCAdir)/atlas_cmv.h \
                   $(INCAdir)/atlas_cmvN.h $(INCAdir)/atlas_cmvT.h
	$(CKC) -o $@ -c $(cKCFLAGS) -DSCPLX $(mySRCdir)/ATL_gemv.c
ATL_zgemv.o : $(mySRCdir)/ATL_gemv.c $(INCAdir)/atlas_zmv.h \
                   $(INCAdir)/atlas_zmvN.h $(INCAdir)/atlas_zmvT.h
	$(ZKC) -o $@ -c $(zKCFLAGS) -DDCPLX $(mySRCdir)/ATL_gemv.c

