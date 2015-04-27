include ../Make.inc
arch=$(ARCH)
mach=$(ARCH)
dcc=default
cc=gcc
skc=$(skcD)
dkc=$(dkcD)
smc=$(smcD)
dmc=$(dmcD)
comps = skc=$(skc) dkc=$(dkc) smc=$(smc) dmc=$(dmc)

adefd = $(BLDdir)/ARCHS/$(mach)

tarfile : $(mach)
	- rm $(mach).tgz $(mach).tar
	tar cvf $(mach).tar $(mach)
	gzip --best $(mach).tar
	mv $(mach).tar.gz $(mach).tgz

$(mach)/kern/$(dcc) :
	mkdir $(mach)/kern/$(dcc)
$(mach)/gemm/$(dcc) :
	mkdir $(mach)/gemm/$(dcc)

$(mach) :
	mkdir $(mach)
	mkdir $(mach)/kern $(mach)/gemm
	$(MAKE) $(mach)/kern/$(skc) dcc=$(skc)
	$(MAKE) $(mach)/kern/$(dkc) dcc=$(dkc)
	$(MAKE) $(mach)/gemm/$(smc) dcc=$(smc)
	$(MAKE) $(mach)/gemm/$(dmc) dcc=$(dmc)

xnegflt : $(SRCdir)/CONFIG/ARCHS/negflt.c
	gcc -O -o xnegflt $(SRCdir)/CONFIG/ARCHS/negflt.c

ssys_put:
	- cp $(adefd)/kern/sMULADD $(SYSdir)/res/.
	- cp $(adefd)/kern/L1CacheSize $(SYSdir)/res/.
ssys :
	- cp $(SYSdir)/res/sMULADD $(adefd)/kern/.
	./xnegflt $(adefd)/kern/sMULADD
	- cp $(SYSdir)/res/L1CacheSize $(adefd)/kern/.
skern_put:
	- cp $(adefd)/kern/$(skc)/[s,c]TRSM_NB $(BLDdir)/tune/blas/level3/res/.
	- cp $(adefd)/kern/$(skc)/[s,c]MVRES $(BLDdir)/tune/blas/gemv/res/.
	- cp  $(adefd)/kern/$(skc)/[s,c]R1RES $(BLDdir)/tune/blas/ger/res/.
	- cp $(adefd)/kern/$(skc)/[s,c]*_SUMM $(L1Tdir)/res/.
skern: xnegflt
	rm -f $(adefd)/kern/$(skc)/[s,c]*
	- cp $(BLDdir)/tune/blas/level3/res/[s,c]TRSM_NB $(adefd)/kern/$(skc)/.
	- cp $(BLDdir)/tune/blas/gemv/res/[s,c]MVRES $(adefd)/kern/$(skc)/.
	- cp $(BLDdir)/tune/blas/ger/res/[s,c]R1RES $(adefd)/kern/$(skc)/.
	./xnegflt $(adefd)/kern/$(skc)/[s,c]R1RES \
                  $(adefd)/kern/$(skc)/[s,c]MVRES
	- cp $(L1Tdir)/res/[s,c]*_SUMM $(adefd)/kern/$(skc)/.
sgemm_put:
	- cp $(adefd)/gemm/$(smc)/[s,c]MMRES $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]gMMRES $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]Clean[M,N,K] $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]best[N,T][N,T]_* $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]uMMRES $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]uClean[M,N,K] $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]uClean[M,N,K]F $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]nreg $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]NB $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]NCNB $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]Xover.h $(INCAdir)/.
	- cp $(adefd)/gemm/$(smc)/[s,c]Xover.h $(MMTdir)/res/.
sgemm: xnegflt
	rm -f $(adefd)/gemm/$(smc)/[s,c]*
	- cp $(MMTdir)/res/[s,c]MMRES $(adefd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]gMMRES $(adefd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]Clean[M,N,K] $(adefd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]best[N,T][N,T]_* $(adefd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]uMMRES $(adefd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]uClean[M,N,K] $(adefd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]uClean[M,N,K]F $(adefd)/gemm/$(smc)/.
	./xnegflt $(adefd)/gemm/$(smc)/*
	- cp $(MMTdir)/res/[s,c]nreg $(adefd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]NB $(adefd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]NCNB $(adefd)/gemm/$(smc)/.
	- cp $(INCAdir)/[s,c]Xover.h $(adefd)/gemm/$(smc)/.
dsys_put:
	- cp $(adefd)/kern/dMULADD $(SYSdir)/res/.
	- cp $(adefd)/kern/L1CacheSize $(SYSdir)/res/.
dsys :
	- cp $(SYSdir)/res/dMULADD $(adefd)/kern/.
	./xnegflt $(adefd)/kern/dMULADD
	- cp $(SYSdir)/res/L1CacheSize $(adefd)/kern/.
dkern_put:
	- cp $(adefd)/kern/$(dkc)/[d,z]TRSM_NB $(BLDdir)/tune/blas/level3/res/.
	- cp $(adefd)/kern/$(dkc)/[d,z]MVRES $(BLDdir)/tune/blas/gemv/res/.
	- cp  $(adefd)/kern/$(dkc)/[d,z]R1RES $(BLDdir)/tune/blas/ger/res/.
	- cp $(adefd)/kern/$(dkc)/[d,z]*_SUMM $(L1Tdir)/res/.
dkern: xnegflt
	rm -f $(adefd)/kern/$(dkc)/[d,z]*
	- cp $(BLDdir)/tune/blas/level3/res/[d,z]TRSM_NB $(adefd)/kern/$(dkc)/.
	- cp $(BLDdir)/tune/blas/gemv/res/[d,z]MVRES $(adefd)/kern/$(dkc)/.
	- cp $(BLDdir)/tune/blas/ger/res/[d,z]R1RES $(adefd)/kern/$(dkc)/.
	./xnegflt $(adefd)/kern/$(dkc)/[d,z]R1RES \
                  $(adefd)/kern/$(dkc)/[d,z]MVRES
	- cp $(L1Tdir)/res/[d,z]*_SUMM $(adefd)/kern/$(dkc)/.
dgemm_put:
	- cp $(adefd)/gemm/$(dmc)/[d,z]MMRES $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]gMMRES $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]Clean[M,N,K] $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]best[N,T][N,T]_* $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]uMMRES $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]uClean[M,N,K] $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]uClean[M,N,K]F $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]nreg $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]NB $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]NCNB $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]Xover.h $(INCAdir)/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]Xover.h $(MMTdir)/res/.
dgemm: xnegflt
	rm -f $(adefd)/gemm/$(dmc)/[d,z]*
	- cp $(MMTdir)/res/[d,z]MMRES $(adefd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]gMMRES $(adefd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]Clean[M,N,K] $(adefd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]best[N,T][N,T]_* $(adefd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]uMMRES $(adefd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]uClean[M,N,K] $(adefd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]uClean[M,N,K]F $(adefd)/gemm/$(dmc)/.
	./xnegflt $(adefd)/gemm/$(dmc)/*
	- cp $(MMTdir)/res/[d,z]nreg $(adefd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]NB $(adefd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]NCNB $(adefd)/gemm/$(dmc)/.
	- cp $(INCAdir)/[d,z]Xover.h $(adefd)/gemm/$(dmc)/.
sys_put : ssys_put dsys_put
kern_put : skern_put dkern_put
gemm_put : sgemm_put dgemm_put
	- cp $(adefd)/gemm/atlas_cacheedge.h $(INCAdir)/.
	- cp $(adefd)/gemm/atlas_??NKB.h $(INCAdir)/.
sys : ssys dsys
kern : skern dkern
gemm : sgemm dgemm
	- cp $(INCAdir)/atlas_cacheedge.h $(adefd)/gemm/
	- cp $(INCAdir)/atlas_??NKB.h $(adefd)/gemm/
ArchNew : $(mach) xnegflt
	$(MAKE) sys kern gemm
	- cp $(BLDdir)/bin/INSTALL_LOG/?PerfSumm.txt $(adefd)/.
	rm -f xnegflt
archput : sys_put kern_put gemm_put
