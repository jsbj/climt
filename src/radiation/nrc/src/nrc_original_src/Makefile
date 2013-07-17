DEFRAD   = radinc.h
DEFLAY   = grid.h

driver:  driver.o setrad.o setspi.o optci.o sfluxi.o\
         gfluxi.o dsolver.o dtridgl.o setspv.o optcv.o\
         sfluxv.o gfluxv.o tpindex.o lagrange.o laginterp.o\
         radsetup.o dustprofile.o jsrchgt.o
	pgf90 -r8  -fast  -o driver driver.o\
               setrad.o setspi.o optci.o sfluxi.o\
               gfluxi.o dsolver.o dtridgl.o setspv.o optcv.o\
               sfluxv.o gfluxv.o tpindex.o lagrange.o laginterp.o\
               radsetup.o dustprofile.o jsrchgt.o

driver.o: driver.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast driver.f
setrad.o:  setrad.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast setrad.f
setspi.o:  setspi.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast setspi.f
setspv.o:  setspv.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast setspv.f
optci.o:   optci.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast optci.f
optcv.o:   optcv.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast optcv.f
sfluxi.o:  sfluxi.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast sfluxi.f
sfluxv.o:  sfluxv.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast sfluxv.f
gfluxi.o:  gfluxi.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast gfluxi.f
gfluxv.o:  gfluxv.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast gfluxv.f
dsolver.o: dsolver.f
	pgf90 -r8 -c -fast dsolver.f
dtridgl.o: dtridgl.f
	pgf90 -r8 -c -fast dtridgl.f
tpindex.o: tpindex.f $(DEFRAD) $(DEFLAY)
	pgf90 -r8 -c -fast tpindex.f
lagrange.o: lagrange.f
	pgf90 -r8 -c -fast lagrange.f
laginterp.o: laginterp.f $(DEFLAY) $(DEFRAD)
	pgf90 -r8 -c -byteswapio -fast laginterp.f
radsetup.o: radsetup.f
	pgf90 -r8 -c -byteswapio -fast radsetup.f
dustprofile.o: dustprofile.f
	pgf90 -r8 -c -fast dustprofile.f
jsrchgt.o: jsrchgt.f
	pgf90 -r8 -c -fast jsrchgt.f
