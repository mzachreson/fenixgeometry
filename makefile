

bin:
	ifort -c fgtypes.f90
	ifort -c fgdata.f90
	ifort -c fgvolumecalc.f90
	ifort -c fgmaskcalc.f90
	ifort -c fgbuilder.f90
	ifort -c fgwrite.f90
	ifort -c fgmain.f90
	ifort -o fenixgeometry10 fgtypes.o fgdata.o fgvolumecalc.o fgmaskcalc.f90 fgbuilder.o fgwrite.o fgmain.o
	rm *.mod
	rm *.o
