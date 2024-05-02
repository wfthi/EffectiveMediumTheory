driver_emt:driver_emt.o maxwell_garnett.o bruggeman.o mg_cde.o br_cde.o hunderi.o hu_cde.o 
	gfortran -o driver_emt driver_emt.o maxwell_garnett.o bruggeman.o mg_cde.o br_cde.o hunderi.o hu_cde.o

driver_emt.o:driver_emt.f
	gfortran -c -O driver_emt.f
maxwell_garnett.o:maxwell_garnett.f
	gfortran -c -O maxwell_garnett.f
bruggeman.o:bruggeman.f
	gfortran -c -O bruggeman.f
hunderi.o:hunderi.f
	gfortran -c -O hunderi.f
mg_cde.o:mg_cde.f
	gfortran -c -O mg_cde.f
br_cde.o:br_cde.f
	gfortran -c -O br_cde.f
hu_cde.o:hu_cde.f
	gfortran -c -O hu_cde.f
