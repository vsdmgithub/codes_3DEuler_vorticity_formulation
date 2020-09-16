# MAKEFILE FOR EULER

# DEFINE VARIABLES
# ---------------------------start-----
prog=euler_3D.f90
out=run_and_output.f90
out_ob=run_and_output.o
sol=solver.f90
sol_ob=solver.o
var=variables_and_arrays.f90
var_ob=variables_and_arrays.o
fft=FFT_mod.f90
fft_ob=FFT_mod.o
vtr=VTR_mod.f90
vtr_ob=VTR_mod.o
sts=STAT_mod.f90
sts_ob=STAT_mod.o
cc=gfortran
cc_lc= -I/usr/local/include
lb_fftw=-L/usr/local/include -lfftw3 -lm
run=./ex
#----------------------------end-------


# MAKEFILE
# ---------------------------start----- 
ex:$(ob)
	$(cc) $(cc_lc) -c $(vtr) 
	$(cc) $(cc_lc) -c $(sts) 
	$(cc) $(cc_lc) -c $(fft) $(lb_fftw) 
	$(cc) $(cc_lc) -c $(var) 
	$(cc) $(cc_lc) -c $(sol) 
	$(cc) $(cc_lc) -c $(out) 
	$(cc) $(cc_lc) -c $(prog) 
	$(cc) $(cc_lc) $(prog) $(var_ob) $(sol_ob) $(out_ob) $(vtr_ob) $(sts_ob) $(fft_ob) $(lb_fftw) -o ex 
	$(run)
#----------------------------end-------

# CLEANING
# ---------------------------start----- 
clean:
	rm ex
	rm *.mod
	rm *.o
cl:
	rm *.mod
	rm *.o
#----------------------------end-------
