MKLROOT=/opt/intel/compilers_and_libraries_2019.0.117/linux/mkl

all: aout

aout: mod_redfield.o redfield.o
	ifort -o aout mod_redfield.o redfield.o -O2 -qopenmp -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include   ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

%.o: %.f90
	ifort -c $<

clean:
	rm *.o aout

