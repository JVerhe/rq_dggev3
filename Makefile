test_qz:
	gfortran -c matrix_module.f90 dggev3_rq.f test_qz.f90
	gfortran test_qz.o dggev3_rq.o matrix_module.o -o test_qz -llapack -lblas

pencils:
	gfortran -c matrix_module.f90 quicksort.f90 dggev3_rq.f pencils.f90
	gfortran pencils.o dggev3_rq.o matrix_module.o quicksort.o -o pencils -llapack -lblas

clean:
	find . -type f -executable -exec rm {} + && rm *.o