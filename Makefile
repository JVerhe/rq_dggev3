test:
	gfortran -c matrix_module.f90 
	gfortran -c dggev3_rq.f test_qz.f90
	gfortran test_qz.o dggev3_rq.o matrix_module.o -o test_qz -llapack -lblas

clean:
	find . -type f -executable -exec rm {} + && rm *.o