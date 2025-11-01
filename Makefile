test:
	gfortran -c test_qz.f90
	gfortran -c dggev3_rq.f
	gfortran test_qz.o dggev3_rq.o -o test_qz -llapack -lblas

clean:
	find . -type f -executable -exec rm {} + && rm *.o