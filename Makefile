FC = gfortran
FFLAGS = -O3 -Wall -Wextra
LIBS = -llapack -lblas

# Module source files (.f90)
MODULE_SRC = modules/matrix_module.f90 modules/quicksort.f90

# Module object files should live in modules/
MODULE_OBJ = modules/matrix_module.o modules/quicksort.o

SRCS = dggev3_rq.f test_qz.f90 compare_methods.f90

all: test_qz compare_methods

# Compile modules → .o & .mod go into modules/
modules/%.o: modules/%.f90
	$(FC) $(FFLAGS) -c $< -J modules -o $@

# Compile regular sources → .o stays in working dir
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -I modules

%.o: %.f
	$(FC) $(FFLAGS) -c $< -I modules

# Executable: test_qz
test_qz: $(MODULE_OBJ) test_qz.o dggev3_rq.o
	$(FC) $(FFLAGS) $(MODULE_OBJ) test_qz.o dggev3_rq.o -o $@ $(LIBS)

# Executable: compare_methods
compare_methods: $(MODULE_OBJ) compare_methods.o dggev3_rq.o
	$(FC) $(FFLAGS) $(MODULE_OBJ) compare_methods.o dggev3_rq.o -o $@ $(LIBS)

clean:
	rm -f *.o *.mod modules/*.o modules/*.mod test_qz compare_methods
