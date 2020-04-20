CC = icc 
#CC = mbuild -v CC="icc -g"
EXTRA_LINK_LIBS = -L/export/intel/compilers_and_libraries_2017.6.256/linux/mkl/lib/intel64_lin/ -L/export/intel/compilers_and_libraries_2017.6.256/linux/compiler/lib/intel64_lin/    

INCLUDE_DIRS = -IincludeH/ -I/export/intel/compilers_and_libraries_2017.6.256/linux/mkl/include/ 

#EXTRA_INCLUDE = -I/usr/local/MATLAB/R2012a/extern/include/ 

CFLAGS = $(INCLUDE_DIRS) -std=c++11  -m64 -D USE_MKL -fopenmp

LFLAGS1 = $(EXTRA_LINK_LIBS)

LIBS = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -liomp5 -ldl -lstdc++

LFLAGS2 = $(LFLAGS1) $(LIBS)


all: Evalulator

Evalulator: main.o model.o graph.o math.o expm.o MarkovChannel.o restart.o modified_sobol.o setup.o cost.o SA.o
	$(CC) $(CFLAGS) -o Evalulator main.o model.o graph.o math.o expm.o MarkovChannel.o restart.o modified_sobol.o setup.o cost.o SA.o $(LFLAGS2)    

main.o: src/main.cc includeH/MarkovChannel.hpp includeH/graph.hpp includeH/cost.hpp includeH/SA.hpp includeH/math.hpp includeH/helper.hpp includeH/restart.hpp 
	$(CC) -c $(CFLAGS) $(LFLAGS2) src/main.cc 
	
model.o: src/model.cc includeH/model.hpp
	$(CC) -c $(CFLAGS) $(LFLAGS2) src/model.cc 
	
graph.o: src/graph.cc includeH/graph.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS2) src/graph.cc 
			
math.o: src/math.cc includeH/math.hpp includeH/modified_sobol.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS2) src/math.cc 
	
MarkovChannel.o: src/MarkovChannel.cc includeH/MarkovChannel.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS2) src/MarkovChannel.cc
	
cost.o: src/cost.cc includeH/cost.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS2) src/cost.cc 

expm.o: src/expm.cc includeH/expm.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS2) src/expm.cc 
		
SA.o: src/SA.cc includeH/SA.hpp includeH/restart.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS2) src/SA.cc 	

restart.o: src/restart.cc includeH/restart.hpp 	
	$(CC) -c  $(CFLAGS) $(LFLAGS2) src/restart.cc
	
modified_sobol.o: src/modified_sobol.cc includeH/modified_sobol.hpp
		$(CC) -c  $(CFLAGS) $(LFLAGS2) src/modified_sobol.cc

setup.o: src/setup.cc includeH/setup.hpp
		$(CC) -c  $(CFLAGS) $(LFLAGS2) src/setup.cc

.PHONY: clean cleanall
	
clean:
	rm -f *.o
	
cleanall: clean
	rm -f Evalulator

# DO NOT DELETE THIS LINE -- make depend needs it
