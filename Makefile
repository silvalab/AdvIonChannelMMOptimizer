#compiler
CC = g++ 
#update the following line for your specifc Intel MKL installation for linking and compilation
EXTRA_LINK_LIBS = -L/export/intel/compilers_and_libraries_2017.6.256/linux/mkl/lib/intel64_lin/ -L/export/intel/compilers_and_libraries_2017.6.256/linux/compiler/lib/intel64_lin/    
#update the following line for Intel MKL headers
INCLUDE_DIRS = -IincludeH/ -I/export/intel/compilers_and_libraries_2017.6.256/linux/mkl/include/ 
#options to send to the compiler
CFLAGS = $(INCLUDE_DIRS) -std=c++11  -m64 -D USE_MKL -fopenmp
#other libraries needed for compilation
LIBS = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -liomp5 -ldl -lstdc++
#linker flags
LFLAGS = $(EXTRA_LINK_LIBS) $(LIBS)


all: Evalulator

Evalulator: main.o model.o graph.o math.o expm.o MarkovChannel.o restart.o modified_sobol.o setup.o cost.o SA.o
	$(CC) $(CFLAGS) -o Evalulator main.o model.o graph.o math.o expm.o MarkovChannel.o restart.o modified_sobol.o setup.o cost.o SA.o $(LFLAGS)    

main.o: src/main.cc includeH/MarkovChannel.hpp includeH/graph.hpp includeH/cost.hpp includeH/SA.hpp includeH/math.hpp includeH/helper.hpp includeH/restart.hpp 
	$(CC) -c $(CFLAGS) $(LFLAGS) src/main.cc 
	
model.o: src/model.cc includeH/model.hpp
	$(CC) -c $(CFLAGS) $(LFLAGS) src/model.cc 
	
graph.o: src/graph.cc includeH/graph.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS) src/graph.cc 
			
math.o: src/math.cc includeH/math.hpp includeH/modified_sobol.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS) src/math.cc 
	
MarkovChannel.o: src/MarkovChannel.cc includeH/MarkovChannel.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS) src/MarkovChannel.cc
	
cost.o: src/cost.cc includeH/cost.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS) src/cost.cc 

expm.o: src/expm.cc includeH/expm.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS) src/expm.cc 
		
SA.o: src/SA.cc includeH/SA.hpp includeH/restart.hpp
	$(CC) -c  $(CFLAGS) $(LFLAGS) src/SA.cc 	

restart.o: src/restart.cc includeH/restart.hpp 	
	$(CC) -c  $(CFLAGS) $(LFLAGS) src/restart.cc
	
modified_sobol.o: src/modified_sobol.cc includeH/modified_sobol.hpp
		$(CC) -c  $(CFLAGS) $(LFLAGS) src/modified_sobol.cc

setup.o: src/setup.cc includeH/setup.hpp
		$(CC) -c  $(CFLAGS) $(LFLAGS) src/setup.cc

.PHONY: clean cleanall
	
clean:
	rm -f *.o
	
cleanall: clean
	rm -f Evalulator

# DO NOT DELETE THIS LINE -- make depend needs it
