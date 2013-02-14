#Make file compatible with gcc version 3.4.4
#To compile and build the program, save makefile in directory above the source directory
#At the Linux command line:
# 1. Browse to the directory where the makefile is saved
# 2. Enter 'make compile build clean'
# g++ source/main.cpp smoking_sim.o sim_exception.o mersenne_class.o -o lbc_smokehist.exe 2> "out.txt"

compile:
	g++ -c -w source/smoking_sim.cpp source/sim_exception.cpp source/mersenne_class.cpp 2> "out.txt"

build:
	g++ source/main.cpp smoking_sim.o sim_exception.o mersenne_class.o -o lbc_smokehist.exe 2> "out.txt"

clean:
	\rm *.o 
	
error:
	cat out.txt
	
