

all: Main

Main: 
	mpic++ -w -std=c++11 -L$(MPI_LIB) -I$(BOOSTPATH)  -I$(GTCOMB_CT_HOME)/include -I$(HDF5_INC) drgep.cpp -c -g -O3 -Wall -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo -lctnumerics -lctmath -lconverters -lctbase -lctcxx -lcvode -lctf2c -lctf2c -lctlapack -lctblas -lm -lhdf5 -lhdf5_hl

clean:
	rm *.o
