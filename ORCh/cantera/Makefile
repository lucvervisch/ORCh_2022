

all: Main

Main: 
	mpic++ -w -std=c++11 -L$(MPI_LIB) -I$(BOOSTPATH) -L$(BOOSTPATH) -I$(GTCOMB_CT_HOME)/include -L$(GTCOMB_CT_HOSTTYPE) -I$(HDF5_INC) -I$(HDF5_INC) -L$(HDF5_LIB) -L$(TENSORFLOW_CAPI_PATH)/lib -I$(TENSORFLOW_CAPI_PATH)/include -L$(OPENCV_PATH)/lib64 -I$(OPENCV_PATH)/include computeMultipleInlet.cpp computePremixedFlames.cpp computeAutoIgnition.cpp Analytic_function.cpp particle.cpp flamemodel.cpp -c -g -O3 -Wall -lcantera -lm -lhdf5 -lhdf5_hl -ltensorflow -ltensorflow_framework
clean:
	rm *.o
