# ORCh
ORCh (Optimised and Reduced Chemistry) is a fully automated method to reduce detailed chemical schemes.

The development of ORCh is reported in two archival papers (make sure to cite those papers when using this program):

- N. Jaouen, L. Vervisch, P. Domingo, G. Ribert (2017) Automatic reduction and optimisation of chemistry for turbulent combustion modeling: Impact of the canonical problem, Combust. Flame, 175: 60-79.
- N. Jaouen, L. Vervisch, P. Domingo (2017) Auto-thermal reforming (ATR) of natural gas: An automated derivation of optimised reduced chemical schemes, Proc. Combust. Inst., 36(3): 3321-3330.

ORCh has been applied to various reactive flow problems:
- K. Bioche, G. Ribert, L. Vervisch (2019) Simulating upstream flame propagation in a narrow channel after wall preheating: Flame analysis and chemistry reduction strategy, Combust. Flame. 200: 219-231.
- C. Locci, L. Vervisch, B. Farcy, P. Domingo, N. Perret (2018) Selective Non-Catalytic Reduction (SNCR) of nitrogen oxide emissions: A perspective from numerical modeling, Flow Turbulence and Combust. 100(2): 301-340.
- K. Bioche, L. Vervisch, G. Ribert (2018) Premixed flame-wall interaction in a narrow channel: Impact of wall thermal conductivity and heat losses, J. Fluid Mech. 856: 5-35
- A. Bouaniche, N. Jaouen, P. Domingo, L. Vervisch (2019) Vitiated high Karlovitz n-decane/air turbulent flames: Scaling laws and micro-mixing modeling analysis, Flow Turbulence and Combust. 102(1): 235–252.
- K. Wan, L. Vervisch, C. Jianga, P. Domingo, Z. Gao, J. Xia, Z. Wang (2020) Development of reduced and optimized reaction mechanism for potassium emissions during biomass combustion based on genetic algorithms, Energy 211: 118565.

ORCh is a preprocessing tool designed to automatically generate and optimise reduced chemistry for specific CFD conditions, including fuel spray, heat losses, etc.

ORCh combines stochastic methods with graph-analysis together with genetic algorithms. Chemical kinetics are probed in laminar flames and turbulent micro-mixing (pairwise interaction or euclidean minimum spanning tree) canonical problems.

ORCh is coupled with HYPE (HYbrid stochastic sectional method for solving Population balance Equations) and machine learning to deal with carbon nanoparticles formation (nanotube and soot) and non-internal solid particles (crystallisation).

ORCh can include the training of ANN-CNN to reduce chemistry and pre-integrate stiff chemical systems or to solve for soot:

- K. Wan, C. Barnaud, L. Vervisch, P. Domingo (2020) Chemistry reduction using machine learning trained from non-premixed micro-mixing modeling: Application to DNS of a syngas turbulent oxy-flame with side-wall effects, Combust. Flame 220: 119-129.
- A. Seltz, P. Domingo, L. Vervisch (2021) Solving the population balance equation for non-inertial particles dynamics using PDF and neural networks: Application to a sooting flame, Phys. Fluids. 33, 013311.
- H.-T. Nguyen, C. Barnaud, P. Domingo, P.-D. Nguyen, L. Vervisch (2023) Large-Eddy Simulation of flameless combustion with neural-network driven chemistry, Application Energy Combust. Sci. 14:100126.

Most recent and on-going applications

Furnaces with Urea DeNOx
High-pressure partial oxydation of methane (hydrogen production)
Combustion in narrow channels with wall heat losses (micro- and meso-scale combustion, fire propagation between battery cells)
Heavy hydrocarbon liquid fuel (kerosene surrogate or pyrolyse)
Chlorofluorocarbon chemistry
Soot precursors
Coal combustion
Combustion of blast furnace gases



=== Pre installation for Cantera 2.4===


The 2.4 version of Cantera , downloaded with ORCh at workdir/orch/Cantera, needs :

* C++ Boost librairies (Version 1.68, the newer don't work with Cantera) to run properly. 
* SCONS compiler

==== SCONS ====

 https://sourceforge.net/projects/scons/files/scons/3.0.3/scons-3.0.3.tar.gz/download

 cd scons-3.0.3
 python setup.py install

==== BOOST LIBRARIES ====

https://www.boost.org/users/history/version_1_68_0.html

 cd boost_1_68_0 
 sh bootstrap.sh
 ./b2  install --libdir=yourPath/boost_1_68_0/lib   --includedir=<yourPath/boost_1_68_0/include

... patience ... the boost installation can be long

==== software used by Cantera 2.4 ====

If you don't already have Cantera 2.4 on your computer, you will need to get several external softwares (fmt, Eigen, googletest and sundials), located in orch/Cantera/ext/ : 

You can either download the ext.tar file and extract it in orch/Cantera/ext/ :

[[File:ext.tar.gz]]

Or directly download each software (could take some time ..)

* fmt

 $ git clone https://github.com/fmtlib/fmt.git
 $ sudo mkdir /usr/local/include/fmt
 $ sudo cp fmt/fmt/format.* /usr/local/include/fmt/
 
* Google test

 get https://github.com/google/googletest/archive/release-1.8.0.tar.gz

 tar xf release-1.8.0.tar.gz
 cd googletest-release-1.8.0
 cmake -DBUILD_SHARED_LIBS=ON .
 make

 then cp -r * ~/orch/Cantera/ext/googletest ...

* eigen

 visit https://eigen.tuxfamily.org/index.php?title=Main_Page
 Download the 3.2.10 tar
 tar -xvf 3.2.10.tar
 cp file into orch/Cantera/ext

* Sundials version 3.1



MPI librairies will be also needed for Stochastic configurations.

=== Installation procedure ===





== Installation ==
The following steps must be completed before you can use ORCh:
* Compile Cantera :
 In the directory workdir/orch/Cantera, modify the paths in the file "cantera.conf" with the appropriate one.
 build with SCONS : "scons build"
 If the build fail, for any reason, use "scons clean" before doing again the previous step, after coping with problems indicated in the output.

* Add to your .bashrc file:

 <code>export BOOSTPATH="/yourPath/boost_1_68_0"</code>

 <code>export GTCOMB_CT_HOME=/home/yourloggin/orch/Cantera</code>
 <code>export GTCOMB_CT_HOSTTYPE=$GTCOMB_CT_HOME/lib</code>
 <code>export GTCOMB_CT_DATA=$GTCOMB_CT_HOME/data</code>

* All the openmpi and gcc-compiler libraries must be available on your machine, typically you need something like below in your .bashrc file (some may already be there no need to duplicate):
 <code>source /opt/intel/composerxe/bin/compilervars.sh intel64</code>
 <code>export INTEL_HOME="/opt/intel/composerxe"</code>
 <code>export INTEL_INC="$INTEL_HOME/include"</code>
 <code>export INTEL_LIB="$INTEL_HOME/lib"</code>
 <code>export INTEL_BIN="$INTEL_HOME/bin"</code>
 <code>export INTEL_MAN="$INTEL_HOME/man"</code>
 <code>export PATH="$INTEL_BIN:$PATH"</code>
 <code>export LIBRARY_PATH="$INTEL_LIB:$LIBRARY_PATH"</code>
 <code>export LD_LIBRARY_PATH="$INTEL_LIB:$LD_LIBRARY_PATH"</code>
 <code>export MANPATH="$INTEL_MAN:$MANPATH"</code>

 <code>export MPI_HOME="/local/openmpi/intel-14.0.2/1.8.1"</code>
 <code>export MPI_INC="$MPI_HOME/include"</code>
 <code>export MPI_LIB="$MPI_HOME/lib"</code>
 <code>export MPI_BIN="$MPI_HOME/bin"</code>
 <code>export MPI_MAN="$MPI_HOME/share/man"</code>
 <code>export PATH="$MPI_BIN:$PATH"</code>
 <code>export LIBRARY_PATH="$MPI_LIB:$LIBRARY_PATH"</code>
 <code>export LD_LIBRARY_PATH="$MPI_LIB:$LD_LIBRARY_PATH"</code>
 <code>export MANPATH="$MPI_MAN:$MANPATH"</code>

 <code>export HDF5_HOME="/local/hdf5/intel-14.0.2/1.8.12"</code>
 <code>export HDF5_INC="$HDF5_HOME/include"</code>
 <code>export HDF5_BIN="$HDF5_HOME/bin"</code>
 <code>export HDF5_LIB="$HDF5_HOME/lib"</code>
 <code>export PATH="$HDF5_BIN:$PATH"</code>
 <code>export LIBRARY_PATH="$HDF5_LIB:$LIBRARY_PATH"</code>
 <code>export LD_LIBRARY_PATH="$HDF5_LIB:$LD_LIBRARY_PATH"</code>

 <code>export PAPI_HOME="/local/papi/intel-14.0.2/5.3.0"</code>
 <code>export LANG=C</code>
 <code>export LC_ALL=C</code>
 <code>export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu/c++/4.8</code>

* Check that the compilation of the ORCh package runs fine by testing one of the several test cases in orch/Tests/


 <code> make clean </code>
 <code> make </code>
 <code> ./mainProgram </code>

