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
