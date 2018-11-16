# OPTIMISTS
OPTIMISTS (Optimized PareTo Inverse Modeling through Integrated STochastic Search) is a data assimilation algorithm aimed at estimating the current conditions of a system based on a model and on recent measurements. It hybridizes concepts from established Bayesian and variational data assimilation methods, such as non-Gaussianity, non-linearity, non-sequential assimilation, and Pareto optimality, to allow for efficient estimation in complex/high-resolution systems.

OPTIMISTS' implementation can be coupled with any model whose state variables can be represented with a vector of real-coded values. The repository includes two specific case studies, where OPTIMISTS has been integrated with the VIC (Variable Infiltration Capacity)  modeling engine for large-scale hydrologic forecasting, and with the DHSVM (Distributed Hydrology Soil Vegetation Model) modeling engine for high-resolution hydrologic forecasting.

Citation: Hern·ndez, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779, https://doi.org/10.5194/hess-22-5759-2018, 2018. (Please cite this article if you use OPTIMISTS.)

OPTIMISTS was developed at the University of Pittsburgh by Felipe Hern√°ndez under the supervision of Xu Liang, with funding from the U.S. Department of Transportation (award OASRTRS-14-H-PIT) and through the William Kepler Whiteford Professorship of the University of Pittsburgh.

Copyright 2018 University of Pittsburgh

OPTIMISTS is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

OPTIMISTS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with OPTIMISTS. If not, see <https://www.gnu.org/licenses/>.
