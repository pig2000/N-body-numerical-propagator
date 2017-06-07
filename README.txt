Open source N-body numerical propagator

Developed by Weichen Xiao (Hong Kong University of Science and Technology), Chit Hong Yam (ispace inc.), Wen Han Chiu (Hong Kong University of Science and Technology).

Users should have the boost C++ library (http://www.boost.org/) and the NAIF-CSPICE (https://naif.jpl.nasa.gov/naif/toolkit_C.html) installed in prior.

Users should also put ephemeride files under the CSPICE directory.

Compile the code using the command "g++ propagator.cpp -Icspice/include cspice/lib/cspice.a".