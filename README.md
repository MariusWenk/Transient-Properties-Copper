# Density Function (DF) Program

This program is part of the result of my [bachelor's thesis](https://mariuswenk.github.io/assets/pdf/BSc_thesis.pdf) in physics at TU Kaiserslautern "Transient optical response and density-evolution of valence electrons in copper after ultrafast laser-excitation". It is based on a program by my supervisor Pascal D. Ndione.

The program predicts how ultrashort laser pulses affect the electronic and optical properties of copper based on the two-temperature model. After laser excitation, conduction electrons rapidly thermalize into a hot Fermi distribution, but band occupation remains out of equilibrium. Using a two-temperature model and rate equations for 4sp- and 3d-electrons, non-equilibrium electron dynamics are simulated. The resulting occupation and temperature data are used to compute the time-resolved dielectric function via the Drude-Lorentz formalism. These simulations help predict transient optical properties such as reflectivity and transmissivity, which can be compared with experimental data.
Studying the predictions, I found a high impact of initial dielectric function parameters and density of states on copper’s relaxation behavior, showing that assuming rapid electron thermalization provides a reasonable optical response model.

## Required libraries
* [GSL](https://www.gnu.org/software/gsl/doc/html/)
* [Boost](https://www.boost.org/)

## Setup
* To run this project, use the following steps:

* 1 $ make clean

* 2 $ make 

* 3 $ bin/df_program

## Python Code
Python scripts in PythonFiles can be used to format some of the input and output data.

Run: $ python3 *`Filename`*.py


### Author 
This project is based on work of Pascal D. Ndione for the two-temperature model in gold. This version was created by Marius Wenk.
Python Files with m_*`Filename`*.py are created by Marius Wenk, the remaining ones by Pascal Ndione.

### Contributors
AG Rethfeld - TU Kaiserslautern (RPTU)

### Contact
- [Marius_Wenk@gmx.de](Marius_Wenk@gmx.de)
- [pndione@protonmail.com](pndione@protonmail.com)

