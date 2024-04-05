PUMPSIM is a Python implementation of a numerical model for 
simulating data from in vitro assays of electrogenic ion 
transport.

The program comprises four time periods which mimic addition 
or removal of substrates during the experimental assays. 
After defining initial conditions, the program iteratively 
calls a subroutine to update time-dependent values at an 
interval of 1 ms.

Code was developed in Python using NumPy and MatPlotLib 
modules and has been tested with Python versions 2.7 and 3.8.


The program is initiated from the command line where 
various initial parameters can be controlled with the 
switches listed below, e.g., 
"pumpsim.py -rate 25 -gk 0.2 -o output.dat". 

Generation of plots can be controlled via the 
"-plot True/False" command line switch. 
Simulated data can be saved to a text-based output 
file via the "-o filename.dat" command line switch 
and plotted with external software, e.g., gnuplot 
or GraphPad Prism. A summary of command line options 
can be obtained with the "-h" switch. 

Command line options:
    -h: help
    -o: output filename, default=None
    -plot: plot data, default=True
    -r: max rate for pump (2nd time interval) [default=32 sec-1]
    -gk: leak conductance [default=0.02 uS/cm**2]
    -gssm: SSM conductance [default=0.05 uS/cm**2]
    -cssm: SSM capacitance (default=0.37 uF/cm**2)
    -r3: rate for inhibited pump (3rd time interval) [def. 0 sec-1]
    -t1: duration of 1st time interval: rate=0 [def. 0 msec]
    -t2: duration of 2nd time interval: pump on [default=3000]
    -t3: duration of 3rd time interval: pump inhib [default=3000 msec]
    -t4: duration of 4th time interval: add ionophore [default = 0 msec]
    -vmin: min Vm of linear voltage dependence for pump rate [default=-120 mV]
    -vmax: max Vm for linear voltage dependence for pump rate [default=-20 mV]
    -tmix: time constant for solution exchange/changing pump rate [default=25 msec]
    -ki: initial internal K concentration [default=10 mM]
    -nai: initial internal Na concentration [default=10 mM]
    -ko: external K conentration [default=10 mM]
    -nao: external Na concentration [default=10 mM]
    -lpr: molar lipid-to-protein ratio [default=2000]
    -rad: vesicle radius [default=1e-5 cm]
    -cm: vesicle membrane capacitance (default=1 uF/cm**2)
    

Copyright (C) 2024  David Stokes
   stokes@nyu.edu
   Dept. of Biochemistry and Molecular Pharmacology
   NYU School of Medicine
   550 First Ave
   New York, NY 10016   USA

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. 
<https://www.gnu.org/licenses/>.

