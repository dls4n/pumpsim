#!/usr/bin/env python
#
# pumpsim: numerical model for electrogenic transport
#    Copyright (C) 2024  David Stokes 
#           stokes@nyu.edu
#           Dept. of Biochemistry and Molecular Pharmacology
#           NYU School of Medicine
#           550 First Ave
#           New York, NY 10016   USA
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#
#
# usage (with default parameters): pumpsim.py
# usage (with list of options): pumpsim.py -h
#
#
#python2.x compatibility
from __future__ import print_function

import platform
import numpy as np
#from scipy.optimize import curve_fit
#from scipy import stats
#from lmfit import Model, Parameters
import matplotlib.pyplot as plt
import sys

# declare globals and initialize parameters
global welcome, license   # help message
global dry_run
dry_run = False
global outfile    # output filename
outfile = ""
global plot
plot = True
# global parameters
global rate_input  # user specified max pumping rate
rate_input = 0.032     # default = 25/sec (msec is stepsize)
global vmaxKdp  # turnover rate specified by user
global GK    # K conductance specified by user input
GK = 0.02   # default in uS/cm**2
global ramp  # controls rate for solution exchange
ramp = 25.0    # 25 msec time constant for changing activity
global GNa
GNa = 0.0       # default in uS/cm**2
global Gssm     # SSM conductance
Gssm = 0.05    # default in uS/cm**2
global Vradius, Varea, Vvol, LPR, Kdp_dens
LPR = 2000    # lipid-to-protein ratio for vesicles
Vradius = 1e-5    # radius of a vescicle (cm)

#global ion concentrations
global ko, nao
ko  = 10.0;    # Bulk K (mM)
nao = 10.0;   # Bulk Na (mM)

#global data arrays (primary model paramters)
# these parameters evolve over time
global rate, time, Vm, k_i, na_i, rate_array
global rate3, tint1, tint2, tint3, tint4
global I_Kdp, I_Kleak, I_total
rate = 0.0   # rate for initial period
rate3 = 0.0  # rate for 3rd period
tint1 = 1000   # time for 1st period (msec)
tint2 = 3000   # time for 2nd period
tint3 = 3000   # time for 3rd period
tint4 = 1     # for valinomycin - off by default
rate_array = [0.0]
time = [0]
Vm = [0.0]
k_i = [10]
na_i = [10]
I_Kdp = [0]
I_Kleak = [0]
I_total = [0]
global rel_rate_array
rel_rate_array = [1]   #for debugging

# Physical constants
global Avogadro, R, T, F, Cm, Cssm
Avogadro = 6.023e23;  # Avogadro's num
R  = 8314.0;     # Gas constant (J/kmol/K)
T  = 310.0;      # Temperature (K)
F  = 96485.0;    # Faraday's constant (coul/mol)
Cm = 1.0;        # Specific membrane capacitance (uF/cm2)
Cssm = 0.37      # Specific SSM capacitance

# for voltage dependence of pump rate 
# linear dependence between Vmin and Vmax with Vslope and Vintercept
# exponential above and below these voltages using decay constants exp_min, exp_max
global Vmin, Vmax, Vslope, Vintercept
global exp_min, exp_max 
Vmin = -120.0   # these default values are in mV 
Vmax = -20.0    # from Gadsby 1989 paper on Na,K-ATPase
Vintercept = 1.0   # this is relative rate on y-axis
rmin = 0.25    # relative rate at Vmin - Gasdby 1989
rmax = 0.85    # relative rate at Vmax - Gadsby 1989
#Vslope = 0.006  # mV^-1  from Gadsby with Vmin=-120,Vmax=-20
Vslope = (rmin-Vintercept)/(Vmin-0)  # assume max turnover at 0 V
exp_min = 60.0  # default decay consts in mV 
exp_max = 25.0  # to fit sigmoidal curve from Gadsby data


##########################################################################
def ProcessCommandLine():
#    print sys.argv
    i = 0
    processing = False
    for c in sys.argv:
        if i:  # this eliminates command itself
            if not processing:
                if c == '-h' or c == '-help': 
                    ProcessCommandOption(c,None)
                    i += 1
                if c == '-n' or c == '-dry_run':
                    ProcessCommandOption(c,None)
                    i += 1
                else:
                    flag = c
                    arg = sys.argv[i+1]
                    ProcessCommandOption(flag,arg)
                    i += 1
                processing = True
            else: # skip to next cmd line item
                processing = False
                i += 1
        else:
            i += 1

def ProcessCommandOption(flag, arg):
    global welcome, license 
    global rate_input, outfile, dry_run, plot
    global rate3, tint1, tint2, tint3, tint4
    global vmaxKdp
    global Vmin, Vmax, Vslope, Vintercept
    global ramp  # controls rate of solution/activity exchange
    global GK
    global GNa
    global Gssm, Cssm
    global ko, nao
    global rate, time, Vm, k_i, na_i
    global I_Kdp, I_Kleak, I_total
    global Avogadro, R, T, F, Cm
    global Vradius, Varea, Vvol, LPR, Kdp_dens, Cm
    #print "ProcessCommandOption(%s,%s)" % (flag,arg)
    if flag == '-h' or flag == '-help':
        print (license)
        print (welcome)
        sys.exit(0)
    elif flag == '-r':
        rate_input = float(arg)/1000.
    elif flag == '-r3':
        rate3 = float(arg)/1000.
    elif flag == '-t1':
        tint1 = int(arg)
        if tint1 == 0: tint1 = 2
    elif flag == '-t2':
        tint2 = int(arg)
    elif flag == '-t3':
        tint3 = int(arg)
    elif flag == '-t4':
        tint4 = int(arg)
        if tint4 == 0: tint4 = 1
    elif flag == '-vmin':
        Vmin = float(arg)
    elif flag == '-vmax':
        Vmax = float(arg)
    elif flag == '-gk':
        GK = float(arg)
    elif flag == "-gssm":
        Gssm = float(arg)
    elif flag == "-cssm":
        Cssm = float(arg)
    elif flag == "-o":
        outfile = arg
    elif flag == "-tmix":
        ramp = float(arg)
    elif flag == "-ki":
        k_i = [float(arg)]
    elif flag == "-nai":
        na_i = [float(arg)]
    elif flag == "-ko":
        ko = float(arg)
    elif flag == "-nao":
        nao = float(arg)
    elif flag == "-lpr":
        LPR = float(arg)
    elif flag == "-rad":
        Vradius = float(arg)
    elif flag == "-cm":
        Cm = float(arg)
    elif flag == "-n" or flag == "-dry_run":
        dry_run = True
        print ("dry_run")
    elif flag == "-plot":
        if arg == '0' or arg == "False": plot = False
    else:
        print ("command line option not recognized: %s %s" % (flag,arg))
        print (welcome)
        sys.exit(1)


#########################################################################
def func_volt_dependence(v):
# for voltage dependence of pump rate 
# linear dependence between Vmin and Vmax with Vslope and Vintercept
# exponential above and below these voltages using decay constants exp_min, exp_max
    global Vmin, Vmax, Vslope, Vintercept
    global exp_min, exp_max
    global rel_rate_array    # for debugging
    
    if v <= Vmax and v >= Vmin:  # linear range
        rel_rate = v * Vslope + Vintercept
    elif v < Vmin:  # depolarized: exp decay
        Rmin = Vmin*Vslope + Vintercept   # rate at Vmin
        rel_rate = Rmin * np.exp((v-Vmin)/exp_min)
    elif v > Vmax:  # near zero volts: exp plateau
        Rmax = Vmax*Vslope + Vintercept  # rate at Vmax
        rel_rate = 1 - (Vintercept-Rmax) * np.exp((Vmax-v)/exp_max)
    rel_rate_array.append(rel_rate)
    return rel_rate   #relative rate between 0 and 1


##########################################################################
def func_kdp_ramp(t0,r0,t):   # t0,r0 are time an rate at outset of current interval
    global vmaxKdp
    global ramp  # controls rate at which pump reaches vmaxKdp
    global GK
    global GNa
    global ko, nao
    global rate, time, Vm, k_i, na_i
    global I_Kdp, I_Kleak, I_total
    global Avogadro, R, T, F, Cm
    global Varea, Vvol
    
    # current values
    v = Vm[-1]
    ki = k_i[-1]
    nai = na_i[-1]
    
    # Calculate equilibrium potentials
    EK   = 1*(R*T/F)*np.log(ko/ki) # mV
    ENa  = 1*(R*T/F)*np.log(nao/nai)
    
    # Calculate the leak
    # note that GK,GNa are in uS/cm2, must convert to mS/cm2
    IKleak = (v - EK)*GK/1000. # uA/cm2: mV * mS/cm2
    INaleak = (v - ENa)*GNa/1000.
    
    # Update rate based on ramp and vmaxKdp
    # this ramp is applied to account for slow response to solution exchange
    tramp = float(t-t0)   # time expired for this period
    if vmaxKdp > r0:
        rate = (vmaxKdp - r0)*(1-np.exp(-tramp/ramp)) + r0
    else:
        rate = (r0 - vmaxKdp)*np.exp(-tramp/ramp) + vmaxKdp
    #print("tramp, rate: %s %s\n" % (tramp,rate))
    #rate = vmaxKdp     # disable ramp
    
    # calculate IKdp
    rate_array.append(func_volt_dependence(v))
    IKdp = func_volt_dependence(v) * rate * Kdp_dens * F * 1e3 # A/cm2: (1/msec)(mol/cm2)(coul/mol)  
    #IKdp = vmaxKdp * Kdp_dens * F * 1e3; # A/cm2: troubleshoot w/ max activity
    IKdp = IKdp * 1e6 # uA/cm2
    Itotal = (IKdp + IKleak + INaleak)


    # change in membrane voltage (time step 1msec)
    #dv = -(1/Cm)* (INaK + Ileak);
    dv = -(1/Cm)* (IKdp + IKleak + INaleak) # (uA/cm2)/(uF/cm2): V/sec=mv/msec
        
    # change in concentrations: mM (for msec interval) (uA/cm2*cm2/F/ul = umoles/ul
    if ki > 0:
        dki = -1*IKdp*Varea/(F*Vvol) - IKleak*Varea/(F*Vvol) 
        #dki = -2*INaK*Varea/(F*Vvol)
    else:
        dki = 0
    
    if nai > 0:
        dnai = -1*INaleak*Varea/(F*Vvol)
    else:
        dnai = 0
        
    # Update global arrays
    time.append(t)
    Vm.append(v + dv)
    k_i.append(ki + dki)
    na_i.append(nai + dnai)
    I_Kdp.append(IKdp)
    I_Kleak.append(IKleak)
    I_total.append(Itotal)
    
    
    
    

#######################################################################
# Main Program

welcome = """usage: pumpsim.py
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
    """

license = """pumpsim: Copyright (C) 2024  David L. Stokes
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it
    under conditions specified in the GNU General Public License version 3
"""

ProcessCommandLine()

# Conductances for liposome (GK) and SSM (Gssm)
#GK = GK/1000.   # change units mS/cm**2 
#GNa = GNa/1000.
#Gssm = Gssm/1000.   # change units to mS/cm**2 

# Voltage dedendence of pump
Vslope = (rmin-Vintercept)/(Vmin-0)  # assume max turnover at 0 V
exp_min = 60.0*np.sqrt(0.006/Vslope)  # empirical decay const in mV
exp_max = 25.0   # to mimic sigmoidal shape from Gadsby 1989

print (license)
print ("max pump rate: %s sec-1, GK: %s uS/cm**2" % (rate_input*1000,GK))
print ('linear decrease of pump rate between: %s and %s mV' % (Vmax, Vmin))
print ("solution exchange time const: %s msec" % ramp)
print ("[K]in, [K]out, [Na]in, [Na]out: %s mM %s mM %s mM %s mM" % (k_i,ko,na_i,nao))
print ("LPR, radius: %s, %s cm" % (LPR, Vradius))
if dry_run: sys.exit(1)


# vescicle parameters
Varea = 4*np.pi*Vradius**2     # surface area of each vescicle (cm^2)
Vvol = 1333.3*np.pi*Vradius**3 # volume of each vescicle (uL)

# Calculate density of Kdp (per cm2)
# For 100 nm vesicle, there are 100 pumps/vesicle, each turns over at 32/s. 
# Note that these calculations are for only one leaflet of the bilayer
#   but pumps facing the wrong way will not be activated by ATP

lipid_dens = 1/63e-16/Avogadro #63 A^2= surf area of lipid mol: mol/cm2
Kdp_dens = lipid_dens/LPR # moles of pumps per cm2: mol/cm2


# First time interval
t0 = 0  # initial time
r0 = 0.0  # initial rate
vmaxKdp = 0.0   # initial period with no activity
tinterval = range(1,tint1)  # time steps in msec
for t in tinterval:
    func_kdp_ramp(t0,r0,t)

# Second time interval with Kdp active
t1 = tinterval[-1]   # continue with time trace
r1 = rate   # rate at end of last period
vmaxKdp = rate_input   # this period with fully active pump
tinterval = range(t1+1,t1+tint2)  
for t in tinterval:
    func_kdp_ramp(t1,r1,t)

# Third time interval with Kdp inhibited
t2 = tinterval[-1]   # continue with time trace
r2 = rate   # rate at end of last period
vmaxKdp = rate3   # 3rd period with inactivated pump
tinterval = range(t2+1,t2+tint3)  
for t in tinterval:
    func_kdp_ramp(t2,r2,t)
    
# Fourth time interval for adding valinomycin
t3 = tinterval[-1]   # continue with time trace
r3 = rate   # rate at end of last period
GK_in = GK
GK = 1000       # uS/cm2 - add valinomycin
vmaxKdp = 0   # switch off pump
tinterval = range(t3+1,t3+tint4)  
for t in tinterval:
    func_kdp_ramp(t3,r3,t)

# Calculate current measured by capacitive circuit
# eqn from Borlinghaus 1987 J. Membr. Biol. 97:161
# I_cap(t) = Cssm/(Cssm+Cm)[Itot(t) - (exp(t/tau)/tau)Integral(0-t)Itot(t)exp(t/tau)dt]
# tau = (Cssm + Cm)/(Gssm + GK) - f = SSM membrane, m = proteoliposome
tau = 1000*(Cssm + Cm)/(Gssm + GK_in + GNa)   # uF/uS => time constant of circuit in sec
I_cap = []
Icmin = []
Iint = []
integral = 0
coeff = Cssm/(Cssm + Cm)
#coeff = 1       # useful for overlaying the current traces
print ("tau = %s msec, scaling coeff = %s " % (tau,coeff))
tstep = 1.0    # tau in ms, so tstep also in msec
for t in range(len(time)):
    integral = integral + I_total[t]*np.exp(time[t]/tau)*tstep
    Iint.append(integral)
    Icmin.append(np.exp(-time[t]/tau)*integral/tau)
    I_cap.append(coeff*(I_total[t] - np.exp(-time[t]/tau)*integral/tau))

# Calculate the integrated charge and peak currents
QKdp = 0
Qtot = 0
Qtot1s = 0
Qcap = 0
Qcap1s = 0
tinterval = range(t1+1,t1+tint2)
for t in tinterval:
    QKdp += I_Kdp[t] * 0.001   # 1msec time intervals
    Qtot += I_total[t] * 0.001
    Qcap += I_cap[t] * 0.001
    if t-(t1+1) < 1000:
        Qtot1s += I_total[t] * 0.001
        Qcap1s += I_cap[t] * 0.001
Peak_Kdp = np.max(I_Kdp)
Peak_total = np.max(I_total)
Peak_cap = np.max(I_cap)
print('Qtot = %.4f, Qtot1s = %.4f, Qcap= %.4f, Qcap1s= %.4f, QKdp = %.4f uC/cm^2' 
          % (Qtot,Qtot1s,Qcap,Qcap1s,QKdp))
print('Peak currents: I_total= %.3f, I_cap= %.3f, I_pump= %.3f uA/cm^2'
          % (-Peak_total,-Peak_cap,-Peak_Kdp))

# output results
if outfile:
    fio = open(outfile,'w')
    for c in sys.argv:
        fio.write("%s " % c)    # print command line
    fio.write("\n")
    fio.write("max pump rate: %s sec-1, inhib pump rate: %s sec-1\n" 
         % (rate_input*1000,rate3*1000))
    fio.write("GK, Gssm: %s, %s uS/cm**2\n" % (GK_in,Gssm))
    fio.write('linear decrease of pump rate between: %s and %s mV\n' 
         % (Vmax, Vmin))
    fio.write("solution exchange time const: %s msec\n" % ramp)
    fio.write("[K]in, [K]out, [Na]in, [Na]out: %s mM %s mM %s mM %s mM\n" 
         % (k_i[0],ko,na_i[0],nao))
    fio.write("LPR, radius: %s, %s cm\n" % (LPR, Vradius))
    fio.write("Integrated charge: Q_tot= %.4f, Q_tot1s= %.4f, Q_cap= %.4f, "
           "Q_cap1s= %.4f, Q_pump= %.4f uC/cm^2\n"
        % (Qtot, Qtot1s, Qcap, Qcap1s, QKdp))
    fio.write('Peak currents: I_total= %.3f, I_cap= %.4f, I_pump= %.3f uA/cm^2\n'
          % (-Peak_total,-Peak_cap,-Peak_Kdp))
    fio.write("%s %s" % ('_'*80,'\n'))
    fio.write('%10s %10s %10s %10s %10s %10s %10s\n'
         % ('time', 'Vm', 'I_pump', 'I_total', 'I_capac', '[K]in', 'rel_rate'))
    fio.write("%s %s" % ('_'*80,'\n'))
    for i in range(len(time)):
        fio.write('%10.3f %10.3f %10.5f %10.5f %10.5f %10.3f %10.3f\n'
              % (time[i],Vm[i],-I_Kdp[i],-I_total[i],-I_cap[i],k_i[i],rate_array[i]))
    fio.close()
    
# plot data
# plot data for both exponential decay and leak
#plt.close()
if not plot: 
    print ("plotting disabled")
    sys.exit()

fig = plt.figure(figsize=(10,7), dpi=100)

# Plot I_total
res = fig.add_subplot(2,2,1)
text_font = {'size':'12', 'color':'black', 'weight':'normal'}
#res.set_title("%s" % filename)
#textstr = '$fit=Ae^{(-t/\\tau)^{\\beta}}+C$\n'
#textstr += '$A= %.2f,$ $C= %.2f$\n$\\tau= %.3f,$ ' % (At,Ct,Bt)
#textstr += '$\\beta= %.2f$\n$initrate= %.2e sec^{-1}$' % (Dt,initrate)
res.set_xlabel('time (msec)')
res.set_ylabel('current (uA/cm^2)')
I_plt_total = -1.0 * np.array(I_total)  #invert sign for K+ pumped out of vesicle
res.plot(time, I_plt_total, 'k-', label='I_total')
res.legend(loc="lower right", frameon=False)
#res.axis('tight')
res.autoscale(enable=True, axis='x', tight=True)
textstr = "turnover=%s/sec,  leak=%s uS/cm^2" % (rate_input*1000,GK_in)
res.text(0.5,1.1,textstr,transform=res.transAxes,verticalalignment='top', **text_font)


#Plot Vm
res2 = fig.add_subplot(2,2,2)
text_font = {'size':'16', 'color':'black', 'weight':'normal'}
res2.set_xlabel('time (msec)')
res2.set_ylabel('Vm (mV)')
res2.plot(time, Vm, 'k-', label='Vm')
res2.legend(loc="upper right",frameon=False)
#res.axis('tight')
res2.autoscale(enable=True, axis='x', tight=True)

res3 = fig.add_subplot(2,2,3)
text_font = {'size':'16', 'color':'black', 'weight':'normal'}
res3.set_xlabel('time (msec)')
res3.set_ylabel('current (uA/cm^2')
I_plt_cap = -1.0 * np.array(I_cap)
res3.plot(time, I_plt_cap, 'k-', label='I_cap')
res3.legend(loc="lower right", frameon=False)
#res.axis('tight')
res3.autoscale(enable=True, axis='x', tight=True)

res4 = fig.add_subplot(2,2,4)
text_font = {'size':'16', 'color':'black', 'weight':'normal'}
res4.set_xlabel('time (msec)')
res4.set_ylabel('[K] in vesicle (mM)')
res4.plot(time, k_i, 'k-', label='[K]in')
res4.legend(loc="upper right",frameon=False)
#res.axis('tight')
res4.autoscale(enable=True, axis='x', tight=True)

plt.show()

#try: dummy = raw_input("hit return to end")
#except NameError: input("\nhit return to end")





