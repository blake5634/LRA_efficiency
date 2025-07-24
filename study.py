import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import subprocess as sub

import control as ctl
import LRAsim as LRA
import heatmap as hm


PULSES = 0
SINUSOID = 1

INPUT = SINUSOID

TRAPEZOIDAL = LRA.TRAPEZOIDAL  # numerical integration method for energy


# Set scan parameters

nwn =   11     # ODD number!!
nzet =  11
fres = 150 # hz
wres = 2*np.pi*fres

if nwn%2 != 1:
    error('Please use an odd number for nwn')

# data about this study
#
studyNameLegends = {'plot':'StudyPlot',
                    'leakage':'Energy Leakage (%)',
                    'eout':'Actuator Energy Output (J)',
                    'ebs':'Skin Damper Dissipation (J)',
                    'esk1':'Energy to Skin (J)',
                    'ee':'Energy Efficiency (%)',
                    'gain':'Amplitude Ratio (x3/x1)',
                    'fgain':'Skin amplitude / LRA force',
                    'emel':'Motor Electric Loss (J)',
                    'srms':'RMS Skin Displacement (mm)'
                        }


studykeys = studyNameLegends.keys()

def error(msg):
    print('\n\nError: ', msg)
    print('Usage: ')
    print('>python3 study [keyword]')
    print('keywords: ')
    for k in studykeys:
        print(f'   {k:15} : {studyNameLegends[k]:10}')
    quit()

sd={}  # study design dictionary

# get command line args:

if len(sys.argv) == 2:
    kwfound = False
    for k in studykeys:
            t1 = k.lower()
            t2 = sys.argv[1].lower()
            if t1.startswith(t2):
                sd['studytype'] = k
                sd['legend'] = studyNameLegends[k]
                kwfound = True
if not kwfound:
    error('unknown argument: ' + sys.argv[1])

#
#  Plot an existing result
#
if sd['studytype']=='plot':
#     fname = input('enter heatmap filename: ')

    #
    #   list the .csv heatmap files
    #
    cmd = 'ls heat*.csv'
    #rawres = sub.check_output(cmd,shell=True)
    try:
        rawres = sub.check_output(cmd,shell=True)
    except sub.CalledProcessError as grepexc:
        print("Sorry, there were no results")
        quit()
    # Parse the output as a list of strings
    lines = rawres.decode("utf-8").splitlines()

    i=0
    print("\n")
    for l in lines:
        i+=1
        print(f'{i:3}  {l}')

    try:
        fnum = int(input('\n Your choice: ').strip())
    except:
        error('Illegal choice: ', fnum)
    fname = lines[fnum - 1]

    FoundInName = False
    for k in studykeys:
        if k in fname:
            studytype = k
            FoundInName = True
    if FoundInName:
        hm.hmap(fname.strip(),studytype=studytype, legend=studyNameLegends[studytype])
    else:
        for k in studykeys:
            print(f'   {k:15} : {studyNameLegends[k]:10}')
        x = input(f'Enter Study Type for {fname}')
        try:
            hm.hmap(fname.strip(),studytype=x,legend=studyNameLegends[x])
        except:
            error('Did you mistype the studytype: ',x)
    quit()



#
#  Simulate and output a new result
#
# get the energy flows

def RepHeatmap(wn,z, fd, heat):
    print(f'{wn:10.2f}, {z:10.5f}, {heat:.3e}', file=fd)

def Report2(wn,z, fd, eso, eca, elm, eld, ebs, esk1, emel, Etot, yp):
    leakage = eso-(eld+ebs)
    plk = 100*leakage/eso
    print(f'{wn:10.2f}, {z:10.5f}, {eso:.3e}, {leakage:.3e}, {plk:8.1f}', file=fd)


wVals    = np.linspace(0.95*wres, 1.05*wres, num=nwn, endpoint=True)
zetaVals = np.linspace(0.01,            1.0, num=nzet, endpoint=True)

dt = LRA.dt

fhz = fres
per = 1/fhz
Tmax = 400*per  # let's model 400 cycles

T = np.arange(0,Tmax,dt)
npts = len(T)

ncontact = 4  # four fingers touching (x4 skin model)
Ain = 1.0  # input amplitude (N)


fname = f'heatmap_{sd['studytype']}_{nwn}x{nzet}.csv'
dataf = open(fname, 'w')

#
#   input force to LRA
#
if INPUT == SINUSOID:
    # Sin wave
    U = Ain*np.sin(wres*T)
elif INPUT == PULSES:
    # short rectangular pulses
    U = np.zeros(len(T))
    pf = 1/fres
    p  = int(pf/dt)

    for i,u in enumerate(U):
        if i%p==0:
            U[i] = Ain
        if (i+1)%p==0:
            U[i] = Ain
        if (i+2)%p==0:
            U[i] = Ain
else:
    error('Unknown Input Type Flag ("SINUSOID" or "PULSES")')

dpar = {}
#
#  for all the omega and zeta possibilities:
#
for w in wVals:
    for zetaLRA in zetaVals:
        print('.',end='',flush=True)
        #
        #   Simulate and analyze energy
        #

        # LRA properties
        Rm = 27.0  # Ohms
        Km = 1.0   # Newtons/Amp
        dpar['Km'] = Km
        dpar['Rm'] = Rm

        # M1 = 0.005 # kg   ()
        M1 = 0.0014 # kg (Lindsay et al 2013)
        dpar['M1'] = M1
        # B1 = 0.03222  # Nsec/m
        K1 = w**2 * M1
        dpar['K1'] = K1

        B1 = zetaLRA*2*np.sqrt(K1*M1)
        dpar['B1'] = B1

        # print(f' Derived K1 ratio: {K1/2800:5f} (ref: Ks= 2800)')
        # K1 = 2800  # N/m # (ref Jack's paper)

        # Case Properties
        Mc = 0.2250   # 225gr
        M2 = Mc
        dpar['M2']=M2

        # Skin Properties
        Kskin = 300 # N/m   (600-1200)
        K2=ncontact*Kskin     #
        dpar['K2'] = K2
        Bsk = (0.75+2.38)/2  # Nsec/m (0.75-2.38)
        # Bl = (2.38)  # Nsec/m (0.75-2.38)
        B3 = ncontact*Bsk
        dpar['B3']=B3

        Msk = 0.01 * M1  # essentially zero
        M3 = ncontact*Msk
        dpar['M3']= M3

        #
        # System Matrix
        a = -K1/M1
        b = -B1/M1
        c = K1/M1
        d = B1/M1

        e = K1/M2
        f = B1/M2
        g = (-K1-K2)/M2
        h = -B1/M2
        i = K2/M2

        j = K2/M3
        k = -K2/M3
        l = -B3/M3

        A = np.array([
            [0,1,0,0,0,0],
            [a,b,c,d,0,0],
            [0,0,0,1,0,0],
            [e,f,g,h,i,0],
            [0,0,0,0,0,1],
            [0,0,j,0,k,l]]  )

        B = np.array([
            [0],
            [-1/M1],
            [0],
            [1/M2],
            [0],
            [0]
            ] )
        C = np.identity(6)  # vector of all states
        D = np.zeros((6,1))           # no input coupling

        sdim = np.shape(A)[0]

        sys = ctl.ss(A,B,C,D)


        #
        #   Simulate the system
        #
        tp, yp = ctl.forced_response(sys, T, U)

        #
        #  RMS skin displ at resonance
        #
        if w == wres and zetaLRA == zetaVals[0]:  # at peak resonance
            _,_,x3RMS = LRA.RMSDisps(yp,U,dpar)
        #
        # compute the energy flows etc.
        #
        eso, eca, elm, eld, ebs, esk1, emel, Etot = LRA.EnergyFlows(yp,U,dpar)
        w_norm = w/wres
        leakage = eso-(eld+ebs)
        a1,a2,a3 = LRA.Amplitudes(yp,U,dpar)


        plk = 100*leakage/eso
        if sd['studytype'] == 'gain':  # skin deformation vs LRA mass deflection
            # print(f'Amplitudes: LRA: {1000*a1:f}mm   Skin: {1000*a2:f}mm')
            heatDataPt = a3/a1
        elif sd['studytype'] == 'fgain':
            heatDataPt = a3/np.max(U)    #  skin amplitude / LRA force
        elif sd['studytype'] == 'ee':
            heatDataPt = 100 * ebs/eso  # skin damper diss / actuator source output
        elif sd['studytype'] == 'leakage':
            heatDataPt = plk
        elif sd['studytype'] == 'eout':
            heatDataPt = eso
        elif sd['studytype'] == 'ebs':
            heatDataPt = ebs
        elif sd['studytype'] == 'esk1':
            heatDataPt = esk1
        elif sd['studytype'] == 'emel':
            heatDataPt = emel
        elif sd['studytype'] == 'srms':
            heatDataPt = 1000 * LRA.RMS(yp[4], LRA.timerange)
        else:
            error('Unknown study type: ',sd['studytype'])
        # store the heatmap square
        RepHeatmap(w_norm, zetaLRA, dataf, heatDataPt)

print('output map:', fname)
dataf.close()

print(f'Studytype: {sd["studytype"]}   legend: {sd["legend"]}')
print('Normalized Omega values: ', np.linspace(0.95, 1.05, num=nwn, endpoint=True)
)

print('RMS Skin displacement at resonance (m):',f'{x3RMS:.3e}')

hm.hmap(fname, studytype=sd['studytype'], legend=sd['legend'])

