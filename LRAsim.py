
import numpy as np
import matplotlib.pyplot as plt
import sys
import control as ctl

#
#   Command line option
#
# if len(sys.argv) == 1:
#     RESONANT = True
# if len(sys.argv) == 2:
#     arg = sys.argv[1].lower()
#     print(f'got arg: ',arg)
#     if 'resonant'.startswith(arg):
#         RESONANT=True
#         zetaLRA = 0.01
#     else:
#         zetaLRA = 1.0
#         RESONANT=False

#
#   Compute energy flows
#

dt = 1.0e-4
TRAPEZOIDAL = True

def RMS(x,tr):
    idx = int((len(x)-1)*(1.0-tr)) # tr: pct of data used
    y = x[idx:]
    return np.sqrt(np.mean(x**2))


def integrate(x, xm1, dt):
    if TRAPEZOIDAL:
        v =  0.5*(x+xm1) * dt
        return v
    else:  # rectangular int.
        return x * dt

timerange = 0.25  # percent data for energy analysis (from end)


def Amplitudes(yp,U,d):
    # for forced_response
    x1  =  yp[0][:]   # MKS to mm
    xd1 =  yp[1][:]   # MKS to mm
    x2  =  yp[2][:]
    xd2 =  yp[3][:]
    x3  =  yp[4][:]
    xd3 =  yp[5][:]

    a1 = np.max(x1)
    a2 = np.max(x2)
    a3 = np.max(x3)

    return a1,a2,a3


def RMSDisps(yp,U,d):
    # for forced_response
    x1  =  yp[0][:]   # MKS to mm
    xd1 =  yp[1][:]   # MKS to mm
    x2  =  yp[2][:]
    xd2 =  yp[3][:]
    x3  =  yp[4][:]
    xd3 =  yp[5][:]

    a1 = np.sqrt(np.mean(x1**2))
    a2 = np.sqrt(np.mean(x2**2))
    a3 = np.sqrt(np.mean(x3**2))
    return a1,a2,a3



def EnergyFlows(yp,U,d):
    # for forced_response
    x1  =  yp[0][:]   # MKS to mm
    xd1 =  yp[1][:]   # MKS to mm
    x2  =  yp[2][:]
    xd2 =  yp[3][:]
    x3  =  yp[4][:]
    xd3 =  yp[5][:]

    # energy accumulators
    eso = 0.0  # source output
    eld = 0.0  # LRA dissipation
    eca = 0.0  # einput to case
    elm = 0.0  # einput to LRA mass
    ebs = 0.0  # skin damping dissipation
    esk1 = 0.0  # e output to skin
    emel = 0.0  # energy loss in motor resistance

    # delayed energy values for trapezoidal integration
    deld = {'eso':0, 'eld':0,'eca':0,'elm':0,'ebs':0,'esk1':0,'emel':0}



    start = int(len(x1)*(1-timerange))
    end = len(x1)-1
    indecesSelected = range(start,end,1)


    last_emel = 0.0

    k=-1
    for i in indecesSelected:
        k+=1

        # print('\n',k)
        # motor electrical loss (i^2R)
        f = U[i]
        im = (1/d['Km'])*f
        pm = im*im*d['Rm']
        emel += integrate(pm, last_emel, dt)
        last_emel = pm

        # energy from source (has two outputs)
        f = U[i]
        dx = xd2[i]*dt
        e1 = f*dx   # energy to M2
        f = -U[i]
        dx = xd1[i]*dt
        e2 = f*dx   # energy to M1
        # eso += f*dx
        eso += integrate(e1+e2, deld['eso'], dt)
        deld['eso'] =e1+e2

        # energy in LRA damper
        dx21 = (xd2[i]-xd1[i])
        f  = dx21*d['B1']   # damping force
        dx = dx21*dt   # length change
        # eld += f*dx
        eld += integrate(f*dx, deld['eld'], dt)
        deld['eld'] = f*dx

        # energy into LRA mass
        f = -U[i]
        dx = xd1[i] * dt
        # elm += f*dx
        elm += integrate(f*dx, deld['elm'], dt)
        deld['elm'] =   f*dx

        # energy to case
        f = U[i]
        dx = xd2[i]*dt
        # eca += f*dx
        eca += integrate(f*dx, deld['eca'], dt)
        deld['eca'] =    f*dx

        # dissipation in Bskin
        f = d['B3']*xd3[i]
        dx = xd3[i]*dt
        # ebs += f*dx
        ebs += integrate(f*dx, deld['ebs'], dt)
        deld['ebs'] =    f*dx

        # energy to skin
        f = d['K2']*(x2[i]-x3[i])
        dx = xd2[i]*dt
        # esk1 += f*dxs
        esk1 += integrate(f*dx, deld['esk1'], dt)
        deld['esk1'] =    f*dx


    # get energy in the state variables (masses, springs)

    Etot = d['K1']*(x2[-1]-x1[-1])**2 + d['K2']*(x2[-1]-x3[-1])**2 + d['M1']*xd1[-1]**2 + d['M2']*xd2[-1]**2 + d['M3']*xd3[-1]**2
    #
    # if TRAPEZOIDAL:
    #     intmeth = 'Trapezoidal'
    # else:
    #     intmeth = 'Rectangular'
    # print(f'EnergyFLows:\n   Energy integration: {intmeth}')
    # print(f'    (time range: last {100*timerange}%)')

    return (eso, eca, elm, eld, ebs, esk1, emel, Etot)



def LeakReport(eso, eca, elm, eld, ebs, esk1, emel, Etot,yp):
    # for forced_response
    x1  =  yp[0][:]   # MKS to mm
    x3  =  yp[4][:]
    #
    #    Reports
    #
    print(f'Oscillation Amplitudes (rms):  LRA: {RMS(1000*x1,timerange):.3e}mm Skin: {1000*RMS(x3,timerange):.3e}mm')

    print(f'        Source energy (eso): {eso:.3e}\nSource flows:')
    print(f'              to Case (eca): {eca:.3e}')
    print(f'          to LRA mass (elm): {elm:.3e}')
    print(f'                      total: {eca+elm:.3e}\n')

    print(f'Energy sinks:')
    print(f'motor resistive loss (emel): {emel:.3e}')
    print(f'      LRA dissipation (eld): {eld:.3e}')
    print(f'     skin dissipation (ebs): {ebs:.3e}')
    print(f'                      total: {eld+ebs:.3e}\n')
    print(f'      output to skin (esk1): {esk1}')
    print(f'   kinetic+potential (Etot): {Etot:.3e}')


    leakage = eso-(eld+ebs)
    print(f'              difference: {leakage:.3e}')
    print(f'          Energy leakage: ({100*(leakage)/(eso):.1f}%)')

