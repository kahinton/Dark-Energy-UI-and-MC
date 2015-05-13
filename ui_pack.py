import numpy as np
import pylab as py

def Sing(x1,x2,y1,y2,c1,c2,u,n,Model):
    a = 10**(-9)
    dlna = 10**(-2)
    z = (1/a)-1
    cons = np.sqrt(6)/2 
    A = np.array(np.log10(a))
    Z = np.array(z)
    xi = np.array([x1, x2])
    yi = np.array([y1, y2])
    ci = np.array([c1, c2])
    Yi = (xi**2)/(yi**2)
    gi = eval(str(Model[1]))
    gipr = eval(str(Model[2]))
    gipr2 = eval(str(Model[3]))
    AYi = eval(str(Model[4]))
    px = gi+(Yi*gipr)
    Lam1 = (88.9*px[0]+.1)**.5
    Lam2 = (2.0**(.5))*px[1]-.1
    Lami = np.array([Lam1, Lam2])
    OmDEi = (xi**2)*(gi+2*Yi*gipr)
    OmegaDE = OmDEi[0] + n*OmDEi[1]
    OmegaRad = (u**2)
    OmegaM = 1-OmegaDE-OmegaRad
    wphi = (((xi[0]**2)*gi[0])+n*((xi[1]**2)*gi[1]))/(((xi[0]**2)*(gi[0]+2*Yi[0]*gipr[0]))+ n*((xi[1]**2)*(gi[1]+2*Yi[1]*gipr[1])))
    hdot = -(3/2) -(3/2)*((xi[0]**2)*gi[0]+n*((xi[1]**2)*gi[1]))-.5*(u**2)
    WPHI = np.array(wphi)
    OMat = np.array(OmegaM)
    ODE = np.array(OmegaDE)
    ORad = np.array(OmegaRad)
    delta = 10**(-3)
    deltap = .5
    deltap2 = ((3/2)*OmegaM*delta)-((hdot+2)*deltap)
    Delta = np.array(delta)
        
    while a < 1:
        f = 3*(((xi[0]**2)*gi[0])+ n*((xi[1]**2)*gi[1]))+(u**2)
        dxi = ((xi/2)*(3+f-2*cons*Lami*xi)+cons*AYi*(Lami*OmDEi-2*cons*xi*(gi+Yi*gipr)))*dlna
        dyi = ((yi/2)*(3+f-2*cons*Lami*xi))*dlna
        du = ((u/2)*(-1+f))*dlna
        xi += dxi
        yi += dyi
        u += du
        Yi = (xi**2)/(yi**2)
        gi = eval(str(Model[1]))
        gipr = eval(str(Model[2]))
        gipr2 = eval(str(Model[3]))
        AYi = eval(str(Model[4]))
        OmDEi = (xi**2)*(gi+2*Yi*gipr)
        OmegaDE = OmDEi[0] + n*OmDEi[1]
        OmegaRad = (u**2)
        OmegaM = 1-OmegaDE-OmegaRad
        wphi = ((((xi[0]**2)*gi[0])+n*((xi[1]**2)*gi[1]))/(((xi[0]**2)*(gi[0]+2*Yi[0]*gipr[0]))+ n*((xi[1]**2)*(gi[1]+2*Yi[1]*gipr[1]))))
        hdot = -(3/2) -(3/2)*((xi[0]**2)*gi[0]+n*((xi[1]**2)*gi[1]))-.5*(u**2)
        delta += deltap*dlna
        deltap += deltap2*dlna
        deltap2 = ((3/2)*OmegaM*delta)-((hdot+2)*deltap)
        Delta = np.append(Delta, delta)
        OMat = np.append(OMat, OmegaM)
        ODE = np.append(ODE, OmegaDE)
        ORad = np.append(ORad, OmegaRad)
        WPHI = np.append(WPHI, wphi)
        a = a*(1.0+dlna)
        z = (1.0/a)-1.0
        Z = np.append(Z,z)
        A = np.append(A,np.log10(a))
        
    # Defining Omega Plot
        
    D = Delta/(delta*(10.0**A))
    D = D/(np.max(D))
    
    return A, Z, OMat, ODE, ORad, WPHI, D
    
def Show(axis1, axis2, Check1, Check2, Check3, Check4, Check5, A, Z, OMat, ODE, ORad, WPHI, D, Model, n):
    if axis1 != axis2:
            X = A*axis1 + Z*axis2
            py.figure(1)
            ylab = ''
            if int(Check1) == 1:
                py.plot(X, OMat, 'y-', label = 'Omega_M', linewidth = 3)
                ylab = ylab+'Omega_M, '
            if int(Check2) == 1:
                py.plot(X, ODE, 'c-', label = 'Omega_DE', linewidth = 3)
                ylab = ylab+'Omega_DE, '
            if int(Check3) == 1:
                py.plot(X, ORad, 'm-', label = 'Omega_Rad', linewidth = 3)
                ylab = ylab+'Omega_Rad, '
            if int(Check4) == 1:
                py.plot(X, WPHI, 'k-', label = 'w_DE', linewidth = 3)
                ylab = ylab+'w_DE, '
            if int(Check5) == 1:
                py.plot(X, D, 'b-', label = 'D/a (normed)', linewidth = 3)
                ylab = ylab+'D/a'
            py.tick_params(length = 7, width = 3, labelsize = 'large', top = 'off' ,right = 'off')
            if axis1 == 1:
                py.xlim((-9,0))
                py.xlabel('Log10(a)',fontsize = 35)
                py.legend(loc = 3)
            if axis2 == 1:
                py.xlim((0,25))
                py.xlabel('z',fontsize = 35)
                py.legend(loc = 4)
            py.ylim((-1,1))
            py.title(Model[0]+' with '+str(n+1)+' fields', fontsize = 50)
            py.ylabel(ylab, fontsize = 35)                    
            py.show()
        
            return "-Succesfully Finished, Ready to Run Another Model-" 
    else:
            return "-Please Select One Axis-"