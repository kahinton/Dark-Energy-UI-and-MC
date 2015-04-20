from __future__ import division
try:
    import numpy as np
except ImportError:
    print('This file requires the numpy package to run properly. Please see the readme for instructions on how to install this package.')
try:
    from scipy.interpolate import InterpolatedUnivariateSpline
except ImportError:
    print('This file requires the scipy package to run properly. Please see the readme for instructions on how to install this package.')
import os
import sys

# Creating the Splines
lok = os.path.dirname(sys.argv[0]) 
zdat = np.genfromtxt(lok+'/Growth/PC_1234.dat','float', usecols = 0, skip_header =1)
a1 = np.genfromtxt(lok+'/Growth/PC_1234.dat','float', usecols = 1, skip_header =1)
a2 = np.genfromtxt(lok+'/Growth/PC_1234.dat','float', usecols = 2, skip_header =1)
spl1 = InterpolatedUnivariateSpline(zdat ,a1, k=3)
spl2 = InterpolatedUnivariateSpline(zdat, a2, k=3)

def ParFile(loc, name, g, g2, g3, Aterm, x1min, x1max,x2min, x2max, y1min, y1max,y2min, y2max,c1min, c1max,c2min, c2max,umin, umax, fields):
    ParamFile = open(loc+'/Parameters/'+name+'.dat', 'w')
    ParamNeeds = [name, str(g), (str(g2)), (str(g3)), (str(Aterm)), x1min, x1max,x2min, x2max, y1min, y1max,y2min, y2max,c1min, c1max,c2min, c2max,
                    umin, umax, fields]
    for entry in ParamNeeds:
        print >> ParamFile,str(entry)
    ParamFile.close()
    
def ResFile(loc, name, fields, Tests):
    InfoFile = open(loc+'/Data/'+name+'_'+str(fields)+'_fields.dat', 'w')
    width = 22
    space = ' '
    items = ['w_0', 'w_a', 'w_p', 'w_phi', 'Omega_M', 'Omega_DE', 'x1', 'x2', 'y1', 'y2', 'c1', 'c2', 'u', 'Lambda_1', 'Lambda_2', 'n']
    line = ''
    for entry in items:
        line += str(entry) + (width - len(str(entry)))*space
    print >> InfoFile, line
    for Final in Tests:
        if Final.Outs.all() != 0:
            line2 = ''
            for value in Final.Outs:
                line2 += str(value) + (width - len(str(value)))*space
            print >> InfoFile, line2       
    InfoFile.close()
    
def GroFile(loc, name, fields, Tests):
    GrowthFile = open(loc+'/Growth/'+name+'Growth'+str(fields)+'.dat', 'w')
    width = 22
    space = ' '
    line = ''
    for entry in Tests[0].Grow:
        line += str(entry) + (width - len(str(entry)))*space
    print >> GrowthFile, line
    for Struct in Tests:
        if Struct.Grow.all() != 0:
            line2 = ''
            for value in Struct.Grow:
                line2 += str(value) + (width - len(str(value)))*space
            print >> GrowthFile, line2
    GrowthFile.close()
    

def Evolution((xi,yi,ci,u,Yi,gi,giS,gipr,giprS,gipr2,gipr2S,AYi,AYiS,n,Lami,OmDEi,OmegaDE,OmegaRad,OmegaM,wphi,hdot,delta,deltap,deltap2,spl1,spl2,Z,WZ,Delta,GrZ,w0,wa,wp)):
    """Evolution acts on a Universe class object and models it from a = 10**-9 to a = 1"""
    Fins = np.array([xi[0],xi[1],yi[0],yi[1],ci[0],ci[1],u,Lami[0],Lami[1],n])
    cons = (6.0**(1/2))/2.0
    a = 10**(-9)
    dlna = 10**(-2)
    z = (1/a)-1 
    while a < 1:
            f = 3*(((xi[0]**2)*gi[0])+ n*((xi[1]**2)*gi[1]))+(u**2)
            dxi = ((xi/2)*(3+f-2*cons*Lami*xi)+cons*AYi*(Lami*OmDEi-2*cons*xi*(gi+Yi*gipr)))*dlna
            dyi = ((yi/2)*(3+f-2*cons*Lami*xi))*dlna
            du = ((u/2)*(-1+f))*dlna
            xi += dxi
            yi += dyi
            u += du
            Yi = (xi**2)/(yi**2)
            gi = eval(giS)
            gipr = eval(giprS)
            gipr2 = eval(gipr2S)
            AYi = eval(AYiS)
            OmDEi = (xi**2)*(gi+2*Yi*gipr)
            OmegaDE = OmDEi[0] + n*OmDEi[1]
            OmegaRad = (u**2)
            OmegaM = 1-OmegaDE-OmegaRad
            wphi = ((((xi[0]**2)*gi[0])+n*((xi[1]**2)*gi[1]))/(((xi[0]**2)*(gi[0]+2*Yi[0]*gipr[0]))+ n*((xi[1]**2)*(gi[1]+2*Yi[1]*gipr[1]))))
            hdot = -(3/2) -(3/2)*((xi[0]**2)*gi[0]+n*((xi[1]**2)*gi[1]))-.5*(u**2)
            delta += deltap*dlna
            deltap += deltap2*dlna
            deltap2 = ((3/2)*OmegaM*delta)-((hdot+2)*deltap)
            
            if z < 25.0:
                Delta = np.append(Delta, delta)
                GrZ = np.append(GrZ, z)
            if z < 3.0:
                WZ = np.append(WZ, wphi)
                Z = np.append(Z,z)
            a = a*(1+dlna)
            z = (1/a)-1
        #Exporting Linear Growth Data
        
    GrZ = np.trim_zeros(GrZ, 'f')
    Delta = np.trim_zeros(Delta, 'f')
    GrZ = np.flipud(GrZ)
    Delta = (np.flipud(Delta))/delta       
    PickZ = np.zeros(20)
    PickD = np.zeros(20)
            
    for p in xrange(20):
            i = 16*p
            PickZ[p] = GrZ[i]
            PickD[p] = Delta[i]
                
    GrZ = PickZ
    Delta = (PickD)/(1/(1+PickZ)) 

    
    # Preparing Data for approximation
    
    Z = np.trim_zeros(Z, 'f')
    WZ = np.trim_zeros(WZ, 'f')
    
    a1s = spl1.__call__(Z)
    a2s = spl2.__call__(Z)
    norm1 = spl1.integral(0,3)
    norm2 = spl2.integral(0,3)
    con1 = 1/norm1
    con2 = 1/norm2
    
    #Defining Integrated Functions.
    
    Eta = (Z)/(1+Z)
    Int1 = np.flipud(con1*WZ*a1s)
    Int2 = np.flipud(con2*WZ*a2s)
    Int3 = np.flipud(con1*a1s*Eta)
    Int4 = np.flipud(con2*a2s*Eta)
    Z = np.flipud(Z)
    alpha1i = InterpolatedUnivariateSpline(Z, Int1, k=3)
    alpha1 = alpha1i.integral(0,3)
    alpha2i = InterpolatedUnivariateSpline(Z, Int2, k=3)
    alpha2 = alpha2i.integral(0,3)
    beta1i = InterpolatedUnivariateSpline(Z, Int3, k=3)
    beta1 = beta1i.integral(0,3)
    beta2i = InterpolatedUnivariateSpline(Z, Int4, k=3)
    beta2 = beta2i.integral(0,3)
    w0 = ((alpha2*beta1)-(alpha1*beta2))/(beta1-beta2)
    wa = (alpha1-alpha2)/(beta1-beta2)
    coef = (-.031)/(.093)
    wp = w0-coef*wa  
    if str(wp) != "nan": 
        END = np.array([w0,wa,wp,wphi,OmegaM,OmegaDE])
        END = np.append(END,Fins)
        GROWTH = np.array([GrZ,Delta])
        return END, GROWTH
    else:
        return np.array([0,0]),np.array([0,0])


class Universe:
    """The Universe Class contains all relevant data for each model universe"""
    yi,ci,u,Yi,gi,gipr,gipr2,Ayi,px,Lam1,n,giS,giprS,gipr2S,AYiS  = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    Lam2,Lami,OmDEi,OmegaDE,OmegaRad,OmegaM,wphi,hdot,delta,deltap,deltap2 = 0,0,0,0,0,0,0,0,0,0,0
    xi = np.zeros(2)
    spl1 = 0
    spl2 = 0
    Z = np.array(0)
    WZ = np.array(0)
    Delta = np.array(0)
    GrZ = np.array(0)
    w0 = 0
    wa = 0
    wp = 0
    def __init__(self,x1, x2, y1, y2, c1, c2, um, N, gi, gipr, gipr2, AYi,S1,S2):   
        self.xi = np.array([x1,x2])
        xi = self.xi
        self.yi = np.array([y1, y2])  
        yi = self.yi      
        self.ci = np.array([c1, c2])
        ci = self.ci
        self.u = um
        u = self.u
        self.n = N
        n = self.n
        self.Yi = (xi**2)/(yi**2)
        Yi = self.Yi
        self.giS = gi
        giS = self.giS
        self.gi = eval(self.giS)
        gi = self.gi
        self.giprS = gipr
        giprS = self.giprS
        self.gipr = eval(self.giprS)
        gipr = self.gipr
        self.gipr2S = gipr2
        gipr2S = self.gipr2S
        self.gipr2 = eval(self.gipr2S)
        gipr2 = self.gipr2
        self.AYiS = AYi
        AYiS = self.AYiS
        self.AYi = eval(self.AYiS)
        AYi = self.AYi
        self.px = gi+(Yi*gipr)
        px = self.px
        self.Lam1 = (88.9*px[0]+.1)**.5
        Lam1 = self.Lam1
        self.Lam2 = (2.0**(.5))*px[1]-.1
        Lam2 = self.Lam2
        self.Lami = np.array([Lam1,Lam2])
        Lami = self.Lami        
        self.OmDEi = (xi**2)*(gi+2*Yi*gipr)
        OmDEi = self.OmDEi
        self.OmegaDE = OmDEi[0] + n*OmDEi[1]
        OmegaDE = self.OmegaDE
        self.OmegaRad = (u**2)
        OmegaRad = self.OmegaRad
        self.OmegaM = 1-OmegaDE-OmegaRad
        OmegaM = self.OmegaM
        self.wphi = (((xi[0]**2)*gi[0])+n*((xi[1]**2)*gi[1]))/(((xi[0]**2)*(gi[0]+2*Yi[0]*gipr[0]))+ n*((xi[1]**2)*(gi[1]+2*Yi[1]*gipr[1])))
        wphi = self.wphi
        self.hdot = -(3/2) -(3/2)*((xi[0]**2)*gi[0]+n*((xi[1]**2)*gi[1]))-.5*(u**2)
        hdot = self.hdot
        self.delta = 10**(-3)
        delta = self.delta
        self.deltap = .5
        deltap = self.deltap
        self.deltap2 = ((3/2)*OmegaM*delta)-((hdot+2)*deltap)
        deltap2 = self.deltap2
        self.spl1 = S1
        spl1 = self.spl1
        self.spl2 = S2
        spl2 = self.spl2
    
    Outs, Grow = 0,0
    def Run(self,l):
        self.Outs,self.Grow = Evolution((self.xi,self.yi,self.ci,self.u,self.Yi,self.gi,self.giS,self.gipr,self.giprS,self.gipr2,self.gipr2S,self.Ayi,self.AYiS,self.n,self.Lami,
                                self.OmDEi,self.OmegaDE,self.OmegaRad,self.OmegaM,self.wphi,self.hdot,self.delta,self.deltap,self.deltap2,self.spl1,
                                self.spl2,self.Z,self.WZ,self.Delta,self.GrZ,self.w0,self.wa,self.wp))
        Outs = self.Outs
        Grow = self.Grow