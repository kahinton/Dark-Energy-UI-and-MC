from __future__ import division
import multiprocessing as mp
import os
import numpy as np
import random as rand
from scipy.interpolate import InterpolatedUnivariateSpline
from sympy import *

##################################################################################
##################################################################################
# To run the Monte Carlo Generator please go to the end of this code and enter all
# neccessary parameters. For a speed reference, 10,000 models split between four 
# processes at 2.9 GHz takes around 20 minutes to complete. Results will be visible
# in the provided UI. If you would like to use this code for some other project
# please include a reference to the original.
##################################################################################
##################################################################################

def Monte((N,x,n,model,loc,mins,maxes)):
    """Monte takes user input and runs many models over distributed processes"""
    zdat = np.genfromtxt(loc+'/Growth/PC_1234.dat','float', usecols = 0, skip_header =1)
    a1 = np.genfromtxt(loc+'/Growth/PC_1234.dat','float', usecols = 1, skip_header =1)
    a2 = np.genfromtxt(loc+'/Growth/PC_1234.dat','float', usecols = 2, skip_header =1)
    spl1 = InterpolatedUnivariateSpline(zdat ,a1, k=3)
    spl2 = InterpolatedUnivariateSpline(zdat, a2, k=3)
    Growth = open(loc+'/Growth/'+str(model[0])+'Growth%01d.dat'%x, 'w')
    Information = open(loc+'/Data/'+str(model[0])+'Information%1d.dat'%x, 'w')
    Start = 0
    while Start < N:        
        x1 = rand.uniform(mins[0],maxes[0])
        #y1 = rand.uniform(mins[1],maxes[1]) # Use this for other Models
        y1 = rand.uniform(.5,1.0)*x1  # Use this for Ghost Condensate
        x2 = rand.uniform(mins[2],maxes[2])
        #y2 = rand.uniform(mins[3],maxes[3]) # Use this for other Models
        y2 = rand.uniform(.5,1.0)*x2  # Use this for Ghost Condensate
        u = rand.uniform(mins[6],maxes[6])
        c1 = rand.uniform(mins[4],maxes[4])
        c2 = rand.uniform(mins[5],maxes[5])
        cons = (6**(1/2))/2 
        xi = np.array([x1, x2])
        yi = np.array([y1, y2])        
        ci = np.array([c1, c2])
        Yi = (xi**2)/(yi**2)
        gi = eval(model[1])
        gipr = eval(model[2])
        gipr2 = eval(model[3])
        px = gi+(Yi*gipr)
        Lam1 = (88.9*px[0]+.1)**.5
        Lam2 = (2.0**(.5))*px[1]-.1
        Lami = np.array([Lam1, Lam2])
        a = 10**(-9)
        dlna = 10**(-2)
        z = (1/a)-1
        AYi = eval(model[4])
        OmDEi = (xi**2)*(gi+2*Yi*gipr)
        OmegaDE = OmDEi[0] + n*OmDEi[1]
        OmegaRad = (u**2)
        OmegaM = 1-OmegaDE-OmegaRad
        wphi = (((xi[0]**2)*gi[0])+n*((xi[1]**2)*gi[1]))/(((xi[0]**2)*(gi[0]+2*Yi[0]*gipr[0]))+ n*((xi[1]**2)*(gi[1]+2*Yi[1]*gipr[1])))
        hdot = -(3/2) -(3/2)*((xi[0]**2)*gi[0]+n*((xi[1]**2)*gi[1]))-.5*(u**2)
        delta = 10**(-3)
        deltap = .5
        deltap2 = ((3/2)*OmegaM*delta)-((hdot+2)*deltap)
        Z = np.array(0)
        WZ = np.array(0)
        Delta = np.array(0)
        GrZ = np.array(0)

        try:
            while a < 1:
                f = 3*(((xi[0]**2)*gi[0])+ n*((xi[1]**2)*gi[1]))+(u**2)
                dxi = ((xi/2)*(3+f-2*cons*Lami*xi)+cons*AYi*(Lami*OmDEi-2*cons*xi*(gi+Yi*gipr)))*dlna
                dyi = ((yi/2)*(3+f-2*cons*Lami*xi))*dlna
                du = ((u/2)*(-1+f))*dlna
                xi += dxi
                yi += dyi
                u += du
                Yi = (xi**2)/(yi**2)
                gi = eval(model[1])
                gipr = eval(model[2])
                gipr2 = eval(model[3])
                AYi = eval(model[4])
                OmDEi = (xi**2)*(gi+2*Yi*gipr)
                OmegaDE = OmDEi[0] + n*OmDEi[1]
                OmegaRad = (u**2)
                OmegaM = 1-OmegaDE-OmegaRad
                wphi = (((xi[0]**2)*gi[0])+n*((xi[1]**2)*gi[1]))/(((xi[0]**2)*(gi[0]+2*Yi[0]*gipr[0]))+ n*((xi[1]**2)*(gi[1]+2*Yi[1]*gipr[1])))
                hdot = -(3/2) -(3/2)*((xi[0]**2)*gi[0]+n*((xi[1]**2)*gi[1]))-.5*(u**2)
                delta += deltap*dlna
                deltap += deltap2*dlna
                deltap2 = ((3/2)*OmegaM*delta)-((hdot+2)*deltap)
                z = (1/a)-1
                if z < 25.0:
                    Delta = np.append(Delta, delta)
                    GrZ = np.append(GrZ, z)
                if z < 3.0:
                    WZ = np.append(WZ, wphi)
                    Z = np.append(Z,z)
                a = a*(1+dlna)
    
        #Exporting Linear Growth Data
        
            GrZ = np.trim_zeros(GrZ, 'f')
            Delta = np.trim_zeros(Delta, 'f')
            GrZ = np.flipud(GrZ)
            Delta = (np.flipud(Delta))/delta       
            PickZ = np.array(0)
            PickD = np.array(0)
            
            for p in xrange(20):
                i = 16*p
                PickZ = np.append(PickZ, GrZ[i])
                PickD = np.append(PickD, Delta[i])
                
            PickZ = np.trim_zeros(PickZ, 'f')
            PickD = np.trim_zeros(PickD, 'f')/(1/(1+PickZ)) 
            Dline = 'D   '
            Zline = 'Z   '
            space = ' '
            width = 25
            
            for m in xrange(20):
                Dline = Dline + str(PickD[m])+(width - len(str(PickD[m])))*space
                Zline = Zline + str(PickZ[m])+(width - len(str(PickZ[m])))*space
            if x ==1:
                if Start == 0:
                    print >> Growth, Zline
            print >> Growth, Dline
    
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
            width = 22
            space = ' '
            checks = [w0, wa, wp, wphi, OmegaM, OmegaDE, x1, x2, y1, y2, c1, c2, u, Lam1, Lam2, n]
            line = ''
            for entry in checks:
                line += str(entry) + (width - len(str(entry)))*space
            if wp != 'nan':
                print >> Information, line
        except:
            pass
        Start +=  1
    Growth.close()
    Information.close()
    

def Cleaner((pros,amount,model,loc,fields)):
    """Cleaner takes the results from all Monte processes and places them in one file, as well as writing the neccessary information for viewing in the UI"""

    InfoFile = open(loc+'/Data/'+str(model[0])+'_'+str(fields)+'_fields.dat', 'w')
    GrowthFile = open(loc+'/Growth/'+str(model[0])+'Growth'+str(fields)+'.dat', 'w')
        
    width = 22
    space = ' '
    items = ['w_0', 'w_a', 'w_p', 'w_phi', 'Omega_M', 'Omega_DE', 'x1', 'x2', 'y1', 'y2', 'c1', 'c2', 'u', 'Lambda_1', 'Lambda_2', 'n']
    line = ''
    for entry in items:
        line += str(entry) + (width - len(str(entry)))*space
    print >> InfoFile, line
    
    for proc in range (1, pros+1):
        Growthinfo = open(loc+'/Growth/'+str(model[0])+'Growth%01d.dat'%proc, 'r')
        Info = open(loc+'/Data/'+str(model[0])+'Information%1d.dat'%proc, 'r')
        amount = sum(1 for line in Info)
        Info.seek(0)
        for n in range (0, amount):
            Results = Info.readline()
            InfoFile.write(Results)          
        if proc == 1:
            amount = amount + 1   
        for b in range (0, amount):
            Growth = Growthinfo.readline()
            GrowthFile.write(Growth)
        Info.close()
        Growthinfo.close()
    GrowthFile.close()
    InfoFile.close()
        
    for item in range(1, pros+1):
        os.remove(loc+'/Data/'+str(model[0])+'Information%1d.dat'%item)
        os.remove(loc+'/Growth/'+str(model[0])+'Growth%01d.dat'%item)
        

####################################################################################
####################################################################################
# Below all the information needed to run the Monte Carlo generator will be entered.
# Each value of Lambda will be automatically chosen in order to 
# Please note that growth data will only be correct for models which in which growth
# is unchanged from the standard GR solutions. If you know how growth should work in
# the model you want to test you should theoretically be able to change the function
# deltap2 in the Monte definition.
####################################################################################
####################################################################################

# Please name the model you wish to run

name = 'Ghost Condensate'

# Set the number of models to test. Note that Monte will only keep results that
# succeed, so the final number of models may be less than this number.

number = 10000

# The line below defines symbols needed for the function g, you should NOT change it
Yi,ci,yi,xi,u,Lami = symbols('Yi ci yi xi u Lami')

# Define the function g, for Quintessence this would be 1-(ci/Yi), for Ghost 
# Condensates the function is -1+(ci*Yi). To define your own please refer to 
# Tsujikawa and Ohashi

g = -1+(ci*Yi)

# Enter the number of late time fields used in each model

fields = 20

# Enter the ranges for each parameter used in your model. The ranges given in the
# UI are the ranges supplied here. The x parameters are analogous to the kinetic
# energy of the fields,the y parameters are analogous to the potential energy,
# the c parameters are the coefficients of the assumed exponential potentials, and
# lastly the u parameters define the initial density of radiation in the universe.
# Note that for the Ghost Condensate model you must change the values of y1 and y2 
# in the definition of the Monte function. The functions that need to be changed 
# are indicated by comments in the definition.

x1min,x1max = .001 , 0.3
y1min,y1max = .001 , 0.3
x2min,x2max = 10**(-14.0) , 10**(-13.0)
y2min,y2max = 10**(-14.0) , 10**(-13.0)
c1min,c1max = 1.0, 2.0
c2min,c2max = 1.0, 2.0
umin,umax = 0.7,0.9

# If you like you can change the number of processes to split the models between, but
# by default this will just be the number of available processors on your computer

pros = mp.cpu_count()

##################################################################################
##################################################################################
# Do not change anything beyond this point!!!!!
##################################################################################
##################################################################################

g2 = simplify(diff(g, Yi))
g3 = simplify(diff(g2, Yi))
Aterm = simplify((g+(5*Yi*g2)+(2*(Yi**2)*g3))**-1.0)
mins = [x1min,y1min,x2min,y2min,c1min,c2min,umin]
maxes = [x1max,y1max,x2max,y2max,c1max,c2max,umax]
mod = [name,str(g),str(g2),str(g3),str(Aterm)]
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
num = int(number/pros)

inputs = [[num,1,fields,mod,__location__,mins,maxes]]
for b in range(1,pros):
    inputs.append(([num, b+1,fields,mod,__location__,mins,maxes]))
    
clean = [(pros,num,mod,__location__,fields)]  
ParamFile = open(__location__+'/Parameters/'+name+'.dat', 'w')
ParamNeeds = [name, str(g), (str(g2)), (str(g3)), (str(Aterm)), x1min, x1max,x2min, x2max, y1min, y1max,y2min, y2max,c1min, c1max,c2min, c2max,
                    umin, umax, fields]
for entry in ParamNeeds:
    print >> ParamFile,str(entry)
ParamFile.close()

if __name__=='__main__':
    mp.freeze_support()
    pool = mp.Pool()
    pool.map(Monte, (inputs))
    pool.map(Cleaner, (clean))
    pool.close()
    pool.join()
