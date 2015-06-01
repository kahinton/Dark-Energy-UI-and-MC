from __future__ import division
import multiprocessing as mp
import os
import sys
from random import uniform
try:
    import numpy as np
except ImportError:
    print('This file requires the numpy package to run properly. Please see the readme for instructions on how to install this package.')
try:
    from sympy import symbols, simplify, diff
except ImportError:
    print('This file requires the sympy package to run properly. Please see the readme for instructions on how to install this package.')
from monte_pack import Universe, e1, e2, ParFile, ResFile, GroFile


###################################################################################
# Welcome to DE_mod_gen. This file allows you to create and test your own models
# of DE in the framework of a multi-field effective field theory. If you are 
# interested in the functions being used they are contained in the universe.py
# file.
# Below all the information needed to run the Monte Carlo generator can be entered.
# Each value of lambda is calculated to minimize the number of failed models.
# Please note that growth data will only be correct for models which in which growth
# is unchanged from the standard GR solutions. 
####################################################################################

# Please name the model you wish to run. This name will also be used in the UI.

name = 'Test_ws'

# Set the number of models to test. Note that only models which produce valid
# results will be presented in the UI, so the final number may be less than this
# amount

number = 1000

# The line below defines symbols needed for the function g, DO NOT CHANGE  IT
Yi,ci,yi,xi,u,Lami = symbols('Yi ci yi xi u Lami')

# Define the function g. This is essentiallly the Lagrangian of the field
# divided by the kinetic energy of the field. For Quintessence this would be 
#1-(ci/Yi), for Ghost Condensates the function is -1+(ci*Yi). To define your 
#own please refer to the readme file or to Tsujikawa and Ohashi.

g = 1-(ci/Yi)

# Enter the total number of fields used in each model

fields = 20

# Enter the ranges for each parameter used in your model. The ranges given in the
# UI are the ranges supplied here. The x parameters are analogous to the kinetic
# energy of the fields,the y parameters are analogous to the potential energy,
# the c parameters are the coefficients of the assumed exponential potentials, and
# lastly the u parameter defines the initial density of radiation in the universe.
# Note that for the Ghost Condensate model you must make the string GC = "True"
# In this case the values of x2 and y2 are defined by x1 and y1 to make sure that
# instabilities do not arise

x1min,x1max = .001 , 0.3
y1min,y1max = .001 , 0.3
x2min,x2max = 10**(-14.0) , 10**(-12.0)
y2min,y2max = 10**(-14.0) , 10**(-12.0)
c1min,c1max = .001, 1.0
c2min,c2max = .001, .002
umin,umax = 0.7,0.9

GC = "False"

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

loc = os.path.dirname(sys.argv[0])

ParFile(loc,name, str(g), (str(g2)), (str(g3)), (str(Aterm)), x1min, x1max,x2min, x2max, y1min, y1max,y2min, y2max,c1min, c1max,c2min, c2max, umin, umax, fields)

Tests = np.array([])

if GC != "True":
    for num in xrange(number):
        Q = Universe(uniform(x1min,x1max),uniform(x2min,x2max),uniform(y1min,y1max),uniform(y2min,y2max),uniform(c1min,c1max),uniform(c2min,c2max),uniform(umin,umax),fields,str(g),str(g2),str(g3),str(Aterm),e1,e2) 
        Tests = np.append(Tests, Q)
else:
    for num in xrange(number):
        x1 = uniform(x1min,x1max)
        x2 = uniform(x2min,x2max)
        Q = Universe(x1,x2,uniform(.5,1.0)*x1,uniform(.5,1.0)*x2,uniform(c1min,c1max),uniform(c2min,c2max),uniform(umin,umax),fields,str(g),str(g2),str(g3),str(Aterm),e1,e2)
        Tests = np.append(Tests, Q)

def Runner(N):
    N.Run(0)
    return N
    
Range = np.arange(0,number,1)
    
if __name__=='__main__':
    mp.freeze_support()
    pool = mp.Pool()
    Tests = pool.map(Runner,Tests)
    pool.close()
    pool.join()
    Data = mp.Process(target = ResFile, args = (loc, name, fields, Tests))
    Grow = mp.Process(target = GroFile, args = (loc, name, fields, Tests))
    Data.start()
    Data.join()
    Grow.start()
    Grow.join()