# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 20:22:41 2019
CLSE 305 Project 2 due 10/25/19 12:30pm
@author: Nathan Lisech
Mixing Rules PREOS (methane, ethane, propane)
"""
from FindZ import FindDen
import numpy as np
from openpyxl import load_workbook

filepath=r"C:\Users\Nathan Lisech\.spyder-py3\Mixing Rules PREOS ExcelWB.xlsx"
wb = load_workbook(filepath)
worksheet = wb.active


y1 = .7 #methane CH4
y2 = .2 #Ethane C2H6
y3 = .1 # Propane C3H8

# pure component a 
a_methane = 1.7846 * 10**(-6)
a_ethane = (5.4 * 10**-6)
a_propane = (1.0125 * 10**-5)

b_methane = (2.68 * 10**-5)
b_ethane = (4.054 *10**-5)
b_propane = (5.63 * 10**-5)

#interaction parameters (Table 9.4-1 in book)
k12 = .022
k13 = .033
k23 = .001

a12 = ((a_methane * a_ethane)**(1/2) * (1-k12)) # = a21
print(a12, "a12")
a13 = ((a_methane * a_propane)**(1/2) * (1-k13)) # = a31 
print(a13, "a13")
a23 = ((a_ethane * a_propane)**(1/2) * (1-k23)) # = a32
print(a23, "a23")

amix = (y1*y1*a_methane) + (y1*y2*a12) + (y1*y3*a13) + (y2*y1*a12) + (y2*y2*a_ethane) + (y2*y3*a23) + (y3*y1*a13) + (y3*y2*a23) + (y3*y3*a_propane)
print(amix, "amix")

bmix = (y1*b_methane) + (y2*b_ethane) + (y3*b_propane)
print(bmix, "bmix")

# Parameters for pure component
print("Initial Conditions : Tstart = 100 K , Pstart = .001 bar")
    #R = float(x=83.14) #R-value = cm3 bar K−1 mol−1 
R = float(x=.0000831447)#R-value = m3 bar K−1 mol−1 
#worksheet['B5'] = R
NNN = 0

for T in range(0,1):
    Tstart = int(x=100) #Starting temperature in Celsius
    start = str(Tstart)
    Tstop = int(Tstart)  # Ending Temperature in Celsius
    Trange = str(Tstop) #str version of Tstop

    for p in range(0,4):
        P_ref = 1 #1 bar
        Pstart = float(x=5) #Starting Pressure in bar
        #worksheet['D2'] = Pstart
        Pstep = int(x=5) #Change in Pressure (bar)
        #worksheet['D4'] = Pstep
        Pstop = int(Pstart + (Pstep * p)) #Ending Pressure in bar
        #worksheet['D3'] = Pstop
        print(Pstop, "P in barr")
        
        A = amix #from Mixing Rules PREOS Program
        #print(A, "A")
        B = bmix
        #print(B, "B")
        AA = (A  * (Pstop) / ( (R * (Tstop + 273 )**2)))
        AAmix = AA
        BB = (B * (Pstop) / ( R * (Tstop + 273)))
        BBmix = BB
        print(AA, "AA")
        print(BB, "BB")
        RA1 = 1
        #print(RA1, "RA1")
        RA2 = (BB - 1)
        #print(RA2, "RA2")
        RA3 = (AA) - (2*(BB)) - (3*(BB)**2)
        #print(RA3, "RA3")
        RA4 = (BB)**3 + (BB)**2 - ((AA)*(BB))
        #print(RA4, "RA4")
        Z = FindDen(RA1,RA2,RA3,RA4,NNN)
        print(Z, " Calculated Zmix value")
        print() #whitespace
        Zmix = Z

#solved from main program and copied down
BBmethane5barr = .0043222 
BBmethane10barr = .0086435
BBmethane15barr = .0129649
BBmethane20barr = .0172862

BBethane5barr = .00653696
BBethane10barr = .0130726
BBethane15barr = .019608
BBethane20barr = .0261439

BBpropane5barr = .00908124
BBpropane10barr = .01860668
BBpropane15barr = .02724009
BBpropane20barr = .03631952
#solved BBmix by placing Amix and Bmix values into main program and solving to get AAmix and BBmix for each pressure
BBmix5barr = .005239420
BBmix10barr = .010478841
BBmix15barr = .0157182608
BBmix20barr = .02095873

BB1ratio5barr = (BBmethane5barr / BBmix5barr) # b1 / bmix @ 5 barr
BB1ratio10barr = (BBmethane10barr / BBmix10barr)
BB1ratio15barr = (BBmethane15barr / BBmix15barr)
BB1ratio20barr = (BBmethane20barr / BBmix20barr)

BB2ratio5barr = (BBethane5barr / BBmix5barr) # b2 / bmix @ 5 barr
BB2ratio10barr = (BBethane10barr / BBmix10barr)
BB2ratio15barr = (BBethane15barr / BBmix15barr)
BB2ratio20barr = (BBethane20barr / BBmix20barr)

BB3ratio5barr = (BBpropane5barr / BBmix5barr) # b3 / bmix @ 5 barr
BB3ratio10barr = (BBpropane10barr / BBmix10barr)
BB3ratio15barr = (BBpropane15barr / BBmix15barr)
BB3ratio20barr = (BBpropane20barr / BBmix20barr)

residual = (Zmix - 1)
ratio1 =( AAmix / (2*np.sqrt(2))*BBmix )
summation1 = 2 * (y1*a_methane + y2*a12 + y3*a13) #comp 1 
summation2 = 2 * (y1*a12 + y2*a_ethane + y3*a23) #comp 2 
summation3 = 2 * (y1*a13 + y2 *a23 + y3 *a_propane) #comp3

ln = np.log(Z - BBmix)
ln2 = np.log(Z + (1 + np.sqrt(2))*BBmix) / (Zmix + (1 - np.sqrt(2))*BBmix)
# vapor fugacity calculation
print( "Fugacity due to component 1, Methane")
fug1_partial5barr = (np.exp((BB1ratio5barr)*(residual) - ln - ratio1*( (summation1 / AAmix) - BB1ratio5barr  )*ln2)* (y1* 5)) #barr
fug1_partial10barr = (np.exp((BB1ratio10barr)*(residual) - ln - ratio1*( (summation1 / AAmix) - BB1ratio10barr  )*ln2) * (y1* 10)) #barr
fug1_partial15barr = (np.exp((BB1ratio15barr)*(residual) - ln - ratio1*( (summation1 / AAmix) - BB1ratio15barr  )*ln2)* (y1* 15)) #barr
fug1_partial20barr = (np.exp((BB1ratio20barr)*(residual) - ln - ratio1*( (summation1 / AAmix) - BB1ratio20barr  )*ln2)* (y1* 20)) #barr
print(fug1_partial5barr, "fugacity1 @ 5barr")
print(fug1_partial10barr, "fugacity1 @ 10barr")
print(fug1_partial15barr, "fugacity1 @ 15barr")
print(fug1_partial20barr, "fugacity1 @ 20barr")
print()

print( "Fugacity due to component 2, Ethane")
fug2_partial5barr = (np.exp((BB2ratio5barr)*(residual) - ln - ratio1*( (summation2 / AAmix) - BB2ratio5barr  )*ln2)* (y2* 5)) #barr
fug2_partial10barr = (np.exp((BB2ratio10barr)*(residual) - ln - ratio1*( (summation2 / AAmix) - BB2ratio10barr  )*ln2) * (y2* 10)) #barr
fug2_partial15barr = (np.exp((BB2ratio15barr)*(residual) - ln - ratio1*( (summation2 / AAmix) - BB2ratio15barr  )*ln2)* (y2* 15)) #barr
fug2_partial20barr = (np.exp((BB2ratio20barr)*(residual) - ln - ratio1*( (summation1 / AAmix) - BB2ratio20barr  )*ln2)* (y2* 20)) #barr
print(fug2_partial5barr, "fugacity2 @ 5barr")
print(fug2_partial10barr, "fugacity2 @ 10barr")
print(fug2_partial15barr, "fugacity2 @ 15barr")
print(fug2_partial20barr, "fugacity2 @ 20barr")
print()

print( "Fugacity due to component 3, Propane")
fug3_partial5barr = (np.exp((BB3ratio5barr)*(residual) - ln - ratio1*( (summation3 / AAmix) - BB3ratio5barr  )*ln2)* (y3* 5)) #barr
fug3_partial10barr = (np.exp((BB3ratio10barr)*(residual) - ln - ratio1*( (summation3 / AAmix) - BB3ratio10barr  )*ln2) * (y3* 10)) #barr
fug3_partial15barr = (np.exp((BB3ratio15barr)*(residual) - ln - ratio1*( (summation3 / AAmix) - BB3ratio15barr  )*ln2)* (y3* 15)) #barr
fug3_partial20barr = (np.exp((BB3ratio20barr)*(residual) - ln - ratio1*( (summation1 / AAmix) - BB3ratio20barr  )*ln2)* (y3* 20)) #barr
print(fug3_partial5barr, "fugacity3 @ 5barr")
print(fug3_partial10barr, "fugacity3 @ 10barr")
print(fug3_partial15barr, "fugacity3 @ 15barr")
print(fug3_partial20barr, "fugacity3 @ 20barr")
print()

#check function
CHECK1 = fug1_partial5barr + fug2_partial5barr + fug3_partial5barr
print(CHECK1, "CHECK: Total Fugacity")
CHECK2 = fug1_partial10barr + fug2_partial10barr + fug3_partial10barr
print(CHECK2, "CHECK: Total Fugacity")
CHECK3 = fug1_partial15barr + fug2_partial15barr + fug3_partial15barr
print(CHECK3, "CHECK: Total Fugacity")
CHECK4 = fug1_partial20barr + fug2_partial20barr + fug3_partial20barr
print(CHECK4, "CHECK: Total Fugacity")