'''This calculation is designed to calculate psi and delta for ellipsometry. 
The formulas are taken from the book 'A User's Guide to Ellipsometry by Harland G.Tompkins'. I have three interference here. 
The top one is air (with refractive index, n0=1), the second layer is my material (with refractive index=n1) I want to measure the thickness and refractive index of, 
and the final layer is the substrate, which is silicon(with refractive index=n3). I have assumed the absorption coefficient to be 0 for both material and silicon since it was close to 0.'''
from __future__ import division
import cmath
import math
from math import sin, cos, tan

'''The first calculation is to calculate refractive index using Cauchy's equation'''
#x is the wavelength
#the formula used here is Cauchy's equation that is the emperical relationship between wavelength and refractive index.
# N0 and N1 are Cauchy's constants in the formula
def ref_index_ormosil():
	N0= 1.550 
	N1= 70.6
	result = []
	for x in range(400,801): #unit is nm
		n1= N0 + N1*100/x**2 
		print x, n1
		result.append((x, n1))
	return result
result1 = ref_index_ormosil()
import matplotlib.pyplot as plt
x = [item[0] for item in result1]
y = [item[1] for item in result1]
plt.plot(x,y)
plt.title('Refractive index vs Wavelength')
plt.xlabel('Wavelength/nm')
plt.ylabel('Refractive index')
plt.show()

'''The second calculation is finding the refractive index of silicon wafer since it also depends on the value of wavelength.'''
def ref_index_si():
	N0= 3.872
	N1= 100
	result = []
	for x in range(400,801):  #unit is nm
		n3= N0 + N1*100/x**2
		print x, n3
		result.append((x, n3))
	return result
result2= ref_index_si()


'''The third calculation is finding reflection_12 in p-polarisation'''
import numpy as np
import scipy as sp 
n1= 1.0
def reflectance_p_12(result1):
	'''The angle of incidence, Q is 70 degree, which is 1.2217 radians'''
	Q1 = np.cos(1.2217)
	result = []
	for item in result1:
		n2= item[1]
		position = item[0]
		Q2 = np.cos(sp.arcsin(np.real_if_close(n1/n2*np.sin(1.2217))))
		r_p_12 = ((n2*Q1)-(n1*Q2))/((n2*Q1)+(n1*Q2))
		print (Q1, Q2, r_p_12)
		result.append((position, r_p_12)) #x is the position
	return result
result3= reflectance_p_12(result1)


'''The fourth calculation is finding reflection_12 in s-polarisation'''
import math 
n1= 1.0
def reflectance_s_12(result1):
	Q1 = np.cos(1.2217)
	result = []
	for item in result1:
		n2= item[1]
		position = item[0]
		'''by snell's law: the relationship between n1 and n2 is: sp.arcsin((n1/n2*np.sin(Q1)))'''
		Q2 = np.cos(sp.arcsin(np.real_if_close(n1/n2*np.sin(1.2217))))
		r_s_12 = ((n1*Q1)-(n2*Q2))/((n1*Q1)+(n2*Q2))
		print r_s_12
		result.append((position, r_s_12))
	return result
result4= reflectance_s_12(result2)

'''The fifth calculation is finding reflection_23 in p-polarisation'''
n1= 1.0
import math 
def reflectance_p_23(result1):
	Q1 = np.cos(1.2217)
	result = []
	for position, item in enumerate(result1):
		x = item[0]
		n2= item[1]
		n3= result2[position][1]
		Q2_1 = sp.arcsin(np.real_if_close(n1/n2*np.sin(1.2217)))
		Q2 = np.cos(Q2_1)
		Q3 = np.cos(sp.arcsin(np.real_if_close(n2/n3*np.sin(Q2_1))))
		r_p_23 = ((n3*Q2)-(n2*Q3))/((n3*Q2)+(n2*Q3))
		print r_p_23
		result.append((x, r_p_23))
	return result
result5= reflectance_p_23(result1)

'''The sixth calculation is finding reflection_23 in s-polarisation'''
n1= 1.0
import math 
def reflectance_s_23(result1):
	Q1 = np.cos(1.2217)
	result = []
	for position, item in enumerate(result1):
		x = item[0]
		n2= item[1]
		n3= result2[position][1]
		Q2_1 = sp.arcsin(np.real_if_close(n1/n2*np.sin(1.2217)))
		Q2 = np.cos(Q2_1)
		Q3 = np.cos(sp.arcsin(np.real_if_close(n2/n3*np.sin(Q2_1))))
		r_s_23 = ((n2*Q2)-(n3*Q3))/((n2*Q2)+(n3*Q3))
		print r_s_23
		result.append((x, r_s_23))
	return result
result6 = reflectance_s_23(result1)

'''The seventh calculation is finding film phase thickness, B'''
n1= 1.0
d= 440.32 #put the thickness of your material, unit is nm
import math 
def film_thickness(d):
	result= []
	Q1 = math.cos(1.2217)
	for item in result1:
		n2= item[1]
		Q2 = np.cos(sp.arcsin(np.real_if_close(n1/n2*np.sin(1.2217))))
		x= item[0]
		B= 2*(math.pi)*(d/x)*n2*Q2
		print x
		print B
		result.append((x,B))
	return result
result7= film_thickness(d)

'''The eigth calculation is calculating total reflection coefficient in p-polarisation'''
def reflectance_P(result7):
	result = []
	for position, item in enumerate(result7):
	    x = item[0]
	    B = item[1]
	    A = np.exp(-2j * B) 
	    r_p_12 = result3[position][1] 
	    r_p_23 = result5[position][1]
	    print A, r_p_12
	    R_P = ((r_p_12)+((r_p_23)*A))/(1+((r_p_12)*(r_p_23)*A))
	    print position, R_P
	    result.append((x, R_P))
	return result
result8 = reflectance_P(result7)

'''The ninth calculation is calculating total reflection coefficient in s-polarisation'''
def reflectance_S(result7):
	result = []
	for position, item in enumerate(result7): 
		x = item[0]
		B = item[1]
		A = cmath.exp(-2j * B)
		r_s_12 = result4[position][1]
		r_s_23 = result6[position][1]
		print A, r_s_12
		R_S = ((r_s_12)+((r_s_23)*A))/(1+((r_s_12)*(r_s_23)*A))
		print position, R_S
		result.append((x, R_S))
	return result
result9 = reflectance_S(result7)

'''The tenth calculation is calculating rho which is the ratio of total reflection in p and s polarisation'''
def rho(result8, result9):
	result = []
	for position, item in enumerate(result8):
		x = item[0]       
		R_P = item[1]
		R_S = result9[position][1]
		rho = (R_P)/(R_S)
		result.append((x,rho))
	return result
result10= rho(result8, result9)

'''The eleventh calculation is calculating psi, which will be in radians'''
def get_psi_radians():
    result = []
    for item in result10:
        x = item[0]
        b = np.abs(item[1])
        psi = np.arctan(b)
        result.append((x, psi))
    return result
result11 = get_psi_radians()

'''The twelfth calculation is calculating psi in degree'''
def get_psi_degree(result11):
    result = []
    for position, item in enumerate(result11):
        x = item[0]
        radians= item[1]
        degree = radians*180/math.pi
        result.append((x, degree))
    return result
result12 = get_psi_degree(result11)

'''Now plotting the graph of psi vs wavelength'''
import matplotlib.pyplot as plt
x = [item[0] for item in result12]
y = [item[1] for item in result12]
plt.plot(x,y)
plt.title('Psi vs Wavelength')
plt.xlabel('Wavelength/nm')
plt.ylabel('Psi/Degree')
plt.show()

'''The thiteenth calculation is calculating delta in radians'''
def get_delta(result11): 
    result = []
    for position, item in enumerate(result11):
        x = item[0]
        psi = item[1]
        rho = result10[position][1]
        C= np.tan(psi)
        delta_0 = (np.log(-rho/C))
        delta_1 = np.imag(delta_0)
        result.append((x, delta_1))
    return result
result13 = get_delta(result11)

'''The fourteenth and final calculation is calculating delta in degree'''
def get_degree(result13):
    result = []
    for position, item in enumerate(result13):
        x = item[0]
        radians= item[1]
        degree = radians*180/math.pi
        result.append((x, degree))
    return result
result14= get_degree(result13)

'''Now plotting the graph of delta vs wavelength'''
import matplotlib.pyplot as plt
x = [item[0] for item in result14]
y = [item[1] for item in result14]
plt.plot(x,y)
plt.title('Delta vs Wavelength')
plt.xlabel('Wavelength/nm')
plt.ylabel('Delta/Degree')
plt.show()
