import sys
import math
import numpy as np
import scipy.linalg as la
import random
import time
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from scipy.interpolate import make_lsq_spline

dict_Bspline = {}       # complex dictionary keys to save values in Bspline
k = 3

def Bsp(x, i, k, t, h):
    if(x < t[i] or x >= t[i+4]): return 0
    elif(t[i] <= x and x < t[i+1]):    
        return (pow(x-t[i],3)/(6*pow(h,3)))
    elif(t[i+1] <= x and x < t[i+2]):
        return (pow(h,3) + 3*pow(h,2)*(x-t[i+1]) + 3*h*pow((x-t[i+1]),2) - 3*pow((x-t[i+1]),3)) / (6*pow(h,3))
    elif(t[i+2] <= x and x < t[i+3]):
        return (pow(h,3) + 3*pow(h,2)*(t[i+3]-x) + 3*h*pow((t[i+3]-x),2) - 3*pow((t[i+3]-x),3)) / (6*pow(h,3))
    else:
        return pow(t[i+4]-x, 3) / (6*pow(h,3))

def B(x, i, k, t):
    '''
    x: data
    i: spline function index
    k: degree of spline
    t: knot sequence
    '''
    if (x,i,k) in dict_Bspline:
        return dict_Bspline.get((x,i,k))
    else:
        if k == 0:
            return 1.0 if t[i] <= x < t[i+1] else 0.0
        c1, c2 = 0.0, 0.0
        if t[i+k] > t[i]:
            c1 = (x - t[i])/(t[i+k] - t[i]) * B(x, i, k-1, t)
        if t[i+k+1] > t[i+1]:
            c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * B(x, i+1, k-1, t)

        #add to dictionary   
        dict_Bspline[(x,i,k)] = c1+c2 
        return c1 + c2

def build_knot_vector_tmp(a, b, k, h_R): 
    x_R = []
    i = 0
    while (a+i*h_R) < b:
        x_R.append(a+i*h_R)
        i = i + 1
    knots = []
    for i in range(k,0, -1): 
        knots.append(a-i*h_R)
    for j in range(len(x_R)):
        knots.append(x_R[j])
    for i in range(k+1):
        knots.append(x_R[-1]+(i+1)*h_R) 
    return knots  

def build_matrix_bspline(x, k, knot, n, m, hR):
    A = np.zeros((n, m))

    for i in range(n):
        j_start = i*(m-2)/(n-1) - 2
        j_end = i*(m-2)/(n-1) + 3
        
        j_start = math.floor(j_start)
        j_end = math.ceil(j_end)

        if(j_start < 0): j_start = 0
        if(j_end > m): j_end = m
        for j in range(j_start,j_end):
            A[i][j] = Bsp(x[i], j, k, knot, hR)
    return A

def get_coeff_appromation(x, y, k, knot, hR):
    '''
    Solve systems of linear equations by method Cholesky for calculate coefficients
    (the least-squares solution to a linear matrix equation)

    '''
    n = len(x)
    m = len(knot)-k-1
    A = build_matrix_bspline(x, k, knot, n, m, hR)
    coeff = []
    if m==n:
        coeff = np.linalg.inv(A) @ y
    elif m < n:
        A_hat = A.T @ A
        y_hat = A.T @ y                  
        L = la.cholesky(A_hat, lower=True)  # lower triangular Cholesky factorization A = L @ L.T
        beta = np.linalg.inv(L) @ y_hat       # L.T @ beta = y_hat
        coeff = np.linalg.inv(L.T) @ beta   # result A @ coeff = Y

    else:
        coeff = A.T @  np.linalg.inv(A @ A.T) @ y
    return coeff


def get_coeff_interpolation(x, y, k, knot, hR):
    n = len(x)
    A = build_matrix_bspline(x, k, knot, n, n, hR)
    # A[-1][-1] = 1
    # coeff = np.linalg.inv(A_hat) @ y_hat
    return np.linalg.inv(A.T @ A) @ (A.T @ y)


def bspline_element(x, k, c, knot, hR):
    '''
     k: degree
     c: coefficients
     knot: knot
    '''
    n = len(knot) - k - 1
    assert (n >= k+1) and (len(c) >= n)
    return sum(c[i] * Bsp(x, i, k, knot, hR) for i in range(n))

def func(x):
    return float(0.5 * np.sin(0.5*math.pi * x) * np.cos(2*x - math.pi/2)) #0.5 * np.sin(0.5*math.pi * x) * np.cos(2*x - math.pi/2))#x*x - 4*x - 3

def plotGraph(n, h, R, x, y, x_vals, bspline):
    plt.figure('B-spline')
    plt.title('В-сплайн третьего порядка (k = 3) \n n = ' + str(n) + '\n за шагом h = ' + str(round(h, 2))
                +'; h_R = ' + str(round(h*R, 2)))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.plot(x, y, 'ro', label='Points')

    # plt.plot(x_vals, np.array([func(i) for i in x_vals]), alpha=0.7, label='Функция f(x)')

    # plt.plot(x_vals, bspline, label='B-spline interpolation')
    plt.plot(x_vals, bspline, lw='2', label='B-spline аpproximation')
    
    plt.grid(True)
    plt.legend()
    plt.show()

def main():
    a, b = 0, 20
    n = 100
    h = round((b - a) / (n-1), 2)
    R = 2
    h_R = h*R 
    
    x = np.linspace(a, b, n)    
    y = []
    # for i in range(len(x)):
    #     y.append(random.randint(-20, 20)) 
    y = np.array([func(i) for i in x])      # 
    # y = np.exp(-x**2) + 0.1 * np.random.randn(n) 
    x_vals = np.linspace(min(x), max(x), 1000)  

    '''B-spline interpolation''' 
    # knots = build_knot_vector_tmp( a, b, k, h) 

    # start = time.time()
    # # print("knots: " + str(len(knots)) + ' ' + str(knots))
    # coeff = get_coeff_appromation(x, y, k, knots)  #get_coeff_interpolation(x, y, k, knots)
    # bspline = [bspline_element(i, k, coeff , knots, hR) for i in x_vals]
    # end = time.time()
    # print("interpolate's time:", end - start)

    '''B-spline аpproximation'''  
    knots1 = build_knot_vector_tmp( a, b, k, h_R)

    start = time.time()
    coeff1 =  get_coeff_appromation(x, y, k, knots1,h_R)
    bspline_1 = [bspline_element(i, k, coeff1 , knots1, h_R) for i in x_vals]

    end = time.time()
    print("approach's time:", end - start)

    plotGraph(n, h, R, x, y, x_vals, bspline_1)

if __name__ == '__main__':
    main()   
