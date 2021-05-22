import sys
import math
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
from dataclasses import dataclass
from itertools import islice
from typing import Callable, List
import warnings
import time
import random

X, Y = [], []
def function(x: float):
    return float(0.5 * np.sin(0.5*math.pi * x) * np.cos(2*x - math.pi/2)) # / 6 + 2* math.pi/0.2)) 
''' 0.5 * np.sin(0.5*math.pi * x) * np.cos(2*x - math.pi/2)
    x * np.sin(x) / (1 + x*x)
    x / (1 + x*x)
    (np.exp(-(x**2)/2))/(2*np.pi)   '''

@dataclass
class Spline:
    a: float = 0.0
    b: float = 0.0
    c: float = 0.0
    d: float = 0.0
    x: float = 0.0
    y: float = 0.0

class CubicSplineInterpolator():
    # Initialize class instance with values
    def __init__(self, left_boundary: float, right_boundary: float,
            epsilon: float, intervals: int, function: Callable ):

        self.left_boundary = left_boundary
        self.right_boundary = right_boundary
        self.epsilon = epsilon
        self.intervals = intervals
        self.function = function

        self.splines = self.buildSpline()
        self.x_vals, self.S = self.interpolate()
    
    # Решить трехдиагональную систему методом прогонки, чтобы опреелить коэфф. c
    def solve_equations_system(self, splines: List[Spline], h: float, Y: List[float]):
        # прогоночные коэфф.
        alpha = np.zeros(self.intervals - 1)
        beta = np.zeros(self.intervals - 1)

        # Прямая подстановка - изменение коэффициентов
        for i in range(1, self.intervals - 1):
            alpha[i] = -1 / (4 + alpha[i-1])
            beta[i] = 1 / (4 + alpha[i-1]) * (6 / h**2 * (Y[i + 1] - 2 * Y[i] + Y[i - 1]) - beta[i - 1])

        # Обратная подстановка - получение решения
        for i in range(self.intervals-2, 0, -1):
            splines[i].c = alpha[i] * splines[i+1].c + beta[i]

    def buildSpline(self):
        h = (self.right_boundary - self.left_boundary) / (self.intervals - 1)

        splines = []
        a_coeff = []    # значения коэфф. ai = Yi

        # for i in x_data:
        #     splines.append(Spline(a=function(i), b=0.0, c=0.0, d=0.0, x=i, y=function(i)))

        for i in range(self.intervals): # (len(x_data)): 
        # for i in x_data:
            value = self.left_boundary + i * h
            result = self.function(value)

            X.append(value)
            Y.append(result)

            a_coeff.append(result)
            splines.append(Spline(a=result, b=0.0, c=0.0, d=0.0, x=value, y=function(value)))

        splines[0].c = 0.0

        self.solve_equations_system(splines, h, a_coeff)

        # Определить значения коэфф. b и d
        for i in range(self.intervals - 1, 0, -1):
            splines[i].d = (splines[i].c - splines[i - 1].c) / h
            splines[i].b = (h / 2 * splines[i].c - h**2 / 6 * splines[i].d + (a_coeff[i] - a_coeff[i - 1]) / h)

        return splines

    def interpolate(self):
        x_values, S = [], []

        for i in range(1, self.intervals):
            x_k = self.splines[i - 1].x
            while x_k <= self.splines[i].x:
                Si = (self.splines[i].a
                    + self.splines[i].b * (x_k - self.splines[i].x)
                    + self.splines[i].c / 2 * (x_k - self.splines[i].x)**2
                    + self.splines[i].d / 6 * (x_k - self.splines[i].x)**3)

                x_values.append(x_k)
                S.append(Si)
                x_k += self.epsilon

        return x_values, S

    def print_coefficients(self):
        print('  k  |    x    |    y    |    a    |    b    |    c    |    d    ')
        print('-----------------------------------------------------------------')

        # splines_slice = islice(self.splines, 15)    # выводить 15 первые значения 

        for k, spline in enumerate(self.splines, start=1):
            print(f' {k:3} | {spline.x:7.3f} | {spline.y:7.3f} | {spline.a:7.3f} | {spline.b:7.3f} | {spline.c:7.3f} | {spline.d:7.3f}')

        
def main():
    a, b = 0, 9
    n = 5000
    x = np.linspace(a, b, n)    
    y = []
    y = np.array([function(i) for i in x]) 

    x_vals = np.linspace(max(x), min(x), 1000)     
    Fx = np.array([function(i) for i in x_vals])      # функция f(x)

    start = time.time()
    # cubic_spline = CubicSplineInterpolator(range_start, range_end, epsilon, intervals, function)
    cubic_spline = CubicSplineInterpolator(a, b, (b-a)/1000, n, function)

    end = time.time()
    print("Time method cubic spline:", end - start)
  
    # print('\nРасчет коэффицентов кубического сплайна:\n')
    # cubic_spline.print_coefficients()

# 0.5 * sin(0.5x * pi) * cos(2x - pi/2)
# xsin(x) / (1 + x*x)
# x / (1 + x*x)     exp(-(x**2)/2) / 2pi
    plt.figure('Интерполяция')
    plt.title('f(x) = 0.5 * sin(0.5x * pi) * cos(2x - pi/2), n = ' + str(n))          #Degree of the polynomial approximation: 
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)

    # Plot the function with points
    plt.plot(x_vals, Fx, alpha=0.7, linewidth=2, label='Функция f(x)') #Funcion
    plt.plot(x, y, 'bo')

    # Plot the cubic spline interpolation кубический сплайн Cubic spline
    plt.plot(cubic_spline.x_vals, cubic_spline.S, 'r', label='Кубический сплайн') #b30086

   
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
