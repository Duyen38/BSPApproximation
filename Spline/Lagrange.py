import sys
import math
import numpy as np
import time

import matplotlib.pyplot as plt

def func(x: float):
    return float(0.5 * np.sin(0.5*math.pi * x) * np.cos(2*x - math.pi/2))

def lagrange_polynomial(x, y, x_vals):
    sum = 0
    n = len(x)
    for i in range(n):
        p = 1    
        for j in range(n):
            if i != j:
                p = p * (x_vals - x[j])/(x[i] - x[j])
        sum += p * y[i]
    return sum

def main():
    a, b = 0, 50
    n = 10
    x = np.linspace(a, b, n)    
    y = []
    y = np.array([func(i) for i in x]) 
    x_vals = np.linspace(max(x), min(x), 1000)        # range X: [max, min]
    Fx = np.array([func(i) for i in x_vals])      # функция f(x)
    
    ''' lagrange_polynomial '''
    start = time.time()
    y_lagrange = np.array([lagrange_polynomial(x, y, x_val) for x_val in x_vals])
    end = time.time()
    print("Time method Lagrange:", end - start)

    plt.figure('Интерполяция')
    plt.title('f(x) = 0.5 * sin(0.5x * pi) * cos(2x - pi/2), n = ' + str(n) 
                +'\nCтепень интерполяционного полинома: '+ str(n - 1), fontsize=10)          #Degree of the polynomial approximation: 
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)

    # Plot the function with points
    plt.plot(x_vals, Fx, alpha=0.7, linewidth=4, label='Функция f(x)') #Funcion
    plt.plot(x, y, 'bo')

    # Plot the interpolation polynomial Polynomial interpolation
    plt.plot(x_vals, y_lagrange, 'k', alpha=0.6, label='Интерполяционный многочлeн Лагрaнжа')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()