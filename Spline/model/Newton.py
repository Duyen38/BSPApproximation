import numpy as np
import math
import time
import matplotlib.pyplot as plt

def func(x: float):
    return float(0.5 * np.sin(0.5*math.pi * x) * np.cos(2*x - math.pi/2))

def divided_diff_coef(x, y):
    '''
    function to calculate the divided
    differences table
    '''
    n = len(y)
    coef = np.zeros([n, n])
    # the first column is y
    coef[:,0] = y
    
    for j in range(1,n):
        for i in range(n-j):
            coef[i][j] = (coef[i+1][j-1] - coef[i][j-1]) / (x[i+j]-x[i])
            
    return coef

def newton_poly(coef, x_data, x):
    '''
    evaluate the newton polynomial 
    at x
    '''
    n = len(x_data) - 1 
    p = coef[n]
    for k in range(1,n+1):
        p = coef[n-k] + (x -x_data[n-k])*p
    return p

def main():
    a, b = 0, 9
    n = 1000
    x = np.linspace(a, b, n)    
    y = []
    y = np.array([func(i) for i in x]) 
    x_vals = np.linspace(max(x), min(x), 1000)        # range X: [max, min]
    Fx = np.array([func(i) for i in x_vals])      # функция f(x)
    
    ''' lagrange_polynomial '''
    start = time.time()
    coef = divided_diff_coef(x, y)[0, :]
    y_newton = newton_poly(coef, x, x_vals)
    end = time.time()
    print("Time method newton:", end - start)

    plt.figure('Интерполяция')
    plt.title('f(x) = 0.5 * sin(0.5x * pi) * cos(2x - pi/2), n = ' + str(n) 
                +'\nCтепень интерполяционного полинома: '+ str(n - 1), fontsize=10)          #Degree of the polynomial approximation: 
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)

    # Plot the function with points
    plt.plot(x_vals, Fx, alpha=0.7, linewidth=2, label='Функция f(x)') #Funcion
    plt.plot(x, y, 'bo')

    # Plot the interpolation polynomial Polynomial interpolation
    plt.plot(x_vals, y_newton, 'k', alpha=0.6, label='Интерполяционный многочлeн Ньютона')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()