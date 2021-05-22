import sys
import math
import numpy as np
import warnings
import time
import matplotlib.pyplot as plt

def function(x: float):
    return float(0.5 * np.sin(0.5*math.pi * x) * np.cos(2*x - math.pi/2)) # / 6 + 2* math.pi/0.2)) 

def lst_squares(x, y, m, x_vals):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', np.RankWarning)
        coeff_appr = np.polyfit(x, y, m)            # coefficients the polynomial of degree M
    poly_approaching = np.poly1d(coeff_appr)
    y_vals = poly_approaching(x_vals)   
    return y_vals
        # poly_approaching = np.poly1d(c_appr)

def main():
    a, b = 0, 50
    n = 100
    m = 20      #degree polynomial
    
    x = np.linspace(a, b, n)    
    y = []
    y = np.array([function(i) for i in x]) 

    x_vals = np.linspace(max(x), min(x), 1000)     
    Fx = np.array([function(i) for i in x_vals])      # функция f(x)

    start = time.time()
    # coeff = lst_squares(x, y, m)
    # poly_approaching = np.poly1d(coeff)
    y_pa = lst_squares(x, y, m, x_vals)
    end = time.time()
    print("Time approch_poly:", end - start)

    plt.figure('Interpolation')
    plt.title('f(x) = 0.5 * sin(0.5x * pi) * cos(2x - pi/2), n = ' + str(n) 
                +'\nCтепень интерполяционного полинома: '+ str(n - 1)       # str(len(c_poly) - 1) Degree of the polynomial interpolation:
                +'\nCтепень аппроксимирующего полинома: '+ str(m), fontsize=10)          #Degree of the polynomial approximation: 
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)

    # Plot the function with points
    plt.plot(x_vals, Fx, alpha=0.7, linewidth=4, label='Функция f(x)') #Funcion
    plt.plot(x, y, 'bo')

    # Plot the interpolation polynomial Polynomial interpolation
    # plt.plot(x_curve, y_lagrange, 'k', alpha=0.6, label='Интерполяционный многочлeн Лагрaнжа')
    # # plt.plot(x_curve, y_poly, 'k', alpha=0.6, label='Интерполяционный полином')#00b300
    # if max(y_poly) > 2 or min(y_poly) < (-2):
    #     plt.ylim(min(X)-2, max(Y)+2, 'ro')
        
    # Plot the graph of a polynomial function approaches    Polynomial approximation
    plt.plot(x_vals, y_pa,'#e68a00', linewidth=2, label='Аппроксимирующий полином')
   
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()