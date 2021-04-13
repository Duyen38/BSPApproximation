import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style
from matplotlib import pyplot as plt
from tkinter import messagebox

import tkinter as tk
from tkinter import ttk     # as css in tkinter ))

import urllib
import json
import pandas as pd     # load data from csv files
import numpy as np
import time
import math
import B_Spline
import CubicSpline
import Lagrange

LARGE_FONT = ("Verdana", 12)
NORM_FONT = ("Verdana", 10)
SMALL_FONT = ("Verdana", 8)
style.use("ggplot")

fig = Figure()      # figsize=(5,5), dpi=100
graph = fig.add_subplot(111)

fig1 = Figure()
graph1 = fig1.add_subplot(111)

fig2 = Figure()
graph2 = fig2.add_subplot(111)

def func(x: float):
    return float(0.5 * np.sin(0.5*math.pi * x) * np.cos(2*x - math.pi/2))

# paneCount = 1

exchange = "B-spline"
DatCounter = 9000
programName = "bsp"

def changeExchange(toWhat, pn):
    global exchange
    global DatCounter
    global programName

    exchange = toWhat
    programName = pn
    DatCounter = 9000

def popupmsg(msg):
    popup = tk.Tk()

    popup.wm_title("!")
    label = ttk.Label(popup, text=msg, font=NORM_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()
    popup.mainloop()

def animate(i):
    # pullData = open("sampleData.txt", "r").read()
    # dataList = pullData.split('\n')
    # xList = []
    # yList = []
    # for eachLine in dataList:
    #     if len(eachLine) > 1:
    #         x, y = eachLine.split(',')
    #         xList.append(float(x))
    #         yList.append(float(y))
    a, b = 0, 50
    n = 200
    x = np.linspace(a, b, n)    
    y = []
    y = np.array([func(i) for i in x]) 
    # y = np.exp(-x**2) + 0.1 * np.random.randn(n) 
    
    x_vals = np.linspace(min(x), max(x), 1000) 

    h = round((b - a) / (n-1), 2)
    R = 5
    h_R = h*R 

    knots1 = B_Spline.build_knot_vector_tmp( a, b, 3, h_R)

    time_start = time.time()
    coeff1 =  B_Spline.get_coeff_appromation(x, y, 3, knots1,h_R)
    bspline_draw = [B_Spline.bspline_element(i, 3, coeff1 , knots1, h_R) for i in x_vals]
    time_end = time.time()
    print("approach's time:", time_end - time_start)

    graph.clear()

    graph.set_title('В-сплайн третьего порядка \n n = ' + str(n) + '\n за шагом h = ' + str(round(h, 2))
                +'; h_R = ' + str(round(h*R, 2)) + "\nApproach's time: " +  str(time_end - time_start))
    graph.plot(x, y, 'ro', label='Points')
    graph.plot(x_vals, bspline_draw, lw='2', label='B-spline аpproximation')
    graph.legend(bbox_to_anchor=(0, 1.02, 1, .102), loc=3, ncol=2, borderaxespad=0)

    # Plot the cubic spline interpolation кубический сплайн 
    time_start = time.time()
    cubic_spline =  CubicSpline.CubicSplineInterpolator(a, b, (b-a)/1000, n, func) 
    # y_lagrange = np.array([Lagrange.lagrange_polynomial(x, y, x_val) for x_val in x_vals])               
    time_end = time.time()
    print("Time method cubic spline:", time_end - time_start)

    graph1.clear()
    graph1.set_title('Cubic Spline Interpolation\n n = ' + str(n)  + "\nTime: " +  str(time_end - time_start))
    graph1.plot(x, y, 'ro', label='Points')
    graph1.plot(cubic_spline.x_vals, cubic_spline.S, 'r', label='Cubic Spline')
    # graph1.plot(x_vals, y_lagrange, 'k', label='Интерполяционный многочлeн Лагрaнжа')
    graph1.legend(bbox_to_anchor=(0, 1.02, 1, .102), loc=3, ncol=2, borderaxespad=0)



class BSPApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        menubar = tk.Menu(container)
        filemenu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="Open file", command=lambda: popupmsg("Not suported just yet!")) 
        filemenu.add_command(label="Save", command=lambda: popupmsg("Not suported just yet!"))       #messagebox.showwarning("Warning","Not suported just yet!") 
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=quit)

        exchangeChoice = tk.Menu(menubar, tearoff=0)
        exchangeChoice.add_command(label="B-Spline", command=lambda: changeExchange("B-spline","bspline"))
        exchangeChoice.add_command(label="Cubic Spline Interpolation", command=lambda: changeExchange("Cubic-spline","cubicspline"))
        exchangeChoice.add_command(label="Polynomial Interpolation", command=lambda: changeExchange("Polynomial","polyinter"))

        menubar.add_cascade(label="Exchange", menu=exchangeChoice)

        # parameter = tk.Menu(menubar, tearoff=1)
        # parameter.add_command(label="Setting x", command=lamb)

        tk.Tk.config(self, menu=menubar)
      
        self.frame = {}

        frame = SettingFrame(container, self)
        self.frame[SettingFrame] = frame
        frame.grid(row=0, column=0, sticky="nsew")

        for Fr in {MainPage, CubicSplinePage, BSplinePage, PolynomialInterPage}:
            frame = Fr(container, self)
            self.frame[Fr] = frame
            frame.grid(row=0, column=1, sticky="nsew")  # north south east west    

        self.show_frame(MainPage)

    def show_frame(self, control):
        frame = self.frame[control]
        frame.tkraise()     # overlap (to stack) frames on top of each other

class MainPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent) 
        label = tk.Label(self, text="Main page", font=LARGE_FONT)
        label.grid(row=0, column=0, columnspan=4, padx=10, pady=10, sticky="snew")

        button1 = ttk.Button(self, text="Cubic Spline", command=lambda: controller.show_frame(CubicSplinePage))
        button1.grid(row=1, column=0, padx=10, pady=10)

        button2 = ttk.Button(self, text="Bspline Page", command=lambda: controller.show_frame(BSplinePage))
        button2.grid(row=1, column=1, padx=10, pady=10)

        button2 = ttk.Button(self, text="Polynomial Interpolation", command=lambda: controller.show_frame(PolynomialInterPage))
        button2.grid(row=1, column=2, padx=10, pady=10)

        button3 = ttk.Button(self, text="exit", command=quit)
        button3.grid(row=1, column=3, padx=10, pady=10)

class SettingFrame(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        # sw = self.parent.winfo_screenwidth()

        label = tk.Label(self, text="Setting parameters", font=LARGE_FONT)
        label.grid(row=0, column=0, columnspan=2, padx=10, pady=10, sticky="WE")

        # button1 = ttk.Button(self, text="Back Home", command=lambda: controller.show_frame(MainPage))
        # button1.grid(row=0, column=0, columnspan=2, padx=10, pady=10, ipady=2, sticky="WE")

        x_lbl = tk.Label(self, text="X: ", font=NORM_FONT)
        x_lbl.grid(row=1, column=0, padx=5, pady=5, sticky="WE")
        x_entry = tk.Entry(self, width="20", bd=2)
        x_entry.grid(row=1, column=1, padx=5, pady=5, sticky="WE")

        y_lbl = tk.Label(self, text="Y: ", font=NORM_FONT)
        y_lbl.grid(row=2, column=0, padx=5, pady=5, sticky="WE")
        y_entry = tk.Entry(self, width="20", bd=2)
        y_entry.grid(row=2, column=1, padx=5, pady=5, sticky="WE")

        def isNumber():
            try:
                int(x_entry.get())
                int(y_entry.get())
                print("par number is added")
            except ValueError:
                print("Validate an Entry! You need to entry a number!")

        btn_add = tk.Button(self, text="Add", command=isNumber)
        btn_add.grid(row=3, column=1, padx=5, pady=5, sticky="WE")

        R_lbl = tk.Label(self, text="R: ", font=NORM_FONT)
        R_lbl.grid(row=4, column=0, padx=5, pady=5, sticky="WE")
        R_scale = tk.Scale(self, from_=0, to=200, orient=tk.HORIZONTAL)
        R_scale.grid(row=4, column=1, padx=5, pady=5, sticky="WE")
        

class CubicSplinePage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Cubic Spline Interpolation Page", font=LARGE_FONT)
        label.pack(padx=10, pady=10)

        button1 = ttk.Button(self, text="Back Home", command=lambda: controller.show_frame(MainPage))
        button1.pack()
        
        canvas = FigureCanvasTkAgg(fig1, self)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class BSplinePage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Bspline Page", font=LARGE_FONT)
        label.pack(padx=10, pady=10)

        button1 = ttk.Button(self, text="Back Home", command=lambda: controller.show_frame(MainPage))
        button1.pack()

        canvas = FigureCanvasTkAgg(fig, self)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class PolynomialInterPage(tk.Frame):

     def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Polynomial Interpolation Page", font=LARGE_FONT)
        label.pack(padx=10, pady=10)

        button1 = ttk.Button(self, text="Back Home", command=lambda: controller.show_frame(MainPage))
        button1.pack()

        canvas = FigureCanvasTkAgg(fig2, self)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

app = BSPApp()
app.title('BSP_APP')
app.geometry("1280x720")
ani = animation.FuncAnimation(fig, animate, interval=5000)
app.mainloop()