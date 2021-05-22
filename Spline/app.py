import sys
import math
import time
import numpy as np
import pandas as pd     # load data from csv files
import json
import urllib
from tkinter.filedialog import Open
from tkinter import filedialog
from tkinter import ttk     # as css in tkinter ))
import tkinter as tk
from tkinter import messagebox
from matplotlib.widgets import Cursor
import matplotlib.widgets as widgets
from matplotlib.widgets import CheckButtons
from matplotlib import pyplot as plt
from matplotlib import style
import matplotlib.animation as animation
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import matplotlib
matplotlib.use("TkAgg")

import model.Least_squares as Least_squares
import model.Lagrange as Lagrange
import model.CubicSpline as CubicSpline
import model.B_Spline as B_Spline
import SnaptoCursor

import view.CubicSplinePage as CubicSplinePage

LARGE_FONT = ("Verdana", 14)
NORM_FONT = ("Verdana", 11)
SMALL_FONT = ("Verdana", 10)
style.use("ggplot")

fig = Figure()      # figsize=(5,5), dpi=100
graph = fig.add_subplot(111)

def func(x: float):
    return float(0.5 * np.sin(0.5*math.pi * x) * np.cos(2*x - math.pi/2))

def getY(x_list):
    return np.array([func(i) for i in x_list])

# -----global variables---
a, b = 0, 40
n = 200
x = []
y = []
x_vals = []
x = np.linspace(a, b, n)
y = np.array([func(i) for i in x])
x_vals = np.linspace(min(x), max(x), 1000)
R = 5   # coefficient for Bspline
h = round((b - a) / (n-1), 2)
h_R = h*R
M = 5

bspline = []
cubic_spline = []
y_poly_approach = []

mode_app = 1
ani = None
namePage = None

def addPointToGraph(x_val, y_val):
    global x, y, b, n
    x = np.append(x, x_val)
    y = np.append(y, y_val)
    b = max(x)
    n = n+1
    build()

def changeCoeffBSP(getvalScaleR):
    global R, h_R
    R = float(getvalScaleR)
    if h*R < (a+b)/2:
        h_R = h*R
        build()
    else:
        messagebox.showerror(title="Error", message="Value R over range! \nNot allow value!")

def changeNumberOfPoint(getvalScaleN):
    global n
    if(int(getvalScaleN) > 10):
        n = int(getvalScaleN)
        build()
    else:
        messagebox.showerror(title="Error", message="Number point invald! \nNot allow value!")

def changeDegreeM(getvalScaleM):
    global M
    if(int(getvalScaleM) > 0):
        M = int(getvalScaleM)
        build()
    else:
        messagebox.showerror(title="Error", message="Invalid value! \nNot allow value!")

def getEntryText(entry):
    try:
        value = float(entry.get())
        print('entry: ', value)
        return value
    except ValueError:
        print("Validate an Entry! You need to entry a number!")
        messagebox.showerror(title="Error", message="Invalid value! \nYou need to entry a number!")

cursor = None
cid = None
time_app_cubic_spline = 0
time_app_b_spline = 0
time_app_polymial = 0

def animate(i):
    global cursor
    global cid
    global graph
    global mode_app
    global ani
    
    st = time.time()

    print('mode_app', mode_app)
    if mode_app == 1:
        print('MainPage')        
        # Main Graph
        graph.clear()
        graph.plot(x, y, 'ro', label='Points')
        l0 = graph.plot(x_vals, bspline, lw='2', color='#5464cc', label='B-spline аpproximation')
        l1= graph.plot(cubic_spline.x_vals, cubic_spline.S, lw='2', color='#419134', label='Cubic Spline')
        l2 = graph.plot(x_vals, y_poly_approach,  color='k', label='Polynomial Approximation')
        graph.legend(bbox_to_anchor=(0, 1.02, 1, .102), loc=3, fontsize = 'medium', ncol=2, borderaxespad=0)
    
    if mode_app == 2:
        print('B-spline Approximation Page')

        graph.clear()
        graph.set_title('Time: ' + str(time_app_b_spline) + '(s)', loc='right')
        graph.plot(x, y, 'ro', label='Points')
        graph.plot(x_vals, bspline, lw='2', color='#5464cc',label='B-spline аpproximation')
        graph.legend(bbox_to_anchor=(0, 1.02, 1, .102), loc=3, fontsize = 'medium', ncol=2, borderaxespad=0)

        '''Defining the cursor'''
        cursor = SnaptoCursor.SnaptoCursor(graph, x_vals, bspline, getY(x_vals))
        cid = fig.canvas.mpl_connect('motion_notify_event', cursor.mouse_move)
    
    if mode_app == 3:
        print('Cubic Spline Interpolation Page')

        graph.clear()
        graph.set_title('Time: ' + str(time_app_cubic_spline)+ '(s)', loc='right')
        graph.plot(x, y, 'ro', label='Points')
        graph.plot(cubic_spline.x_vals, cubic_spline.S, color='#419134', label='Cubic Spline')
        graph.legend(bbox_to_anchor=(0, 1.02, 1, .102),
                        loc=3, fontsize = 'medium', ncol=2, borderaxespad=0)

        '''Defining the cursor'''
        cursor = SnaptoCursor.SnaptoCursor(graph, cubic_spline.x_vals, cubic_spline.S, getY(cubic_spline.x_vals))
        cid = fig.canvas.mpl_connect('motion_notify_event', cursor.mouse_move)

    if mode_app == 4:
        print('Polynomial Approximation Page')

        graph.clear()
        graph.set_title('Time: ' + str(time_app_polymial)+ '(s)', loc='right')
        graph.plot(x, y, 'ro', label='Points')
        graph.plot(x_vals, y_poly_approach, color='k', label='Polynomial Approximation')
        graph.legend(bbox_to_anchor=(0, 1.02, 1, .102), loc=3, fontsize = 'medium', ncol=2, borderaxespad=0)

        '''Defining the cursor'''
        cursor = SnaptoCursor.SnaptoCursor(graph, x_vals, y_poly_approach, getY(x_vals))
        cid = fig.canvas.mpl_connect('motion_notify_event', cursor.mouse_move)
    if ani is not None:
        ani.pause()
    end = time.time()
    print('time', end-st)

def build():
    global x, y, x_vals
    global ani
    global time_app_cubic_spline, time_app_b_spline, time_app_polymial
    x = np.linspace(a, b, n)
    y = np.array([func(i) for i in x])
    x_vals = np.linspace(min(x), max(x), 1000) 

    global bspline, cubic_spline, y_poly_approach
    start = time.time()
    knots1 = B_Spline.build_knot_vector_tmp(a, b, 3, h_R)
    coeff1 = B_Spline.get_coeff_appromation(x, y, 3, knots1, h_R)
    bspline = [B_Spline.bspline_element(i, 3, coeff1, knots1, h_R) for i in x_vals]
    end = time.time()
    time_app_b_spline = round((end - start),3)
    
    start = time.time()
    cubic_spline = CubicSpline.CubicSplineInterpolator(a, b, (b-a)/1000, n, func)
    end = time.time()
    time_app_cubic_spline = round((end - start),3)

    start = time.time()
    y_poly_approach = Least_squares.lst_squares(x, y, M, x_vals)
    end = time.time()
    time_app_polymial = round((end - start),3)
    
    if ani is not None:
        ani.resume()

class BSPApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        container.grid_rowconfigure(1, weight=5)
        container.grid_columnconfigure(1, weight=1)

        menubar = tk.Menu(container)
        filemenu = tk.Menu(menubar, tearoff=0)

        self.container = container

        menubar.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="Open file", command=self.openFile)
        filemenu.add_command(label="Save", command=self.saveFile)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=quit)

        menubar.add_cascade(label="Help")

        tk.Tk.config(self, menu=menubar)

        # self.current_frame = None
        # self.frame = {}
        frame = SettingFrame(container, self)
        # self.frame[SettingFrame] = frame
        frame.grid(row=0, column=0, sticky="nsew")
        # for Fr in {MainPage, CubicSplinePage, BSplinePage, PolynomialApprochPage}:
        #     frame = Fr(container, self)
        #     self.frame[Fr] = frame
        #     frame.grid(row=0, column=1, sticky="nsew")  # north south east west
        #     frame.grid_remove()

        frame = MainPage(self.container, self)
        frame.grid(row=0, column=1, sticky="nsew")  # north south east west
        
        self.show_frame(MainPage, 1)

    def openFile(self):
        ftypes = [('JSON files', '*.json'), ('All files', '*')]
        dlg = Open(self, filetypes=ftypes)
        fl = dlg.show()
        if fl != '':
            text = self.readFile(fl)
            if text != '':
                data = json.loads(text)
                global a, b, n, R, M
                a = data["a"]
                b = data["b"]
                R = data["R"]
                M = data["M"]
                build()
        else:            
            messagebox.showerror(title="Error", message="Erorr file!")

    def readFile(self, filename):
        f = open(filename, "r")
        text = f.read()
        if text == '':
            messagebox.showerror(title="Error", message="Empty file!")
            return ''
        else:
            return text

    def saveFile(self):
        file = filedialog.asksaveasfile(initialdir="C:\\Users\\duyenNH\\OneDrive\\Documents\\D_Python\\Spline\\data",
                                        mode='w',
                                        defaultextension='.json',
                                        filetype=[("JSON file", ".json"), ("Text file", ".txt"), ("All files", ".*")])
        if file is None:
            messagebox.showerror(title="Error", message="Save file unsuccessful!")
            return
        data = {
            'a': a,
            'b': b,
            'n': n,
            'R': R,
            'M': M
        }
        text = json.dumps(data)
        print(text)
        file.write(text)
        file.close()

    def show_frame(self, control, mode):
        # 1 - MainPage, 2 - B-spline, 3 - Cubic spline, 4 - InterPoly, 5 - ApprPoly
        global mode_app
        mode_app = mode

        global namePage
        if mode_app == 1:
            namePage = "Main Page"
        elif mode_app == 2:
            namePage = "B-spline Page"
        elif mode_app == 3:
            namePage = "Cubic Spline Page"
        else:
            namePage = "Polynomial Approximation page"
        if ani is not None:
            ani.resume()
        # frame = self.frame[control]
        # frame.tkraise()     # overlap (to stack) frames on top of each other
        # if self.current_frame is not None:
        #     self.current_frame.grid_remove()
        #     self.current_frame = None
        # frame = control(self.container, self)
        # frame.grid(row=0, column=1, sticky="nsew")  # north south east west
        # self.current_frame = frame

class MainPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        
        # label = tk.Label(self, text=namePage, font=LARGE_FONT)          
        # label.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        button = ttk.Button(self, text="Main Page",
                             command=lambda: controller.show_frame(self, 1))
        button.pack(side=tk.TOP, fill=tk.X, expand=True)

        button1 = ttk.Button(self, text="B-spline",
                             command=lambda: controller.show_frame(self, 2))
        button1.pack(side=tk.TOP, fill=tk.X, expand=True)

        button2 = ttk.Button(self, text="Cubic Spline",
                            command=lambda: controller.show_frame(self, 3))
        button2.pack(side=tk.TOP, fill=tk.X, expand=True)

        button3 = ttk.Button(self, text="Polynomial Approximation",
                             command=lambda: controller.show_frame(self, 4))
        button3.pack(side=tk.TOP, fill=tk.X, expand=True)

        canvas = FigureCanvasTkAgg(fig, self)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

class SettingFrame(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        label = tk.Label(self, text="Setting parameters", font=LARGE_FONT)
        label.grid(row=0, column=0, columnspan=3, padx=10, pady=10, sticky="WE")

        label = tk.Label(self, text="Boundaries:", font=SMALL_FONT)
        label.grid(row=1, column=0, padx=5, pady=5, sticky="WE")
        
        def getLeftBoundary(self):
            global a, x, y, x_vals
            a = getEntryText(a_entry)            
            x = np.linspace(a, b, n)
            y = np.array([func(i) for i in x])
            x_vals = np.linspace(min(x), max(x), 1000) 
            build()

        def getRightBoundary(self):
            global b, x, y, x_vals
            b = getEntryText(b_entry)            
            x = np.linspace(a, b, n)
            y = np.array([func(i) for i in x])
            x_vals = np.linspace(min(x), max(x), 1000)
            build()

        a_entry = tk.Entry(self, width=10, justify=tk.CENTER)
        a_entry.insert(tk.END, str(a))
        a_entry.bind('<Return>', getLeftBoundary)
        a_entry.grid(row=1, column=1)
        b_entry = tk.Entry(self, width=10, justify=tk.CENTER)
        b_entry.insert(tk.END, str(b))
        b_entry.bind('<Return>', getRightBoundary)
        b_entry.grid(row=1, column=2)

        N_lbl = tk.Label(self, text="Number points: ", font=SMALL_FONT)
        N_lbl.grid(row=3, column=0, padx=5, pady=5, sticky="WE")
        N_scale = tk.Scale(self, from_=1, to=1000, orient=tk.HORIZONTAL, 
                        activebackground="#32a88f", resolution=10, command=changeNumberOfPoint)
        N_scale.set(n)
        N_scale.grid(row=3, column=1, columnspan=2, padx=5, pady=5, sticky="WE")

        x_lbl = tk.Label(self, text="X: ", font=SMALL_FONT)
        x_lbl.grid(row=4, column=0, padx=5, pady=5, sticky="WE")
        x_entry = tk.Entry(self, width="20", bd=2)
        x_entry.grid(row=4, column=1, padx=5, pady=5, sticky="WE")

        y_lbl = tk.Label(self, text="Y: ", font=SMALL_FONT)
        y_lbl.grid(row=5, column=0, padx=5, pady=5, sticky="WE")
        y_entry = tk.Entry(self, width="20", bd=2)
        y_entry.grid(row=5, column=1, padx=5, pady=5, sticky="WE")

        def isNumber():
            try:
                float(x_entry.get())
                float(y_entry.get())
                print("par number is added")
                addPointToGraph(float(x_entry.get()), float(y_entry.get()))
                build()
            except ValueError:
                print("Validate an Entry! You need to entry a number!")
                messagebox.showerror(title="Error", message="Invalid value! \nYou need to entry a number!")

        btn_add = ttk.Button(self, text="Add", command=isNumber)
        btn_add.grid(row=6, column=1, padx=5, pady=5, sticky="WE")

        label = tk.Label(self, text="Coefficient Approximation B-spline", font=NORM_FONT)
        label.grid(row=7, column=0, columnspan=3, padx=10, pady=10)

        hr_lbl = tk.Label(self, text="Delta hR: "+ str(h_R), font=SMALL_FONT)
        hr_lbl.grid(row=9, column=0, padx=5, sticky="SNWE")
        def changeCoeffBSP1(valScaleR):
            global R, h_R
            R = float(valScaleR)
            if h*R < (a+b)/2:
                h_R = h*R
                hr_lbl.config(text='Delta hR: ' + str(round(h_R, 2)))
                build()
            else:
                messagebox.showerror(title="Error", message="Value R over range! \nNot allow value!")
        
        R_lbl = tk.Label(self, text="R: ", font=SMALL_FONT)
        R_lbl.grid(row=8, column=0, padx=5, sticky="WE")
        R_scale = tk.Scale(self, from_=1, to=50, orient=tk.HORIZONTAL,
                    activebackground="#32a88f", tickinterval=10, command=changeCoeffBSP1)
        R_scale.set(R)
        R_scale.grid(row=8, column=1, columnspan=2, padx=5, sticky="WE")

        # parameter for polynomial interpolation
        label = tk.Label(self, text="Coefficient Approximation Polynomial", font=NORM_FONT)
        label.grid(row=10, column=0, columnspan=3,padx=10, pady=10, sticky="NSWE")
        M_lbl = tk.Label(self, text="M: ", font=SMALL_FONT)
        M_lbl.grid(row=11, column=0, padx=5, sticky="WE")
        M_scale = tk.Scale(self, from_=1, to=50, orient=tk.HORIZONTAL,
                           activebackground="#32a88f", tickinterval=10, command=changeDegreeM)
        M_scale.set(M)
        M_scale.grid(row=11, column=1, columnspan=2, padx=5, sticky="WE")

build()
app = BSPApp()
app.title('BSP_APP')
app.geometry("1000x650")  # 1280x720

ani = animation.FuncAnimation(fig, animate, interval=100)
ani.pause()

app.mainloop()