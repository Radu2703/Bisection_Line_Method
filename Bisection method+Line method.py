import sys
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import figure
from sympy import *
import numpy as np
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication, QLabel, QWidget, QPushButton, QVBoxLayout, QMainWindow, QGridLayout, QLineEdit, QRadioButton, QButtonGroup
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtCore import QRegExp
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import math
from sympy import Symbol, diff, lambdify
from scipy.optimize import fsolve
from time import perf_counter, process_time, time, sleep
from matplotlib.animation import FuncAnimation, writers
import pandas as pd
import json
matplotlib.use('Qt5Agg')

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=5, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
        self.axes.axhline(y=0, color='black', linestyle='-')
        
def my_bisection(f, a, b, tol, r, win): 
    # approximates a root, R, of f bounded 
    # by a and b to within tolerance 
    # | f(m) | < tol with m the midpoint 
    # between a and b Recursive implementation
    # check if a and b bound a root
    if np.sign(f(a)) == np.sign(f(b)):
        win.lab6.setText("The scalars a and b do not bound a single root")
        raise Exception("The scalars a and b do not bound a single root")
    # get midpoint
    m = (a + b)/2  #=x0
    win.sc.axes.plot(m,f(m), 'o', color='b', linewidth=1)
    win.sc.axes.plot(m,0, 'o', color='r', linewidth=1)
    if len(win.vbsx.text())>0:
       win.vbsx.setText(win.vbsx.text()+"; "+str(m))
       win.vbsy.setText(win.vbsy.text()+"; "+str(f(m)))
    else :
        win.vbsx.setText(str(m))
        win.vbsy.setText(str(f(m)))
    if np.abs(m-r[0]) < tol:
        # stopping condition, report m as root
        return m
    elif np.sign(f(a)) == np.sign(f(m)):
        # case where m is an improvement on a. 
        # Make recursive call with a = m
        return my_bisection(f, m, b, tol, r, win)
    elif np.sign(f(b)) == np.sign(f(m)):
        # case where m is an improvement on b. 
        # Make recursive call with b = m
        return my_bisection(f, a, m, tol, r, win)
    
def my_bisection_fig(f, a, b, tol, r, axis):
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception("The scalars a and b do not bound a single root")
    # get midpoint
    m = (a + b)/2
    axis.plot(m,f(m), 'o', color='b', linewidth=1)
    axis.plot(m,0, 'o', color='r', linewidth=1)
    while np.abs(m-r[0])>=tol:
      if np.sign(f(a)) == np.sign(f(m)):
       a=m
       m = (a + b)/2
       axis.plot(m,f(m), 'o', color='b', linewidth=1)
       axis.plot(m,0, 'o', color='r', linewidth=1)
      elif np.sign(f(b)) == np.sign(f(m)):
       b=m
       m = (a + b)/2
       axis.plot(m,f(m), 'o', color='b', linewidth=1)
       axis.plot(m,0, 'o', color='r', linewidth=1)
    
def my_bisection_iter(f, a, b, win, iter, count=0): 
    # approximates a root, R, of f bounded 
    # by a and b to within tolerance 
    # | f(m) | < tol with m the midpoint 
    # between a and b Recursive implementation
    # check if a and b bound a root
    if np.sign(f(a)) == np.sign(f(b)):
        win.lab6.setText("The scalars a and b do not bound a single root")
        raise Exception("The scalars a and b do not bound a single root")
    # get midpoint
    m = (a + b)/2
    win.sc.axes.plot(m, f(m), 'o', color='b', linewidth=1)
    win.sc.axes.plot(m, 0, 'o', color='r', linewidth=1)
    if len(win.vbsx.text())>0:
       win.vbsx.setText(win.vbsx.text()+"; "+str(m))
       win.vbsy.setText(win.vbsy.text()+"; "+str(f(m)))
    else :
        win.vbsx.setText(str(m))
        win.vbsy.setText(str(f(m)))
    if count==iter:
        # stopping condition, report m as root
        return m
    elif np.sign(f(a)) == np.sign(f(m)):
        # case where m is an improvement on a. 
        # Make recursive call with a = m
        count+=1
        return my_bisection_iter(f, m, b, win, iter, count)
    elif np.sign(f(b)) == np.sign(f(m)):
        # case where m is an improvement on b. 
        # Make recursive call with b = m
        count+=1
        return my_bisection_iter(f, a, m, win, iter, count)
    
def my_bisection_iter_fig(f, a, b, axis, iter, count=0): 
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception("The scalars a and b do not bound a single root")
    # get midpoint
    m = (a + b)/2
    axis.plot(m, f(m), 'o', color='b', linewidth=1)
    axis.plot(m, 0, 'o', color='r', linewidth=1)
    while count!=iter:
      if np.sign(f(a)) == np.sign(f(m)):
        count+=1
        a=m
        m=(a + b)/2
        axis.plot(m, f(m), 'o', color='b', linewidth=1)
        axis.plot(m, 0, 'o', color='r', linewidth=1)
      elif np.sign(f(b)) == np.sign(f(m)):
        count+=1
        b=m
        m=(a + b)/2
        axis.plot(m, f(m), 'o', color='b', linewidth=1)
        axis.plot(m, 0, 'o', color='r', linewidth=1)

def my_line(f,a,b,tol,r,win): #f=string
    x=Symbol('x')
    f1=diff(f,x)
    f2=diff(f1,x)
    f=lambdify(x,f)
    f2=lambdify(x,f2)
    if np.sign(f(a)) == np.sign(f(b)):
        win.lab6.setText("The scalars a and b do not bound a single root")
        raise Exception("The scalars a and b do not bound a single root")
    if f(a)*f2(a)<0:
        x=a
        win.sc.axes.plot(x, f(x), 'o', color='b', linewidth=1)
        win.sc.axes.plot(x, 0, 'o', color='r', linewidth=1)
        win.vbsx.setText(str(x))
        win.vbsy.setText(str(f(x)))
        while abs(x-r[0])>tol:
            x=x-(f(x)/(f(x)-f(b)))*(x-b)
            win.sc.axes.plot(x, f(x), 'o', color='b', linewidth=1)
            win.sc.axes.plot(x, 0, 'o', color='r', linewidth=1)
            win.vbsx.setText(win.vbsx.text()+"; "+str(x))
            win.vbsy.setText(win.vbsy.text()+"; "+str(f(x)))
    else:
        x=b
        win.sc.axes.plot(x, f(x), 'o', color='b', linewidth=1)
        win.sc.axes.plot(x, 0, 'o', color='r', linewidth=1)
        win.vbsx.setText(str(x))
        win.vbsy.setText(str(f(x)))
        while abs(x-r[0])>tol:
            x=x-(f(x)/(f(x)-f(a)))*(x-a)
            win.sc.axes.plot(x, f(x), 'o', color='b', linewidth=1)
            win.sc.axes.plot(x, 0, 'o', color='r', linewidth=1)
            win.vbsx.setText(win.vbsx.text()+"; "+str(x))
            win.vbsy.setText(win.vbsy.text()+"; "+str(f(x)))
    return x

def my_line_fig(f,a,b,tol,r,axis):
    x=Symbol('x')
    f1=diff(f,x)
    f2=diff(f1,x)
    f=lambdify(x,f)
    f2=lambdify(x,f2)
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception("The scalars a and b do not bound a single root")
    if f(a)*f2(a)<0:
        x=a
        axis.plot(x, f(x), 'o', color='b', linewidth=1)
        axis.plot(x, 0, 'o', color='r', linewidth=1)
        while abs(x-r[0])>tol:
            x=x-(f(x)/(f(x)-f(b)))*(x-b)
            axis.plot(x, f(x), 'o', color='b', linewidth=1)
            axis.plot(x, 0, 'o', color='r', linewidth=1)
    else:
        x=b
        axis.plot(x, f(x), 'o', color='b', linewidth=1)
        axis.plot(x, 0, 'o', color='r', linewidth=1)
        while abs(x-r[0])>tol:
            x=x-(f(x)/(f(x)-f(a)))*(x-a)
            axis.plot(x, f(x), 'o', color='b', linewidth=1)
            axis.plot(x, 0, 'o', color='r', linewidth=1)

def my_line_iter(f,a,b,count,win): #f=string
    x=Symbol('x')
    f1=diff(f,x)
    f2=diff(f1,x)
    f=lambdify(x,f)
    f2=lambdify(x,f2)
    if np.sign(f(a)) == np.sign(f(b)):
        win.lab6.setText("The scalars a and b do not bound a single root")
        raise Exception("The scalars a and b do not bound a single root")
    if f(a)*f2(a)<0:
        x=a
        win.sc.axes.plot(x, f(x), 'o', color='b', linewidth=1)
        win.sc.axes.plot(x, 0, 'o', color='r', linewidth=1)
        win.vbsx.setText(str(x))
        win.vbsy.setText(str(f(x)))
        for i in range(count):
            x=x-(f(x)/(f(x)-f(b)))*(x-b)
            win.sc.axes.plot(x, f(x), 'o', color='b', linewidth=1)
            win.sc.axes.plot(x, 0, 'o', color='r', linewidth=1)
            win.vbsx.setText(win.vbsx.text()+"; "+str(x))
            win.vbsy.setText(win.vbsy.text()+"; "+str(f(x)))
    else:
        x=b
        win.sc.axes.plot(x, f(x), 'o', color='b', linewidth=1)
        win.sc.axes.plot(x, 0, 'o', color='r', linewidth=1)
        win.vbsx.setText(str(x))
        win.vbsy.setText(str(f(x)))
        for i in range(count):
            x=x-(f(x)/(f(x)-f(a)))*(x-a)
            win.sc.axes.plot(x, f(x), 'o', color='b', linewidth=1)
            win.sc.axes.plot(x, 0, 'o', color='r', linewidth=1)
            win.vbsx.setText(win.vbsx.text()+"; "+str(x))
            win.vbsy.setText(win.vbsy.text()+"; "+str(f(x)))
    return x

def my_line_iter_fig(f,a,b,count,axis):
    x=Symbol('x')
    f1=diff(f,x)
    f2=diff(f1,x)
    f=lambdify(x,f)
    f2=lambdify(x,f2)
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception("The scalars a and b do not bound a single root")
    if f(a)*f2(a)<0:
        x=a
        axis.plot(x, f(x), 'o', color='b', linewidth=1)
        axis.plot(x, 0, 'o', color='r', linewidth=1)
        for i in range(count):
            x=x-(f(x)/(f(x)-f(b)))*(x-b)
            axis.plot(x, f(x), 'o', color='b', linewidth=1)
            axis.plot(x, 0, 'o', color='r', linewidth=1)
    else:
        x=b
        axis.plot(x, f(x), 'o', color='b', linewidth=1)
        axis.plot(x, 0, 'o', color='r', linewidth=1)
        for i in range(count):
            x=x-(f(x)/(f(x)-f(a)))*(x-a)
            axis.plot(x, f(x), 'o', color='b', linewidth=1)
            axis.plot(x, 0, 'o', color='r', linewidth=1)

class Window(QMainWindow):
    #constructor
    def __init__(self,parent=None):
     super(Window,self).__init__(parent)
     self.setWindowTitle("Approximating solutions of non-linear equations (Bisection method and Line method)")
     self.resize=(500,600)
     widget=QWidget()
     widget.resize=(500,1000)
     layout=QGridLayout()
     self.lab1=QLabel("f(x):[a,b] = ")
     self.lab2=QLabel("a = ")
     self.lab3=QLabel("b = ")
     self.lab4=QLabel("Error toleration = ")
     self.lab5=QLabel("Number of iterations = ")
     self.lab6=QLabel()
     self.lab7=QLabel("Approximative value of xn = ")
     self.lab8=QLabel("Approximative value of yn = ")
     self.labfx=QLabel("Exact value of xn (fsolve)= ")
     self.labfy=QLabel("Exact value of yn (fsolve)= ")
     self.lab9=QLabel("Error for abs(xn-fsolve) = ")
     self.lab10=QLabel("Error for abs(yn-fsolve) = ")
     self.lab11=QLabel("Performance time for approximative xn = ")
     self.lab12=QLabel("Performance time for exact xn = ")
     self.lab13=QLabel("Process time for approximative xn = ")
     self.lab14=QLabel("Process time for for exact xn = ")
     self.lab15=QLabel("Performance time for approximative yn = ")
     self.lab16=QLabel("Performance time for exact yn = ")
     self.lab17=QLabel("Process time for approximative yn = ")
     self.lab18=QLabel("Process time for exact yn = ")
     self.labvbs=QLabel("Set of the tested variables (x and y)")
     self.vbsx=QLineEdit()
     self.vbsy=QLineEdit()
     layout.addWidget(self.lab1,0,0)
     layout.addWidget(self.lab2,1,0)
     layout.addWidget(self.lab3,2,0)
     layout.addWidget(self.lab4,4,1)
     layout.addWidget(self.lab5,5,1)
     layout.addWidget(self.labvbs,6,0)
     layout.addWidget(self.vbsx,6,1)
     layout.addWidget(self.vbsy,6,2)
     layout.addWidget(self.lab6,9,6)
     layout.addWidget(self.lab7,7,0)
     layout.addWidget(self.lab8,8,0)
     layout.addWidget(self.labfx,9,0)
     layout.addWidget(self.labfy,10,0)
     layout.addWidget(self.lab9,11,0)
     layout.addWidget(self.lab10,12,0)
     layout.addWidget(self.lab11,13,0)
     layout.addWidget(self.lab12,14,0)
     layout.addWidget(self.lab13,15,0)
     layout.addWidget(self.lab14,16,0)
     layout.addWidget(self.lab15,17,0)
     layout.addWidget(self.lab16,18,0)
     layout.addWidget(self.lab17,19,0)
     layout.addWidget(self.lab18,20,0)
     self.fct1=QLineEdit()
     self.fct1.setPlaceholderText("Enter function")
     self.fct2=QLineEdit()
     self.fct2.setPlaceholderText("Enter start of interval")
     self.fct3=QLineEdit()
     self.fct3.setPlaceholderText("Enter end of interval")
     self.rb_b=QRadioButton()
     self.rb_b.setText("BISECTION METHOD")
     self.rb_l=QRadioButton()
     self.rb_l.setText("LINE METHOD")
     self.er=QLineEdit()
     self.er.setPlaceholderText("Enter tolerated error")
     self.iter=QLineEdit()
     self.iter.setPlaceholderText("Enter number of iterations")
     validator1=QRegExpValidator(QRegExp(r".{0,}x.{0,}"))
     self.fct1.setValidator(validator1)
     validator2=QRegExpValidator(QRegExp(r"\-?\d+\.\d*"))
     validator3=QRegExpValidator(QRegExp(r"\d+\.\d*"))
     validator4=QRegExpValidator(QRegExp(r"\d*"))
     self.fct2.setValidator(validator2)
     self.fct3.setValidator(validator2)
     self.er.setValidator(validator3)
     self.iter.setValidator(validator4)
     self.er_b=QRadioButton()
     self.er_b.setText("Error method")
     self.iter_b=QRadioButton()
     self.iter_b.setText("Iterative method")
     self.xn=QLineEdit()
     self.yn=QLineEdit()
     self.fx=QLineEdit()
     self.fy=QLineEdit()
     self.xe=QLineEdit()
     self.ye=QLineEdit()
     self.pemx=QLineEdit()
     self.pefx=QLineEdit()
     self.prmx=QLineEdit()
     self.prfx=QLineEdit()
     self.pemy=QLineEdit()
     self.pefy=QLineEdit()
     self.prmy=QLineEdit()
     self.prfy=QLineEdit()
     self.file=QLineEdit()
     self.file.setPlaceholderText("Enter file path and name (with extension)")
     self.xn.setEnabled(False)
     self.yn.setEnabled(False)
     self.fx.setEnabled(False)
     self.fy.setEnabled(False)
     self.xe.setEnabled(False)
     self.ye.setEnabled(False)
     self.pemx.setEnabled(False)
     self.pefx.setEnabled(False)
     self.prmx.setEnabled(False)
     self.prfx.setEnabled(False)
     self.pemy.setEnabled(False)
     self.pefy.setEnabled(False)
     self.prmy.setEnabled(False)
     self.prfy.setEnabled(False)
     self.vbsx.setEnabled(False)
     self.vbsy.setEnabled(False)
     self.method_group=QButtonGroup(widget)
     self.method_group.addButton(self.rb_b)
     self.method_group.addButton(self.rb_l)
     self.type_group=QButtonGroup(widget)
     self.type_group.addButton(self.er_b)
     self.type_group.addButton(self.iter_b)
     layout.addWidget(self.fct1,0,1)
     layout.addWidget(self.fct2,1,1)
     layout.addWidget(self.fct3,2,1)
     layout.addWidget(self.rb_b,3,0)
     layout.addWidget(self.rb_l,3,2)
     layout.addWidget(self.er,4,2)
     layout.addWidget(self.iter,5,2)
     layout.addWidget(self.er_b,4,0)
     layout.addWidget(self.iter_b,5,0)
     layout.addWidget(self.xn,7,2)
     layout.addWidget(self.yn,8,2)
     layout.addWidget(self.fx,9,2)
     layout.addWidget(self.fy,10,2)
     layout.addWidget(self.xe,11,2)
     layout.addWidget(self.ye,12,2)
     layout.addWidget(self.pemx,13,2)
     layout.addWidget(self.pefx,14,2)
     layout.addWidget(self.prmx,15,2)
     layout.addWidget(self.prfx,16,2)
     layout.addWidget(self.pemy,17,2)
     layout.addWidget(self.pefy,18,2)
     layout.addWidget(self.prmy,19,2)
     layout.addWidget(self.prfy,20,2)
     self.bplot=QPushButton('PLOT')
     layout.addWidget(self.bplot,7,6)
     self.sc = MplCanvas(self, width=5, height=5, dpi=100)
     layout.addWidget(self.sc,6,6)
     self.showtime=QPushButton('SHOW/HIDE PLOT')
     layout.addWidget(self.showtime,8,6)
     self.ran=QPushButton('GENERATE RANDOM NUMBERS (INTERVAL: BETWEEN -1000 AND 1000; b-a<=10)')
     layout.addWidget(self.ran,0,2)
     layout.addWidget(self.file,1,2)
     self.fread=QPushButton('READ FILE (.txt<->separate lines; .csv<->header: f,v1,v2; .json<->elements: f,v1,v2)')
     layout.addWidget(self.fread,2,2)
     widget.setLayout(layout)
     self.bplot.clicked.connect(self.plotFcn)
     self.showtime.clicked.connect(self.gph)
     self.ran.clicked.connect(self.randnbs)
     self.fread.clicked.connect(self.fileread)
     self.setCentralWidget(widget)
     self.setGeometry(0,0,10000,1000)
     self.show()
    def plotFcn(self):
     try:
      myfun=self.fct1.text()
      x=Symbol('x') #your symbolic value
      myfunl=lambdify(x,myfun)
      if len(self.fct1.text())>0 and "x" in self.fct1.text() and len(self.fct2.text())>0 and len(self.fct3.text())>0 and float(self.fct2.text())<float(self.fct3.text()):
        if (self.rb_b.isChecked()==True or self.rb_l.isChecked()==True) and ((self.er_b.isChecked()==True and len(self.er.text())>0 and float(self.er.text())>0) or (self.iter_b.isChecked()==True and len(self.iter.text())>0 and int(self.iter.text())>0)):
          #metadata1=dict(title="Bisection Error Method",artist='Matplotlib',comment='Bisection Error Method')
          #Writer=writers['ffmpeg']
          #writer1=Writer(fps=1,metadata=metadata1)
          #metadata2=dict(title="Bisection Iterative Method",artist='Matplotlib',comment='Bisection Iterative Method')
          #writer2=Writer(fps=1,metadata=metadata2)
          #metadata3=dict(title="Line Error Method",artist='Matplotlib',comment='Line Error Method')
          #writer3=Writer(fps=1,metadata=metadata3)
          #metadata4=dict(title="Line Iterative Method",artist='Matplotlib',comment='Line Iterative Method')
          #writer4=Writer(fps=1,metadata=metadata4)
          self.sc.axes.cla()
          plt.clf()
          self.vbsx.setText("")
          self.vbsy.setText("")
          self.xn.setText("")
          self.yn.setText("")
          self.fx.setText("")
          self.fy.setText("")
          self.xe.setText("")
          self.ye.setText("")
          self.pemx.setText("")
          self.pefx.setText("")
          self.prmx.setText("")
          self.prfx.setText("")
          self.pemy.setText("")
          self.pefy.setText("")
          self.prmy.setText("")
          self.prfy.setText("")
          pex2=perf_counter()
          prx2=process_time()
          r=fsolve(myfunl,(float(self.fct2.text())+float(self.fct3.text()))/2)
          pex2=perf_counter()-pex2
          prx2=process_time()-prx2
          pey2=perf_counter()
          pry2=process_time()
          q=myfunl(r[0].astype(float))
          pey2=perf_counter()-pey2
          pry2=process_time()-pry2
          a=np.linspace(float(self.fct2.text()),float(self.fct3.text()),10)
          self.sc.axes.plot(a, myfunl(a))
          self.sc.axes.set_xlim(float(self.fct2.text()), float(self.fct3.text()))
          self.sc.axes.set_xticks(a)
          self.sc.axes.axhline(y=0, color='black', linestyle='-')
          #mov=plt.figure()
          axis=plt.axes()
          axel=np.linspace(float(self.fct2.text()),float(self.fct3.text()),10)
          axis.plot(axel, myfunl(axel))
          axis.set_xlim(float(self.fct2.text()), float(self.fct3.text()))
          axis.set_xticks(axel)
          axis.axhline(y=0, color='black', linestyle='-')
          plt.xlabel('Interval: ['+self.fct2.text()+"; "+self.fct3.text()+"]")
          plt.ylabel("Function: "+self.fct1.text())
          if self.rb_b.isChecked()==True and self.er_b.isChecked()==True and len(self.er.text())>0 and float(self.er.text())>0:
              pex1=perf_counter()
              prx1=process_time()
              rm=my_bisection(myfunl, float(self.fct2.text()), float(self.fct3.text()), float(self.er.text()), r, self)
              pex1=perf_counter()-pex1
              prx1=process_time()-prx1
              pey1=perf_counter()
              pry1=process_time()
              qm=myfunl(rm)
              pey1=perf_counter()-pey1
              pry1=process_time()-pry1
              plt.title('Bisection Error Method ('+self.er.text()+')')
              my_bisection_fig(myfunl, float(self.fct2.text()), float(self.fct3.text()), float(self.er.text()), r, axis)
              #anim=FuncAnimation(mov, func=my_bisection_fig(myfunl, float(self.fct2.text()), float(self.fct3.text()), float(self.er.text()), r, axis), interval=10, frames=np.arange(1,int(self.vbsx.text().count(';')+2)))
              #anim.save('BisectionErrorMethod.mp4',writer1)
          if self.rb_b.isChecked()==True and self.iter_b.isChecked()==True and len(self.iter.text())>0 and int(self.iter.text())>0:
              pex1=perf_counter()
              prx1=process_time()
              rm=my_bisection_iter(myfunl, float(self.fct2.text()), float(self.fct3.text()), self, int(self.iter.text()))
              pex1=perf_counter()-pex1
              prx1=process_time()-prx1
              pey1=perf_counter()
              pry1=process_time()
              qm=myfunl(rm)
              pey1=perf_counter()-pey1
              pry1=process_time()-pry1
              plt.title('Bisection Iterative Method ('+self.iter.text()+')')
              my_bisection_iter_fig(myfunl, float(self.fct2.text()), float(self.fct3.text()), axis, int(self.iter.text()))
              #anim=FuncAnimation(mov, func=my_bisection_iter_fig(myfunl, float(self.fct2.text()), float(self.fct3.text()), axis, int(self.iter.text())), interval=10, frames=np.arange(1,int(self.iter.text())+1))
              #anim.save('BisectionIterativeMethod.mp4',writer2)
          if self.rb_l.isChecked()==True and self.er_b.isChecked()==True and len(self.er.text())>0 and float(self.er.text())>0:
              pex1=perf_counter()
              prx1=process_time()
              rm=my_line(myfun, float(self.fct2.text()), float(self.fct3.text()), float(self.er.text()), r, self)
              pex1=perf_counter()-pex1
              prx1=process_time()-prx1
              pey1=perf_counter()
              pry1=process_time()
              qm=myfunl(rm)
              pey1=perf_counter()-pey1
              pry1=process_time()-pry1
              plt.title('Line Error Method ('+self.er.text()+')')
              my_line_fig(myfun, float(self.fct2.text()), float(self.fct3.text()), float(self.er.text()), r, axis)
              #anim=FuncAnimation(mov, func=my_line_fig(myfun, float(self.fct2.text()), float(self.fct3.text()), float(self.er.text()), r, axis), interval=10, frames=np.arange(1,int(self.vbsx.text().count(';')+2)))
              #anim.save('LineErrorMethod.mp4',writer3)
          if self.rb_l.isChecked()==True and self.iter_b.isChecked()==True and len(self.iter.text())>0 and float(self.iter.text())>0:
              pex1=perf_counter()
              prx1=process_time()
              rm=my_line_iter(myfun, float(self.fct2.text()), float(self.fct3.text()), int(self.iter.text()), self)
              pex1=perf_counter()-pex1
              prx1=process_time()-prx1
              pey1=perf_counter()
              pry1=process_time()
              qm=myfunl(rm)
              pey1=perf_counter()-pey1
              pry1=process_time()-pry1
              plt.title('Line Iterative Method ('+self.iter.text()+')')
              my_line_iter_fig(myfun, float(self.fct2.text()), float(self.fct3.text()), int(self.iter.text()), axis)
              #anim=FuncAnimation(mov, func=my_line_iter_fig(myfun, float(self.fct2.text()), float(self.fct3.text()), int(self.iter.text()), axis), interval=10, frames=np.arange(1,int(self.iter.text())+1))
              #anim.save('LineIterativeMethod.mp4',writer4)
          self.xn.setText(str(rm))
          self.yn.setText(str(qm))
          self.fx.setText(str(r[0].astype(float)))
          self.fy.setText(str(q))
          self.xe.setText(str(math.fabs(rm-r[0].astype(float))))
          self.ye.setText(str(math.fabs(qm-q)))
          self.pemx.setText(str(pex1))
          self.pefx.setText(str(pex2))
          self.prmx.setText(str(prx1))
          self.prfx.setText(str(prx2))
          self.pemy.setText(str(pey1))
          self.pefy.setText(str(pey2))
          self.prmy.setText(str(pry1))
          self.prfy.setText(str(pry2))
          self.vbsx.setEnabled(True)
          self.vbsy.setEnabled(True)
          if np.isnan(rm):
              self.lab6.setText("The function's number of roots in the specified interval is different from one")
              self.sc.axes.cla()
              self.vbsx.setText("")
              self.vbsy.setText("")
              self.xn.setText("")
              self.yn.setText("")
              self.fx.setText("")
              self.fy.setText("")
              self.xe.setText("")
              self.ye.setText("")
              self.pemx.setText("")
              self.pefx.setText("")
              self.prmx.setText("")
              self.prfx.setText("")
              self.pemy.setText("")
              self.pefy.setText("")
              self.prmy.setText("")
              self.prfy.setText("")
          else: self.lab6.setText("")
          self.sc.draw()
        else:
            if self.rb_b.isChecked()==False and self.rb_l.isChecked()==False:
               self.lab6.setText("Method wasn't chosen")
            elif self.er_b.isChecked()==False and self.iter_b.isChecked()==False:
               self.lab6.setText("Condition wasn't chosen")
            elif (self.er_b.isChecked()==True and len(self.er.text())<=0) or (self.iter_b.isChecked()==True and len(self.iter.text())<=0):
               self.lab6.setText("Condition value wasn't specified")
            elif (self.er_b.isChecked()==True and len(self.er.text())>0 and float(self.er.text())<=0) or (self.iter_b.isChecked()==True and len(self.iter.text())>0 and int(self.iter.text())<=0):
               self.lab6.setText("Condition value needs to be strictly bigger than 0")
      else :
          if len(self.fct1.text())<=0:
             self.lab6.setText("Data value wasn't specified for function")
          elif len(self.fct2.text())<=0:
             self.lab6.setText("Data value wasn't specified for start of interval")
          elif len(self.fct3.text())<=0:
             self.lab6.setText("Data value wasn't specified for end of interval")
          elif float(self.fct2.text())>=float(self.fct3.text()):
             self.lab6.setText("Start of interval needs to be strictly smaller than end of interval")
          elif "x" not in self.fct1.text():
             self.lab6.setText("Function doesn't have the x variable")
             
     except:
      self.lab6.setText("Function is wrongly written, the function cannot be applied to the interval, or the number of roots in the specified interval is different from one ")
      self.sc.axes.cla()
      plt.clf()
      self.vbsx.setText("")
      self.vbsy.setText("")
      self.xn.setText("")
      self.yn.setText("")
      self.fx.setText("")
      self.fy.setText("")
      self.xe.setText("")
      self.ye.setText("")
      self.pemx.setText("")
      self.pefx.setText("")
      self.prmx.setText("")
      self.prfx.setText("")
      self.pemy.setText("")
      self.pefy.setText("")
      self.prmy.setText("")
      self.prfy.setText("")
      
    def gph(self):
      if self.sc.isVisible()==True:
          self.sc.hide()
      else:
          self.sc.show()
    
    def randnbs(self):
      self.fct2.setText(str(np.random.randint(-1000,1000)+np.random.rand(1)[0]))
      ranvb=float(self.fct2.text())
      self.fct3.setText(str(ranvb+np.random.randint(0,10)+np.random.rand(1)[0]))
      
    def fileread(self):
     try:
      if len(self.file.text())>0 and (".txt" in self.file.text() or ".csv" in self.file.text() or ".json" in self.file.text()):
       filee=repr(self.file.text().strip())[1:-1]
       file1=open(filee,"r")
       if ".txt" in self.file.text():
         lines=file1.readlines()
         lines[0].replace("\n","")
         lines[0]=lines[0].strip()
         lines[1].replace("\n","")
         lines[1]=lines[1].strip()
         self.fct1.setText(lines[0])
         self.fct2.setText(lines[1])
         self.fct3.setText(lines[2])
       elif ".csv" in self.file.text():
         my_data=pd.read_csv(file1)
         fcsv=my_data['f']
         fcsv=np.array(fcsv)
         acsv=my_data['v1']
         acsv=np.array(acsv)
         bcsv=my_data['v2']
         bcsv=np.array(bcsv)
         self.fct1.setText(fcsv[0])
         self.fct2.setText(str(acsv[0]))
         self.fct3.setText(str(bcsv[0]))
       elif ".json" in self.file.text():
         data=json.load(file1)
         fcsv=data['f']
         fcsv=np.array(fcsv)
         acsv=data['v1']
         acsv=np.array(acsv)
         bcsv=data['v2']
         bcsv=np.array(bcsv)
         self.fct1.setText(data['f'])
         self.fct2.setText(str(data['v1']))
         self.fct3.setText(str(data['v2']))
       file1.close()
      else: self.lab6.setText("File path wasn't entered, or file isn't of a specified type")
     except:
      self.lab6.setText("File doesn't exist")
      
if __name__=='__main__':
 app = QtWidgets.QApplication(sys.argv)
 main = Window()
 main.show()
 sys.exit(app.exec_())