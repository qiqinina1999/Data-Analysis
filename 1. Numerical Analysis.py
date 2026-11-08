"""
Project: Numerical Analysis
Description:
This project implements various numerical methods using Python as part of academic coursework.
"""

# Least-square methods

import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *


X=np.array([1,2,3,4,5,6,7,8,9])
Y=np.array([1,1.5,2,3,4,5,8,10,13])


X_ave=np.mean(X)
Y_ave=np.mean(Y)

num=0
den=0
for i in range(len(X)):
  num += (X[i]-X_ave)*(Y[i]-Y_ave)
  den += (X[i]-X_ave)**2
m=num/den
c=Y_ave-m*X_ave
print(m,c)
Yp = m * X + c
plt.scatter(X,Y)
plt.plot(X,Yp,color='black')
plt.title('linear regression')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

#Newtons

nth = int(input('nth_order_difference='))

x = np.array([-8,-6.22222222,-4.44444444,-2.66666667,-0.88888889,0.88888889,2.66666667,4.44444444,6.22222222,8])
y = np.array([5.94060012e-11,7.94829537e-08,2.61023666e-05,2.10401088e-03,4.16273041e-02,2.02148629e-01,2.40949220e-01,7.04924261e-02,5.06199370e-03,8.92201505e-05])

for i in range(nth):
    i = 0
    Newtons =[0]*len(x)
    while i < nth:
        j = nth
        while j > i:
            if i == 0:
                Newtons[j]=((y[j]-y[j-1])/(x[j]-x[j-1]))
            else:
                Newtons[j] = (Newtons[j]-Newtons[j-1])/(x[j]-x[j-1-i])
            j -= 1
        i+= 1
print(Newtons)
def function(data):
    return y[0]+Newtons[1]*(data-x[0])+Newtons[2]*(data-x[0])*(data-x[1])+\
    Newtons[3]*(data-x[0])*(data-x[1])*(data-x[2])+\
    Newtons[4]*(data-x[0])*(data-x[1])*(data-x[2])*(data-x[3])+\
    Newtons[5]*(data-x[0])*(data-x[1])*(data-x[2])*(data-x[3])*(data-x[4])*(data-x[5])+\
    Newtons[6]*(data-x[0])*(data-x[1])*(data-x[2])*(data-x[3])*(data-x[4])*(data-x[5])*(data-x[6])+\
    Newtons[7]*(data-x[0])*(data-x[1])*(data-x[2])*(data-x[3])*(data-x[4])*(data-x[5])*(data-x[6])*(data-x[7])+\
    Newtons[8]*(data-x[0])*(data-x[1])*(data-x[2])*(data-x[3])*(data-x[4])*(data-x[5])*(data-x[6])*(data-x[7])*(data-x[8])+\
    Newtons[9]*(data-x[0])*(data-x[1])*(data-x[2])*(data-x[3])*(data-x[4])*(data-x[5])*(data-x[6])*(data-x[7])*(data-x[8])*(data-x[9])


def calculate_data(x,Newtons):
    returnData=[];
    for data in x:
        returnData.append(function(data))
    return returnData

et=(function(8)-8.92201505e-05)/8.92201505e-05
print('Newtons interpolating polynomial=',function(8))
print('Newtons interpolating polynomial error=',abs(et))


plt.plot(x,calculate_data(x,Newtons),color="orange")
plt.scatter(x,y)

#Lagrange

x = np.array([-8,-6.22222222,-4.44444444,-2.66666667,-0.88888889,0.88888889,2.66666667,4.44444444,6.22222222,8])
y = np.array([5.94060012e-11,7.94829537e-08,2.61023666e-05,2.10401088e-03,4.16273041e-02,2.02148629e-01,2.40949220e-01,7.04924261e-02,5.06199370e-03,8.92201505e-05])


def Lagrange(x,y):
    n = len(x)
    m = len(y)
    p = np.poly1d([0])
    for k in range(n):
        L = np.poly1d([1])
        for i in range(n):
            if i!=k:
              L = L*np.poly1d([1.0/(x[k]-x[i]),-float(x[i]/(x[k]-x[i]))])
        p = p+y[k]*L
    return p

p = Lagrange(x,y)

xi = 4
print('f(4)的值為',p(xi))
plt.plot(x,y)
plt.scatter(xi,p(xi),color='red')
plt.title('3th order Lagrange interpolating polynomial')

#Cubic spline

x = np.array([-8,-6.22222222,-4.44444444,-2.66666667,-0.88888889,0.88888889,2.66666667,4.44444444,6.22222222,8])
y = np.array([5.94060012e-11,7.94829537e-08,2.61023666e-05,2.10401088e-03,4.16273041e-02,2.02148629e-01,2.40949220e-01,7.04924261e-02,5.06199370e-03,8.92201505e-05])
n=len(x)-1
h=np.zeros(n)
for i in range(1,n+1):
  h[i-1]=x[i]-x[i-1]

u=np.zeros(n-1)
l=np.zeros(n-1)
for i in range(1,n):
  u[i-1]=h[i-1]/(h[i-1]+h[i])
  l[i-1]=1-u[i-1]

d=np.zeros(n)
for i in range(2):
  j=n
  while j>i:
    if i==0:
      d[j-1]=((y[j]-y[j-1])/(x[j]-x[j-1]))
    else:
      d[j-1]=(d[j-1]-d[j-2])/(x[j]-x[j-1-i])
    j-=1

d1=np.zeros(n-1)
for i in range(n-1):
  d1[i]=6*d[i+1]

A=np.zeros([n-1,n-1])
for i in range(1,n-2):
  A[i,i-1]=u[i]
  A[i,i]=2
  A[i,i+1]=l[i]
A[0,0]=A[2,2]=2
A[0,1]=l[0]
A[2,1]=u[2]
A[7,7]=2

M0=np.array(n-1)
M0=np.linalg.solve(A,d1)
M=np.zeros(n+1)
for i in range(1,n):
  M[i]=M0[i-1]

F=[0]*len(x)
m=symbols('x')
for i in range(1,n+1):
  S=M[i-1]*((x[i]-m)**3)/(6*h[i-1])+M[i]*((m-x[i-1])**3)/(6*h[i-1])+(y[i-1]-(h[i-1]**2)*M[i-1]/6)*(x[i]-m)/h[i-1]+(y[i]-(h[i-1]**2)*M[i]/6)*(m-x[i-1])/h[i-1]
  F[i-1]=expand(S)
  print(x[i-1],'<',m,'<',x[i])
  print(F[i-1])

b=4.07191685267853e-6*27**3 + 0.00219893973214286*27**2 - 0.271907924107143*27 + 13.6209107142857
et2=(b-7.986)/7.986

print('cubic spline=',b)
print('cubic spline error=',abs(et2))

x1=np.linspace(-8,-6.22222222,100)
y1=-2.25032613751899e-5*x1**3 - 0.000540078273004558*x1**2 - 0.00424946008944587*x1 - 0.0109523410079665
x2=np.linspace(-6.22222222,-4.44444444,100)
y2=0.00011713368537716*x2**3 + 0.00206647806544172*x2**2 + 0.0119691126773164*x2 + 0.0226861802740453
x3=np.linspace(-4.44444444 ,-2.66666667,100)
y3=-8.54572986725475e-5*x3**3 - 0.000634735052519827*x3**2 - 3.62789460628508e-5*x3 + 0.00490041492386174
x4=np.linspace(-2.66666667,-0.88888889,100)
y4=0.00652397790759864*x4**3 + 0.052240746663744*x4**2 + 0.140965005806892*x4 + 0.13023489041649
x5=np.linspace(-0.88888889, 0.88888889,100)
y5=-0.0111398783016572*x5**3 + 0.00513713004684886*x5**2 + 0.0990951243173148*x5 + 0.117828999589256
x6=np.linspace(0.88888889, 2.66666667,100)
y6=-0.00516319015123525*x6**3 - 0.0108007050408654*x6**2 + 0.113262088857436*x6 + 0.11363138046101
x7=np.linspace(2.66666667, 4.44444444,100)
y7=0.0162129956041378*x7**3 - 0.181810191297611*x7**2 + 0.569287386112124*x7 - 0.291724439827629
x8=np.linspace(4.44444444, 6.22222222,100)
y8=-0.00375303700844996*x8**3 + 0.0844035766040112*x8**2 - 0.613884914489694*x8 + 1.46112341116296
x9=np.linspace(6.22222222, 8.0,100)
y9=-0.002690041085006*x9**3 + 0.0645609860401441*x9**2 - 0.510783215174026*x9 + 1.33175287049656

plt.plot(x1,y1)
plt.plot(x2,y2)
plt.plot(x3,y3)
plt.plot(x4,y4)
plt.plot(x5,y5)
plt.plot(x6,y6)
plt.plot(x7,y7)
plt.plot(x8,y8)
plt.plot(x9,y9)

"""微分"""

#Euler’s method

def Euler(f,x0,xf,y0=2,n=100):
    x = np.linspace(x0,xf,n)
    y = list([y0])
    for i in range(n-1):
        h = x[i+1]-x[i]
        y.append(y[-1]+h*f(x[i],y[-1]))

    y = np.array(y)
    return y

def f(x,y):
  return 4*np.exp(0.8*x) -0.5*y

y_Euler = Euler(lambda x,y:f(x,y),0,3,1,100)


#the analytical solution
def F(x):
  return 4*(np.exp(0.8*x)-np.exp(-0.5*x))/1.3 +2*np.exp(-0.5*x)

x1 = np.linspace(0,3,100)
a = F(x1)

plt.plot(x1,a,'b',x1,y_Euler,'r')

#Runge-Kutta method

def RK4(f,x0,xf,y0=2,n=100):
  x = np.linspace(x0,xf,n)
  y = list([y0])
  for i in range(n-1):
    h = x[i+1]-x[i]
    k1 = h*f(x[i],y[-1])
    k2 = h*f(x[i]+h/2.0,y[-1]+k1/2.0)
    k3 = h*f(x[i]+h/2.0,y[-1]+k2/2.0)
    k4 = h*f(x[i]+h,y[-1]+k3)
    y.append(y[-1]+(k1+2*k2+2*k3+k4)/6)

  y = np.array(y)
  return y

def f(x,y):
  return 4*np.exp(0.8*x) -0.5*y

y_RK4 = RK4(lambda x,y: f(x,y),0,3,1,100)

#the analytical solution
def F(x):
  return 4*(np.exp(0.8*x)-np.exp(-0.5*x))/1.3 +2*np.exp(-0.5*x)

x1 = np.linspace(0,3,100)
a = F(x1)

plt.plot(x1,a,'b',x1,y_RK4,'r')

"""積分"""

#simpon 1/3
def simpson(f,a,b,N):
  if N/2==0 :
    n=N
  else:
    n=N+1
  h=float(b-a)/n
  s=float(f(a))

  for i in range(1,n,2):
    x  = a+h*i
    s += 4*f(x)
  for i in range(2,n-1,2):
    x =  a+h*i
    s += 2*f(x)
  s += f(b)
  s =  h*s/3

  return s
def f(x):
  return x**2
a=0
b=1

simpson(lambda x:f(x),a,b,10**3)

def f1(x):
   f1 = x**4.0
   return f1;

def f2(x):
   f2 = np.exp(-x*x)
   return f2;

def f3(x):
    return np.sinc(x / np.pi)

# trapezoidal rule
def trapezoid(f,a,b,N):
    h   = (b-a)/N
    xi  = np.linspace(a,b,N+1)
    fi  = f(xi)
    s   = 0.0
    for i in range(1,N):
        s = s + fi[i]
    s = (h/2)*(fi[0] + fi[N]) + h*s
    return s

# romberg method
def romberg(f,a,b,eps,nmax):
# f     ... function to be integrated
# [a,b] ... integration interval
# eps   ... desired accuracy
# nmax  ... maximal order of Romberg method
    Q         = np.zeros((nmax,nmax),float)
    converged = 0
    for i in range(0,nmax):
        N      = 2**i
        Q[i,0] = trapezoid(f,a,b,N)
        for k in range(0,i):
            n        = k + 2
            Q[i,k+1] = 1.0/(4**(n-1)-1)*(4**(n-1)*Q[i,k] - Q[i-1,k])
        if (i > 0):
            if abs(Q[i,i] - Q[i-1,i-1]) < eps:
               converged = 1
               break
    print(Q[i,i],N,converged)
    return Q[i,k+1],N,converged

# main program
a  = 0.0;b = 1.0  # integration interval [a,b]
romberg(f1,a,b,1.0e-12,10)
romberg(f2,a,b,1.0e-12,10)
a  = 0.0;b = 20.0*np.pi  # integration interval [a,b]
romberg(f3,a,b,1.0e-12,10)

"""線性"""

#(a)
det=np.linalg.det
inv=np.linalg.inv

A=np.asmatrix([[15,-3,-1],[-3,18,-6],[-4,-1,12]])
B=np.asmatrix(([[3300],[1200],[2400]]))

def inverse(arr):
  if round(det(arr),2)==0:
    return 'no inverse'
  return inv(arr)
print(inverse(A))

x=inverse(A)*B
print(x)
