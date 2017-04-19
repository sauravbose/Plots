# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 13:29:06 2017

@author: User
"""
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.mlab import griddata
import numpy as np
from sympy import *
import math
import xlsxwriter


def read_data():
    global df, n1, n2, A_AT,nT, rho1star, rho2star, bias_box, full_box
    
    full_box = 6.0
    bias_box = 3.0 
    
    df = pd.read_csv('C:\Users\User\Dropbox\Saurav Bose\A_LJ_Plot.csv')
    
    nT = df['nT']
    nT=np.array(nT)
    
    n1 = df['n1']
    n1=np.array(n1)
    rho1star = n1/(bias_box**3)
    
    n2 = df['n2']
    n2 = np.array(n2)
    rho2star = n2/((full_box**3)-(bias_box**3))
    
    A_AT = df['A-AT']
    A_AT = np.array(A_AT)


def xiao_data_1_2():
    global rho1_mu, rho2_mu, mu_xiao_12,mu_id
    
    p00 = -0.420667966848815
    p10 = -0.561825932101550
    p01 = 2.23115817186585
    p20 = 1.25065441851983
    p11 = -10.3948786207038
    p02 = -14.7336522151647
    p30 = 0.889816425429155
    p21 = 4.28099293535048
    p12 = 14.5228143374848
    p03 = 10.7054914580212        

    x = Symbol('x')
    y = Symbol('y')
    z = p00 + p10*x + p01*y + p20*x**2 + p11*x*y + p02*y**2 + p30*x**3 + p21*x**2*y + p12*x*y**2 + p03*y**3
    
    MU = lambdify((x,y),z)
    
    rho1_mu = np.array([])
    rho2_mu = np.array([])
    mu_xiao_12 = np.array([])
        
    for X in rho1star:
        for Y in rho2star:
            rho1_mu = np.insert(rho1_mu,len(rho1_mu),X)
            rho2_mu = np.insert(rho2_mu,len(rho2_mu),Y)
            n_x = X*27.0
            n_y = Y*189.0
            
            deltax = 1e-4*n_x
            nT = n_x+n_y
            
            A1=-1.2*0.010324*math.log((1/(math.gamma(nT-n_x+1)*math.gamma(n_x+1)))*((((3*0.0000000003405)**3)/(0.0000000000221356**3))**n_x)*(((6*0.0000000003405)**3-(3*0.0000000003405)**3)/(0.0000000000221356**3))**(nT-n_x))
            A2=-1.2*0.010324*math.log((1/(math.gamma(nT-(n_x+deltax)+1)*math.gamma((n_x+deltax)+1)))*((((3*0.0000000003405)**3)/(0.0000000000221356**3))**(n_x+deltax))*(((6*0.0000000003405)**3-(3*0.0000000003405)**3)/(0.0000000000221356**3))**(nT-(n_x+deltax)))

            mu_id = (A2-A1)/deltax
            # mu_id = -nT*(1.2*0.010324)*math.log((12*0.0000000003405)**3)+3*nT*(1.2*0.010324)*math.log(0.0000000000221356)+(1.2*0.010324)*(nT*math.log(nT)-nT)
            #mu_id = (1.3*0.010324)*math.log((nT/(6*0.0000000003405)**3)*0.0000000000221356**3)
            mu_T = 1.2*0.010324*MU(X,Y) + mu_id            
            #mu_T = 1.2*0.010324*MU(X,Y) + mu_id            
           
            mu_xiao_12 = np.insert(mu_xiao_12,len(mu_xiao_12),mu_T)
            
    
  
def surf_fit():
    global F,n1_mu,n2_mu,mu
    n1_mu = []
    n1_mu = np.array(n1_mu)
    n2_mu = []
    n2_mu = np.array(n2_mu)
    mu = []
    mu = np.array(mu)
    
    x = Symbol('x')
    
    p00 = 0.07806  
    p10 = 0.03209  
    p01 = -0.004618  
    p20 = 0.0002499  
    p11 = -0.0006212  
    p02 = 7.192e-05  
    p30 = 2.618e-05  
    p21 = 9.237e-07  
    p12 = 7.367e-07  
    p03 = -1.369e-07          

   
    #z = p00 + p10*x + p01*y + p20*x**2 + p11*x*y + p02*y**2 + p30*x**3 + p21*x**2*y + p12*x*y**2 + p03*y**3

    #F = lambdify((x, y),z, 'numpy')
    nT_min = nT.min()
    nT_max = nT.max()
    n1_min = n1.min()
    n1_max = n1.max()
   
    for c in range(nT_min,nT_max):
        z = p00 + p10*x + p01*(c-x) + p20*x**2 + p11*x*(c-x) + p02*(c-x)**2 + p30*x**3 + p21*x**2*(c-x) + p12*x*(c-x)**2 + p03*(c-x)**3
        zprime = z.diff(x)
        F = lambdify(x,z,'numpy')
        MU = lambdify(x,zprime,'numpy')
        y = c-x
        Y = lambdify(x,y,'numpy')
        ''' 
        for X in range(int(n1_min+1),int(n1_max)):
            n1_mu = np.insert(n1_mu,len(n1_mu),X)
            n2_mu = np.insert(n2_mu,len(n2_mu),Y(X))
            mu = np.insert(mu,len(mu),MU(X))
        '''  
        for X in n1:
            n1_mu = np.insert(n1_mu,len(n1_mu),X)
            n2_mu = np.insert(n2_mu,len(n2_mu),Y(X))
            mu = np.insert(mu,len(mu),MU(X))


    
    
            
def plot_init():
    global fig,ax
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    plt.hold(True)
    
def plot_surf(X,Y,Z,map_type):       
   
    surf = ax.plot_trisurf(X, Y, Z, linewidth=0, antialiased=False,cmap=map_type)
    ax.scatter(X, Y, Z,s=0.3,color = 'black')
    #plt.xlim(0.8,0)  
    ax.set_title('3 sigma system')
    ax.set_xlabel('rho1star')
    ax.set_ylabel('rho2star')
    ax.set_zlabel('Mu')
    
    fig.colorbar(surf)
    plt.show()
    
    
if __name__ == '__main__':
    read_data()
    plot_init()
    xiao_data_1_2()
    surf_fit()
    
    
    plot_surf(rho1_mu,rho2_mu,mu_xiao_12,cm.coolwarm)
    #plot_surf(rho1_mu,rho2_mu,mu_xiao_11,cm.rainbow)
    
   # plot_surf(n1_mu/27.0,n2_mu/189.0,mu,cm.rainbow)
