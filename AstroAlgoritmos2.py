# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 22:46:11 2013

@author: usuario
"""
from numpy import *
from AstroAlgoritmos import *
def nutation(D,M,A):
    """
    Esta función calcula  la nutación del eje terrestre
    @author: Fernando Darío Mazzone
    """
    ddt=dt(A)
    JDE=gre2jul(D,M,A)+ddt/86400 
    T=(JDE-2451545)/36525 
    L=280.4665+36000.7698*T 
    Lp=218.3165+481267.8813*T 
    ome=125.04452-1934.136261*T 
    ome=ome*pi/180
    L=L*(pi/180)
    Lp=Lp*(pi/180)
    deps=(9.2*cos(ome)+.57*cos(2*L)+.1*cos(2*Lp)-.09*cos(2*ome))/3600
    epsilon0=23.43929111-46.815*T/3600-0.00059*T**2/3600+0.001813*T**3/3600
    dpsi=(-17.2*sin(ome)-1.32*sin(2*L)-0.23*sin(2*Lp)+0.21*sin(2*ome))/3600 
    eps=epsilon0+deps
    return dpsi,eps,T