from numpy import *
def gre2jul(D,M,A):
    """
    Esta función calcula el tiempo en días julianos para una determinada
    fecha. Es valida a partir del cambio al calendario gregoriano, que fue
    hace mucho. Igual esto no importa mucho. El tiempo es medido en tiempo
    universal. D es medido en días, con decimales si es necesario. M es el número 
    del mes y A el año. Así, el 23 de febrero de 1997 a las 12hs corresponde a
    D=23.5,   M=2, A=1997. El resultadeo es JD indicando el tiempo en días 
    julianos. 
    @author: Fernando Darío Mazzone
    """
    I=M<3
    A[I]=A[I]-1.0;
    M[I]=M[I]+12.0;
    C=floor(A/100.0);
    B=2-C+floor(C/4.0);
    JD=floor(365.25*(A+4716))+floor(30.6001*(M+1))+D-1524.5+B;
    return JD
def jul2gre(J):
    """
    Dado J en días julianos retorna año, mes y día en el calendario gregoriano
    @author: Fernando Darío Mazzone
    """
    J=J+0.5
    Z=floor(J)
    F=J-Z
    alpha=floor((Z-1867216.25)/36524.25)
    A=Z+1+alpha-floor(alpha/4.0)
    Ind=Z<2299161
    A[Ind]=Z[Ind]
    B=A+1524
    C=floor((B-122.1)/365.25)
    D=floor(365.25*C)
    E=floor((B-D)/30.6001)
    Dia=B-D-floor(30.6001*E)+F
    Ind1=E<14
    Mes=zeros(len(J))
    Mes[Ind1]=E[Ind1]-1
    Ind2=E>=14
    Mes[Ind2]=E[Ind2]-13
    Anio=C-4715
    Ind=Mes>2
    Anio[Ind]=C[Ind]-4716
    return Anio,Mes,Dia
def lst(JD,longitud):
    """
    Dado JD en días julianos y longitud en grados retorna la hora sideral local
    @author: Fernando Darío Mazzone
    """   
    if not longitud==0:
        longitud=-longitud+floor(180/longitud)*360
    T=(JD-2451545)/36525
    LSTG=280.46061837+360.98564736629*(JD-2451545)+0.000387933*T**2-(T**3/38710000)
    LST1=(LSTG-floor(LSTG/360)*360-longitud)
    LST=LST1-floor(LST1/360)*360
    return LST
def dt(a):
    """
    Esta función calcula dt-UT la diferencia entre el tiempo dinámico y
    universal, el argumento es el  dia juliano
    @author: Fernando Darío Mazzone
    """
    dt=-15+0.00325*(a-1810)**2
    return dt
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