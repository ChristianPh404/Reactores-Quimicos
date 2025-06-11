clear; clc;
// RCFP-1d.sce
// A  <=> B
// Adiabático
// Estado estacionario
// 2 reactores con enfriamiento intermedio

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdtau = f(tau,x)
    // Variables diferenciales
    CA = x(1)
    CB = x(2)
    T  = x(3)
    // Ecuación de Arrhenius
    kd = kd0*exp(-Ed/(R*T))
    //
    ki = Ki0*exp(-Ei/(R*T))
    // Velocidad de reacción
    // r = rd - ri = kd*CA - ki*Cb 
    r = kd*CA - CB*ki
    // Balance de materia para A
    // RDMP: d(V*CA)dt = -r*V 
    dCAdtau = -r
    // Balance de materia para B
    // RDMP: d(V*Cb)dt = r*V 
    dCBdtau = r
    // Balance de energía
    // RDMP: d(V*RHO*CP*T)dt = -H*r*V
    dTdtau = -H*r/(RHO*CP) 
    // Derivadas
    dxdtau(1) = dCAdtau
    dxdtau(2) = dCBdtau
    dxdtau(3) = dTdtau
endfunction

// CONSTANTES
kd0 =3.23E13; // cal/mol
Ki0 = 1.04E18; //cal/mol
Ed = 4.5E4 ;// cal/mol
Ei = 6.2E4; // cal/mol
H = Ed-Ei; // cal/mol
R =  1.9872; // cal/(mol·K)
RHO = 500; // g/L
CP = 1.5; // cal/(g·K)
F = 0.5; // L/min
L1 = 50; // dm
D1 = 2 ; // dm
L2 = 50 ; // dm
D2 = 3; // dm

// *********
// REACTOR 1
// *********

// CONSTANTES

V1 = %pi/4*D1^2*L1 // L
TAU1 = V1/F // h

// ENTRADA
CA0 = 1; CB0 = 0; // mol/L
T0 = 600; // K
x01 = [CA0;CB0;T0]; // Entrada del reactor 1

// TIEMPO DE RESIDENCIA
N =600; tau1 = 0:TAU1/N:TAU1; // h
l1 = 0:L1/N:L1; // dm

// RESOLVER
x1 = ode(x01,0,tau1,f);
CA1 = x1(1,:); CA1s = CA1($);
CB1 = x1(2,:); CB1s = CB1($);
T1  = x1(3,:); T1s  = T1($);
XA1 = 1 - CA1/CA0; XA1s = XA1($);
disp('Conversion del reactor 1 es: ', XA1s, ' y la temperatura es: ', T1s, ' K')


// *********
// REACTOR 2
// *********

// CONSTANTES
V2 = %pi/4*D2^2*L2 // L
TAU2 = V2/F // h

// ENTRADA
x02 = [CA1s;CB1s;T0];  // Enfriamiento: T0 puesto que en ambos son 600k
// TIEMPO DE RESIDENCIA
N = 600; tau2 = 0:TAU2/N:TAU2; // h
l2 = 0:L2/N:L2; // dm

// RESOLVER
x2 = ode(x02,0,tau2,f);
CA2 = x2(1,:); CA2s = CA2($) 
CB2 = x2(2,:); CB2s = CB2($) 
T2  = x2(3,:); T2s  = T2($)
XA2 = 1 - CA2/CA0; XA2s = XA2($) //! se pone CA2/Ca0 por que es la conversion respecto a la entrada inicial para ver la conversion total del reactor y no de forma individual   
disp('Conversion del reactor 2 es: ', XA2s, ' y la temperatura es: ', T2s, ' K')


// GRÁFICAS
scf(1); clf(1); 
plot(l1,XA1,'m',l2,XA2,'b--');
xgrid; xlabel('l,dm'); legend('XA1','XA2',-2,%f); 
xtitle("Conversion frente a longitud del reactor");

scf(2); clf(2); 
plot(l1,T1,'r',l2,T2,'g--');
xgrid; xlabel('l,dm'); legend('T1','T2',-2,%f);
xtitle("Temperatura frente a longitud del reactor", "l , dm", "T , K")

scf(3); clf(3);
plot(T1,XA1,'r--',T2,XA2,'b--');
xgrid; xlabel('T'); legend('XA1','XA2',-2,%f); 
xtitle("Conversion frente a temperatura", "T , k", "XA")

