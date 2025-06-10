clear; clc;
// 2024-06-10-B.sce 

// RCFP
// A => B
// No adiabático

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdtau = f(tau,x)
    // Variables diferenciales
    CA = x(1)
    T  = x(2)
    // Ecuación de Arrhenius
    k = 1.2E6*exp(-5000/T) // min-1
    // Velocidad de reacción
    r = k*CA
    // Balance de materia para A
    // RDMP: d(V*CA)dt = -r*V 
    dCAdtau = -r
    // Balance de energía
    // RDMP: d(V*RHO*CP*T)dt = -H*r*V - Q = -H*r*V - U*A*(T-TJ)
    // RDMP: dTdt = -H*r/(RHO*CP) - U*A*(T-TJ)/(V*RHO*CP)
    // A/V = %pi*D*L / (%pi/4*D^2*L) = 4/D 
    dTdtau = -H*r/(RHO*CP) - 4*U*(T-TJ)/(D*RHO*CP)
    // Derivadas
    dxdtau(1) = dCAdtau
    dxdtau(2) = dTdtau
endfunction

// CONSTANTES
L = 30; // m
D = 0.25; // m
V = %pi/4*D^2*L // m3
F = 0.50; // m3/min 
TAU = V/F // min
H = -4.1E5 // J/mol
RHO = 950; // kg/m3
CP = 3.5E3; // J/(kg*K)
TJ = 280; // K
U = 6.5E4; // J/(min*m2*K)

// ENTRADA
CA0 = 1500; // mol/m3
T0 = 310; // K
x0 = [CA0;T0];

// TIEMPO DE RESIDENCIA
tau = 0:TAU/1000:TAU; // min
l = 0:L/1000:L; // m

// RESOLVER
x = ode(x0,0,tau,f);
CA = x(1,:); CAs = CA($) 
T  = x(2,:); Ts  = T($)
XA = 1 - CA/CA0; XAs = XA($)

// GRÁFICAS
scf(1); clf(1); 
plot(l,XA);
xgrid; xlabel('l'); legend('XA',-2,%f);

d2XA = diff(XA,2);
indexXAi = find(d2XA(1:$-1).*d2XA(2:$)<0) + 1;
lXAi = l(indexXAi)
XAi = XA(indexXAi)
plot(lXAi,XAi,'o');

scf(2); clf(2); 
plot(l,T,'r');
xgrid; xlabel('l'); legend('T',-2,%f);

Tobj = 420; // K
indexTobj = find(T > Tobj);
lTobj = L/1000*length(indexTobj)
plot(l(indexTobj),T(indexTobj),'r.');
