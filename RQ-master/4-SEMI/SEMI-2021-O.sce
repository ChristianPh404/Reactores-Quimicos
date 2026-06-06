clear; clc; 
// 2021-06-04-B.sce
// 1) A + B => C
// 2) A + B => D

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt=f(t, x)
    // Variables diferenciales
    V = x(1)
    NA = x(2)
    NB = x(3)
    NC = x(4)
    ND = x(5)
    VT = x(6)
    // Concentraciones y temperatura
    CA = NA/V
    CB = NB/V
    CC = NC/V
    CD = ND/V
    T = VT/V
    // Ecuación de Arrhenius
    k1 = k10*exp(-E1/(R*T))
    k2 = k20*exp(-E2/(R*T))
    // Velocidad de reacción
    r1 = k1*CA*CB
    r2 = k2*CA*CB
    // Caudal de alimentación
    // F = Fini - Fini/tfin*t
    F = Fini*(1-t/tfin)
    // Balance de materia global
    // d(V*RHO)dt = F*RHO
    dVdt = F
    // Balance de materia para A
    dNAdt = -(r1+r2)*V
    // Balance de materia para B
    dNBdt = F*CB0 - (r1+r2)*V
    // Balance de materia para C
    dNCdt = r1*V
    // Balance de materia para D
    dNDdt = r2*V
    // Balance de energía
    // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - (H1*r1+H2*r2)*V
    dVTdt = F*T0 - (H1*r1+H2*r2)*V/(RHO*CP)
    // Derivadas
    dxdt(1) = dVdt
    dxdt(2) = dNAdt
    dxdt(3) = dNBdt
    dxdt(4) = dNCdt
    dxdt(5) = dNDdt
    dxdt(6) = dVTdt
endfunction

// CONSTANTES
k10 = 1.2E7; // L/(mol*h)
k20 = 1.5E7; // L/(mol*h)
E1 = 10200; // cal/mol
E2 = 9800; // cal/mol
H1 = -5E4; // cal/mol
H2 = -2E4; // cal/mol
R = 1.987; // cal/(mol*K)
RHO = 1000; // g/L
CP = 1; // cal/(g*K)

Fini = 40; // L/h
CB0 = 1; // mol/L
T0 = 290; // K

// CONDICIONES INICIALES
Vini = 500; // L
CAini = 1; // mol/L
NAini = Vini*CAini; NBini = 0; NCini = 0; NDini = 0; // mol
Tini = 335; // K
xini = [Vini; NAini; NBini; NCini; NDini; Vini*Tini];

// TIEMPO
tfin = 25; dt = 0.01; t = 0:dt:tfin; // h

// RESOLVER
x = ode(xini,0,t,f);

// (b)
V = x(1,:); Vfin = V($)

scf(1); clf(1);
plot(t,V);
xgrid; xtitle('B.sce','t','V');

Vobj = 800; // L
indexV = find(V>Vobj,1);
tV = t(indexV)
plot(tV,Vobj,'ro');

// (c)
NA = x(2,:); NAfin = NA($)
NB = x(3,:); NBfin = NB($)
NC = x(4,:); NCfin = NC($)
ND = x(5,:); NDfin = ND($)

scf(2); clf(2);
plot(t,NA,t,NB,t,NC,t,ND);
xgrid; xtitle('B.sce','t','NA(azul), NB(verde), NC(rojo), ND(cian)');

indexNA = find(NA>NC & NA<ND);
tNA = dt*length(indexNA)
plot(t(indexNA),NA(indexNA),'b.');

// (d)
VT = x(6,:); T = VT./V; Tfin = T($)

scf(3); clf(3);
plot(t,T);
xgrid; xtitle('B.sce','t','T');

dT = diff(T);
indexdT = find(dT<0);
tdT = dt*length(indexdT)
plot(t(indexdT),T(indexdT),'b.');
