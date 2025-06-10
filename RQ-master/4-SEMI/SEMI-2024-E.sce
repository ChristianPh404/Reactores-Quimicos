clear; clc;
// 2024-07-04-B.sce

// SEMI
// 1) A => B
// 2) A => C
// 3) B => D

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = f(t,x)
    // Variables diferenciales
    V  = x(1)
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
    T  = VT/V
    // Ecuación de Arrhenius
    k1 = 18*exp(-1800/T); // min
    k2 = 983*exp(-3200/T); // min
    k3 = 2670*exp(-3100/T); // min
    // Velocidad de reacción
    r1 = k1*CA
    r2 = k2*CA
    r3 = k3*CB
    // Caudal de alimentación    
    F = Fini  - Fini*t/tfin
    // Balance de materia global
    // d(V*RHO)dt = F*RHO
    dVdt = F
     // Balance de materia para A
    dNAdt = -r1*V - r2*V
    // Balance de materia para B
    dNBdt = F*CB0 + r1*V - r3*V
    // Balance de materia para C
    dNCdt = r2*V
    // Balance de materia para D
    dNDdt = r3*V
    // Balance de energía
    // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - (H1*r1+H2*r2+H3*r3)*V
    dVTdt = F*T0 - (H1*r1+H2*r2+H3*r3)*V/(RHO*CP)
    // Derivadas
    dxdt(1) = dVdt
    dxdt(2) = dNAdt
    dxdt(3) = dNBdt
    dxdt(4) = dNCdt
    dxdt(5) = dNDdt
    dxdt(6) = dVTdt
endfunction

// CONSTANTES
Fini = 10; // L/min
CB0 = 1.18; // mol/L  (Tantear para NDfin = 450)
T0 = 300; // K
RHO = 1; // kg/L
CP = 4.18; // kJ/(kg*K)
H1 = -376; // kJ/mol
H2 = -334; // kJ/mol
H3 = -209; // kJ/mol

// CONDICIONES INICIALES
Vini = 1000; // L
NAini = 500; NBini = 0; NCini = 0; NDini = 0; // mol
Tini = 280; // K
xini = [Vini; NAini; NBini; NCini; NDini; Vini*Tini];

// TIEMPO
tfin = 30; dt = 0.01; t = 0:dt:tfin; // min

// RESOLVER
x = ode(xini,0,t,f);
V  = x(1,:); Vfin = V($)
NA = x(2,:); NAfin = NA($)
NB = x(3,:); NBfin = NB($)
NC = x(4,:); NCfin = NC($)
ND = x(5,:); NDfin = ND($)
VT = x(6,:); T = VT./V; Tfin = T($)

[NBmax,indexNBmax] = max(NB)
tNBmax = t(indexNBmax)

dTdt = diff(T)/dt;
dTdtobj = 2; // K/min
indexdTdt = find(dTdt > dTdtobj);
tdTdt = dt*length(indexdTdt)

// GRÁFICAS
scf(1); clf(1);
plot(t,V);
xgrid; xlabel('t'); legend('V',-2,%f);

scf(2); clf(2);
plot(t,NA,t,NB,t,NC,t,ND);
plot(tNBmax,NBmax,'go');
xgrid; xlabel('t'); legend('NA','NB','NC','ND',-2,%f);

scf(3); clf(3);
plot(t,T,'r');
plot(t(indexdTdt),T(indexdTdt),'r.');
xgrid; xlabel('t'); legend('T',-2,%f);

scf(4); clf(4);
plot(t(1:$-1),dTdt,'r');
plot(t(indexdTdt),dTdt(indexdTdt),'r.');
xgrid; xlabel('t'); legend('dTdt',-2,%f);
