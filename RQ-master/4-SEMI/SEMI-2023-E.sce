clear; clc;
// 2023-07-13-B.sce

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = f(t,x)
    // Variables diferenciales
    V  = x(1)
    NA = x(2)
    NB = x(3)
    NC = x(4)
    VT = x(5)
    // Concentraciones y temperatura
    CA = NA/V
    CB = NB/V
    CC = NC/V
    T  = VT/V
    // Ecuación de Arrhenius
    k = k0*exp(-E/(R*T))
    // Velocidad de reacción: A + B => C 
    r = k*CA*CB    
    // Caudal de alimentación    
    if t < tc then 
        F = F1; CB0 = CB01; T0 = T01;  // Etapa 1
    else 
        F = F2; CB0 = CB02; T0 = T02;  // Etapa 2
    end
    // Balance de materia global
    // d(V*RHO)dt = F*RHO
    dVdt = F
    // Balance de materia para A
    dNAdt = -r*V
    // Balance de materia para B
    dNBdt = F*CB0 - r*V
    // Balance de materia para C
    dNCdt = r*V
    // Balance de energía
    // d(V*RHO*CP*T)dt = F*RHO*CP*TO - H*r*V
    dVTdt = F*T0 - H*r*V/(RHO*CP)
    // Derivadas
    dxdt(1) = dVdt
    dxdt(2) = dNAdt
    dxdt(3) = dNBdt
    dxdt(4) = dNCdt
    dxdt(5) = dVTdt
endfunction

// CONSTANTES
kTa = 0.578; // L/(mol*h))
Ta = 300; // K
kTb = 6.25; // L/(mol*h))
Tb = 350; // K
R = 1.987; // cal/(mol*K)
E = -R*(log(kTa)-log(kTb))/(1/Ta-1/Tb)
k0 = kTa/exp(-E/(R*Ta))

F1 = 50; F2 = 20; // L/h
CB01 = 1; CB02 = 0.5; // mol/L
T01 = 300; T02 = 278; // K

RHO = 1000; // g/L
CP = 1; // cal/(g*K)
H = -8E4; // cal/mol

// CONDICIONES INICIALES
Vini = 500; // L
CAini = 1; // mol/L
NAini = Vini*CAini; NBini = 0; NCini = 0; // mol
Tini = 300; // K
xini = [Vini; NAini; NBini; NCini; Vini*Tini];

// TIEMPO
tc = 9.07; // Tantear para Tc > 330
tfin = 23.96; // Tantear para NAfin < 10
dt = 0.01; t = 0:dt:tfin; // h

// RESOLVER
x = ode(xini,0,t,f);
V  = x(1,:); Vc = V(t==tc), Vfin = V($)
NA = x(2,:); NAc = NA(t==tc), NAfin = NA($)
NB = x(3,:); NBc = NB(t==tc), NBfin = NB($)
NC = x(4,:); NCc = NC(t==tc), NCfin = NC($)
VT = x(5,:); T = VT./V; Tc = T(t==tc), Tfin = T($)

// GRÁFICAS
scf(1); clf(1);
plot(t,V);
plot(tc,Vc,'o');
plot(tfin,Vfin,'.');
xgrid; xlabel('t'); legend('V',-2,%f);

scf(2); clf(2);
plot(t,NA,t,NB,t,NC);
plot(tc,NAc,'o',tc,NBc,'o',tc,NCc,'o');
plot(tfin,NAfin,'.',tfin,NBfin,'.',tfin,NCfin,'.');
xgrid; xlabel('t'); legend('NA','NB','NC',-2,%f);

scf(3); clf(3);
plot(t,T,'r');
plot(tc,Tc,'ro');
plot(tfin,Tfin,'r.');
xgrid; xlabel('t'); legend('T',-2,%f);
