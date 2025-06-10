clear; clc; 
// 2023-06-20-A.sce

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = f(t,x)
    // Variables diferenciales
    CA = x(1)
    T  = x(2)
    // Ecuación de Arrhenius
    k = k0*exp(-E/(R*T))
    // Velocidad de reacción
    r = k*CA
    // Calor transferido del reactor a la camisa
    Q = UA*(T-TJ)
    // Balance de materia para A
    // d(V*CA)dt = -r*V
    dCAdt = -r                             
    // Balance de energía
    // d(V*RHO*CP*T)dt = -H*r*V - Q
    dTdt = -H*r/(RHO*CP) -  Q/(V*RHO*CP)   
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dTdt
endfunction

// CONSTANTES
kT1 = 0.0231; // min-1
T1 = 280; // K
kT2 = 0.822; // min-1
T2 = 350; // K
R = 8.314E-3; // kJ/(mol*K)
// kT1 = k0*exp(-E/(R*T1)) => log(kT1) = log(k0) - E/(R*T1)
// kT2 = k0*exp(-E/(R*T2)) => log(kT2) = log(k0) - E/(R*T2)
// log(kT1)-log(kT2) = -E/R*(1/T1-1/T2)
E = -R*(log(kT1)-log(kT2))/(1/T1-1/T2)
k0 = kT1/exp(-E/(R*T1))

RHO = 1; // kg/L
CP = 4.18; // kJ/(kg*K)
V = 2500; // L
TJ = 280; // K
UA = 300; // kJ/(min*K)

H = -393; // kJ/mol  (Tantear para obtener Tmax = 290)

// CONDICIONES INICIALES
CAini = 0.25; // mol/L
Tini = 280; // K
xini = [CAini; Tini];

// TIEMPO
tfin = 100; dt = 0.1; t = 0:dt:tfin; // min

// RESOLVER
x = ode(xini,0,t,f);
CA = x(1,:); CAfin = CA($)
T = x(2,:); Tfin = T($)

[Tmax,indexTmax] = max(T)
tTmax = t(indexTmax)

Ta = 284; Tb = 288;
indexTab = find(T > Ta & T < Tb);
tTab = dt*length(indexTab)

// GRÁFICAS
scf(1); clf(1); 
plot(t,CA);
xgrid; xlabel('t'); legend('CA',-2,%f);

scf(2); clf(2); 
plot(t,T,'r-',tTmax,Tmax,'ro');
plot(t(indexTab),T(indexTab),'r.');
xgrid; xlabel('t'); legend('T',-2,%f);

scf(3); // clf(3); 
plot(H,Tmax,'ro');
xgrid; xlabel('H'); ylabel('Tmax');  