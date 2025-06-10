clear; clc;
// 2024-06-10-A.sce

// RDMP-MULT
// 1) 2 A => B
// 2)   A => C 
// No adiabático

// SISTEMA DE ECUACIONES ALGEBRAICAS
function dxdt = f(t,x)
    // Variables 
    CA = x(1)
    T  = x(2)
    // Ecuaciones de Arrhenius
    k1 = 1.2E16*exp(-14500/T) // L/(h·mol)
    k2 = 9.4E12*exp(-10600/T) // h-1
    // Velocidades de reacción
    r1 = k1*CA^2
    r2 = k2*CA
    // Balance de materia para A
    // d(V*CA)dt = - r1*V - r2*V
    dCAdt = - r1 - r2
    // Calor transferido del reactor a la camisa
    Q = UA*(T-TJ)
    // Balance de energía
    // d(V*RHO*CP*T)dt = - H1*r1*V - H2*r2*V - Q
    dTdt =  - (H1*r1 + H2*r2)/(RHO*CP) - Q/(V*RHO*CP)
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dTdt
endfunction

// CONSTANTES
H1  = -29000; // cal/mol
H2  = -23000; // cal/mol
CP  =  1; // cal/(g*K)
RHO =  1000; // g/L
V   =  0.1 // L
TJ  =  320; // K
UA  =  74.5; // cal/(h*K)  (Tantear para Tmax = 335 K)

// CONDICIONES INICIALES
CAini = 5; // mol/L
Tini  = 320; // K
xini  = [CAini;Tini];

// TIEMPO
tfin = 50; dt = 0.01; t = 0:dt:tfin; // h

// RESOLVER
x = ode(xini,0,t,f);
CA = x(1,:); CAfin = CA($)
T  = x(2,:); Tfin = T($)

CAobj = 1;
indexCAobj = find(CA<CAobj,1);
tCAobj = t(indexCAobj)

[Tmax,indexTmax] = max(T)
tTmax = t(indexTmax)

scf(1); clf(1); 
plot(t,CA,tCAobj,CAobj,'bo');
xgrid; xlabel('t'); legend('CA',-2,%f);

scf(2); clf(2); 
plot(t,T,'r',tTmax,Tmax,'ro');
xgrid; xlabel('t'); legend('T',-2,%f);

scf(3); 
plot(UA,Tmax,'ro');
xgrid; xlabel('UA'); ylabel('Tmax');
