clear; clc; 
// 2021-04-23-A.sce

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt=f(t, x)
    // Variables
    CA = x(1)
    T = x(2)
    TJ = x(3)
    // Ecuación de Arrhenius
    k = k0*exp(-E/(R*T))
    // Velocidad de reacción
    r = k*CA
    // Calor transferido del reactor a la camisa
    Q = UA*(T-TJ)
    // Balance de materia para A
    // d(V*CA)dt = F*CA0 - F*CA - r*V
    dCAdt = F*(CA0-CA)/V - r
    // Balance de energía en el reactor
    // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - F*RHO*CP*T -H*r*V - Q
    dTdt = F*(T0-T)/V - H*r/(RHO*CP) - Q/(V*RHO*CP)
    // Balance de energía en la camisa
    // d(VJ*RHOJ*CPJ*TJ)dt = FJ*RHOJ*CPJ*TJ0 - FJ*RHOJ*CPJ*TJ + Q
    dTJdt = FJ*(TJ0-TJ)/VJ + Q/(VJ*RHOJ*CPJ)
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dTdt
    dxdt(3) = dTJdt
endfunction

// CONSTANTES
k0 = 1.5E12; // min-1
E = 8.5E4; // J/mol
R = 8.314 // J/(mol*K)
H = -3.1E5; // J/mol
V = 1000; // L
F = 1200; // L/min
CA0 = 2; // mol/L
T0 = 290; // K
VJ = 100; // L
UA = 2.5E6; // J/(min*K)
FJ = 600; // L/min
RHOJ = 1000; // g/L
CPJ = 4.18; // J/(g*K)
TJ0 = 280; //// K
RHO = 800; // g/L
CP = 3.5; // J/(g*K)

// CONDICIONES INICIALES
CAini = 0; // mol/L
Tini = 350; // K
TJini = TJ0; // K
xini = [CAini; Tini; TJini];

// TIEMPO
tfin = 10; dt = 0.01; t = 0:dt:tfin; // min

// RESOLVER
x = ode(xini,0,t,f);
CA = x(1,:); CAee = CA($)
T = x(2,:); Tee = T($)
TJ = x(3,:); TJee = TJ($)

// GRÁFICAS
scf(1); clf(1);
plot(t,CA,'b-');
xgrid; xtitle('A.sce','t','CA');

scf(2); clf(2);
plot(t,T,'r-',t,TJ,'g-');
xgrid; xtitle('A.sce','t','T(rojo), TJ(verde)');

scf(3); //clf(3);
plot(Tini,CAee,'bo');
xgrid; xtitle('A.sce','Tini','CAee');

scf(4); //clf(4);
plot(Tini,Tee,'ro',Tini,TJee,'go');
xgrid; xtitle('A.sce','Tini','Tee(rojo), TJee(verde)');

DT = T-TJ;
scf(5); clf(5);
plot(t,DT);
xgrid; xtitle('A.sce','t','DT');

// (b)
DTobj = 5; // K
indexDTobj = find(DT<DTobj,1);
tDTobj = t(indexDTobj)
TDTobj = T(indexDTobj)
TJDTobj = TJ(indexDTobj)
plot(tDTobj,DTobj,'ro');

// (c)
[DTmin,indexDTmin] = min(DT)
tDTmin = t(indexDTmin)
plot(tDTmin,DTmin,'ro');
