clear; clc; 
// 2023-07-13-A.sce

// (a) Estado estacionario

// SISTEMA DE ECUACIONES ALGEBRAICAS
function dxdt = f(x)
    // Variables diferenciales
    CA = x(1)
    CB = x(2)
    CC = x(3)
    T  = x(4)
    // Ecuación de Arrhenius
    kd = kd0*exp(-E/(R*T))
    // Ecuación de Van't Hoff
    Keq = Keq0*exp(-H/(R*T))
    // Velocidad de reacción: A + B <=> C
    // r = rd - ri = kd*CA*CB - ki*CC = kd*CA*CB - kd*CC/Keq
    r = kd*(CA*CB - CC/Keq)
    // Calor transferido del reactor a la camisa
    Q = UA*(T-TJ)
    // Balance de materia para A
    // d(V*CA)dt = F*CA0 - F*CA - r*V
    dCAdt = F*(CA0-CA)/V - r
    // Balance de materia para B
    // d(V*CB)dt = F*CB0 - F*CB - r*V
    dCBdt = F*(CB0-CB)/V - r
    // Balance de materia para C
    // d(V*CC)dt = F*CC0 - F*CC + r*V
    dCCdt = F*(CC0-CC)/V + r
    // Balance de energía
    // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - F*RHO*CP*T - H*r*V - Q
    dTdt = F*(T0-T)/V - H*r/(RHO*CP) - Q/(V*RHO*CP) 
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dCBdt
    dxdt(3) = dCCdt
    dxdt(4) = dTdt
endfunction

// CONSTANTES
kd0 = 960; // m3/(mol*h)
Keq0 = 1.1E-4; // m3/mol
E = 3.2E4; // J/mol
H = -1.7E5; // J/mol
R = 8.314; // J/(mol*K)
RHO = 1000; // kg/m3
CP = 4180; // J/(kg*K)

V = 1; // m3
F = 1; // m3/h
CA0 = 1000; CB0 = 1500; CC0 = 0; // mol/m3
T0 = 300; // K
TJ = 290; // K
UA = 9E5; // J/(K*h)

// SOLUCIÓN SUPUESTA
CAeeguess = 1000; CBeeguess = 1500; CCeeguess = 0; // mol/m3
Teeguess = 300; // K
xeeguess = [CAeeguess; CBeeguess; CCeeguess; Teeguess];

// RESOLVER
[xee,fxee,info] = fsolve(xeeguess,f)
CAee = xee(1)
CBee = xee(2)
CCee = xee(3)
 Tee = xee(4)

// ESTABILIDAD
J = numderivative(f,xee) // Jacobiano
lambda = spec(J)  // Valores propios
Estable = and(real(lambda) < 0)


// (c) Dinámica

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = g(t,x)
    dxdt = f(x)
endfunction

// CONDICIONES INICIALES
CAini = 1000; CBini = 1500; CCini = 0; // mol/m3
Tini = 300; // K
xini = [CAini; CBini; CCini; Tini];

// TIEMPO
dt = 0.01; tfin = 5; t = 0:dt:tfin; // h

// RESOLVER
x = ode(xini,0,t,g);
CA = x(1,:); CAfin = CA($)
CB = x(2,:); CBfin = CB($)
CC = x(3,:); CCfin = CC($)
T  = x(4,:);  Tfin  = T($)

indexCCCA = find(CC > CA,1);
tCCCA = t(indexCCCA)

indexCCCB = find(CC > CB,1);
tCCCB = t(indexCCCB)

[Tmax,indexTmax] = max(T)
tTmax = t(indexTmax)

// GRÁFICAS
scf(1); clf(1); 
plot(t,CA,t,CB,t,CC); 
plot(tCCCA,CC(indexCCCA),'ro'); 
plot(tCCCB,CC(indexCCCB),'ro'); 
xgrid; xlabel('t'); legend('CA','CB','CC',-2,%f);

scf(2); clf(2); 
plot(t,T,'r');
plot(tTmax,Tmax,'ro');
xgrid; xlabel('t'); legend('T',-2,%f);
