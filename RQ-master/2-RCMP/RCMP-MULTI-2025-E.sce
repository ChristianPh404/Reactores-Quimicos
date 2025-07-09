clear; clc; 
// 2024-07-04-A.sce

// RCMP
// 1) A => B
// 2) B => C

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = f(t,x)
    // Variables
    CA = x(1)
    CB = x(2)
    CC = x(3)
    T  = x(4)
    TJ = x(5)
    // Ecuación de Arrhenius
    k1 = k10*exp(-20/(R*T)) 
    k2 = k20*exp(-50/(R*T))
    // Velocidad de reacción
    r1 = k1*CA
    r2 = k2*CB
    // Balance de materia para A
    // d(V*CA)dt = F*CA0 - F*CA - r1*V
    dCAdt = F*(CA0-CA)/V - r1
    // Balance de materia para B
    // d(V*CB)dt = F*CB0 - F*CB + r1*V - r2*V
    dCBdt = F*(CB0-CB)/V + r1 - r2
    // Balance de materia para C
    // d(V*CC)dt = F*CC0 - F*CC + r2*V
    dCCdt = F*(CC0-CC)/V + r2
    // Calor transferido del reactor a la camisa
    Q = UA*(T-TJ)
    // Balance de energía
    // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - F*RHO*CP*T - H1*r1*V - H2*r2*V - Q
    dTdt  = F*(T0-T)/V - (H1*r1+H2*r2)/(RHO*CP) - Q/(V*RHO*CP)
    dTjdt= Fj*(TJ0-TJ)/vj+Q/(vj*RHOj*cpj)
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dCBdt
    dxdt(3) = dCCdt
    dxdt(4) = dTdt
    dxdt(5) = dTjdt
endfunction

// CONSTANTES
H1 = -80;  // kj/mol
H2 = -60;  // kj/mol
RHO = 1000; // kg/m3
CP = 4.186; // kj/(kg*K)
cpj = 4.186; // kj/(kg*K)
RHOj = 1000; // kg/m3
V = 5; //m3
F = 1; // m3/h
T0 = 300; // K
CA0 = 1000; CB0 = 0; CC0 = 0; // mol/m^3
UA = 400; // kj/(h*K)
k10 = 3e3; // h-1
k20 = 2.5e7; // h-1
Fj = 0.25 ; // m3/h
TJ0 = 285; // K
vj = 1; // m3
R = 8.314e-3; // kj/(mol*K)
// CONDICIONES INICIALES
CAini = 1000; CBini = 0; CCini = 0; TJini = 285; // mol/m^3
Tini = 300; // K
xini = [CAini; CBini; CCini; Tini; TJini];

// TIEMPO
tee = 200; tfin = 100; dt = 0.01; t = 0:dt:tfin; // h

// RESOLVER
x = ode(xini,0,t,f);
CA = x(1,:); CAf = CA($);
CB = x(2,:); CBf = CB($);
CC = x(3,:); CCf = CC($);
T  = x(4,:); Tf  = T($);
TJ = x(5,:); TJf = TJ($);
disp(' concentraciones en Estado estacionario:');
disp('CAf = %f, CBf = %f, CCf = %f, Tf = %f, TJf = %f', CAf, CBf, CCf, Tf, TJf);

// GRÁFICAS
scf(1); clf(1);
plot(t,CA,t,CB,t,CC);
xgrid; xlabel('t'); legend('CA','CB','CC',-2,%f);

scf(2); clf(2);
plot(t,T,);
xgrid; xlabel('t'); legend('T',-2,%f);

//B 
tol = 1e-6; // Tolerancia para el método de Newton
for length(CA)
    Caee = diff(CA);
    Cbee = diff(CB);
    Ccee = diff(CC);
    Tee = diff(T);
    Tjee = diff(TJ);
    if abs(Caee) < tol && abs(Cbee) < tol && abs(Ccee) < tol && abs(Tee) < tol && abs(Tjee) < tol
        disp('Estado estacionario alcanzado en t = %f', t(length(CA)));
        disp('Concentraciones: CA = %f, CB = %f, CC = %f', CAf, CBf, CCf);
        disp('Temperaturas: T = %f, TJ = %f', Tf, TJf);
        break;
    end
//! piden a)	Obtener las condiciones del estado estacionario 
//* b)	Obtener la dinámica del proceso si el reactor se pone en marcha en las mismas condiciones que la alimentación


