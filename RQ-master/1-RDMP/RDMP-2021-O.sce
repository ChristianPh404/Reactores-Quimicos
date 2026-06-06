clear; clc; 
// A + B <=> C
// Adiabático

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt=f(t, x)
    // Variables diferenciales
    CA = x(1)
    CB = x(2)
    CC = x(3)
    T = x(4)
    // Ecuación de Arrhenius
    kd = kd0*exp(-E/(R*T))
    // Ecuación de Van't Hoff
    Keq = Keq0*exp(-H/(R*T))
    // Velocidad de reacción
    // r = rd - ri = kd*CA*CB - ki*CC = kd*CA*CB - kd*CC/Keq
    r = kd*(CA*CB - CC/Keq)
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
    // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - F*RHO*CP*T -H*r*V - Q
    dTdt = F*(T0-T)/V - H*r/(RHO*CP)
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dCBdt
    dxdt(3) = dCCdt
    dxdt(4) = dTdt
endfunction

// CONSTANTES
H = -1.4E5; // J/mol
RHO = 1.2; // Kg/L
CP = 3850; // J/(Kg*K)
kd0 = 2.1E8; // L/(mol*s)
E = 6.3E4; // J/mol
Keq0 = 1.8E-21; // L/mol
R = 8.314; // J/(mol*K)
V = 1000; // L
F = 10; // L/s
CA0 = 1; // mol/L
CB0 = 1.2; // mol/L
CC0 = 0; // mol/L
T0 = 300; // K

// CONDICIONES INICIALES
CAini = 1; CBini = 1.2; CCini = 0; // mol/L
Tini = 300; // K
xini = [CAini; CBini; CCini; Tini];

// TIEMPO
dt = 0.01; tfin = 1000; t = 0:dt:tfin; // s

tol = 1E-5; ng = 1;
function z=g(t, x)
    z = max(abs(f(t,x))) - tol
endfunction

// RESOLVER
[x,rd] = ode("root", xini, 0, t, f, ng, g);
tee = rd(1)
t = t(1:find(t>tee,1));

CA = x(1,:); CAee = CA($)
CB = x(2,:); CBee = CB($)
CC = x(3,:); CCee = CC($)
T = x(4,:); Tee = T($)

dCAdt = diff(CA)/dt;
dCBdt = diff(CB)/dt;
dCCdt = diff(CC)/dt;
dTdt = diff(T)/dt;

// GRÁFICAS
// a)
scf(1); clf(1);
plot(t,CA,t,CB,t,CC);
xgrid; xtitle('Ejercicio B','t','CA(azul), CB(verde), CC(rojo)');

scf(2); clf(2);
plot(t,T);
xgrid; xtitle('Ejercicio B','t','T');

// b)

scf(3); clf(3);
plot(t(1:$-1),abs(dCAdt),t(1:$-1),abs(dCBdt),t(1:$-1),abs(dCCdt),t(1:$-1),abs(dTdt));
plot(tee,tol,'ro');
xgrid; xtitle('Ejercicio B','t','|dCAdt| (azul), |dCBdt| (verde), |dCCdt| (rojo), |dTdt| (cian)');
a3 = gca;
a3.log_flags = "nln";

// c)
scf(4); clf(4);
plot(T,CC,'ro');
xgrid; xtitle('Ejercicio B','T','CCee');
