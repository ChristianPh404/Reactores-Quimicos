// RCMP-2.sce
// A + B <=> C
// Isotermo

// (a) Estado estacionario

// SISTEMA DE ECUACIONES ALGEBRAICAS
function dxdt = f(x)
    // Variables
    CA = x(1)
    CB = x(2)
    CC = x(3)
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
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dCBdt
    dxdt(3) = dCCdt
endfunction

// CONSTANTES
F = 10; // L/min
CA0 = 1; CB0 = 1.5; CC0 = 0; // mol/L
V = 150; // L
kd = 0.1; // L/(mol*min)
Keq = 10; // L/mol

// SOLUCIÓN SUPUESTA
CAeeguess = 1; CBeeguess = 1.5; CCeeguess = 0; // mol/L
xeeguess = [CAeeguess; CBeeguess; CCeeguess];

// RESOLVER
[xee,fxee,info] = fsolve(xeeguess,f)
CAee = xee(1);
CBee = xee(2);
CPee = xee(3);
disp("CAEe = " + string(CAee) + " mol/L")
disp("CBEe = " + string(CBee) + " mol/L")
disp("CPEe = " + string(CPee) + " mol/L")
disp(info)
//dinamica

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = g(t,x)
    dxdt = f(x)
endfunction
// CONDICIONES INICIALES
CAini = 0; CBini = 0; CCini = 0; // mol/L
xini = [CAini; CBini; CCini];
// TIEMPO
tfin = 150; dt = 0.5; t = 0:dt:tfin; // min

// RESOLVER
x = ode(xini,0,t,g);
CA = x(1,:)
CB = x(2,:)
CC = x(3,:)
// Graficar
scf(1); clf(1);
plot(t,CA,'r',t,CB,'g',t,CC,'b')
xlabel("Tiempo (min)")
ylabel("Concentración (mol/L)")
legend("CA","CB","CC")
title("Dinámica del reactor")
