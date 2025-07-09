clear; clc;
// 2024-06-10-B.sce 

// RCFP
// A + B => C
// No adiabático

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdtau = f(tau,x)
    // Variables diferenciales
    CA = x(1)
    CB = x(2)
    CC = x(3)
    T  = x(4)
    // Ecuación de Arrhenius
    k = 5E6*exp(-4800/T) //  L/(mol·h)
    // Velocidad de reacción
    r = k*CA*CB
    // Balance de materia para A
    // RDMP: d(V*CA)dt = -r*V 
    dCAdtau = -r
    dCBdtau = -r
    dCCdtau = +r
    // Balance de energía
    dTdtau = -H*r/(RHO*CP) - 4*U*(T-TJ)/(D*RHO*CP)
    // RDMP: d(V*RHO*CP*T)dt = -H*r*V - Q = -H*r*V - U*A*(T-TJ)
    // RDMP: dTdt = -H*r/(RHO*CP) - U*A*(T-TJ)/(V*RHO*CP)
    // A/V = %pi*D*L / (%pi/4*D^2*L) = 4/D 

    // Derivadas
    dxdtau(1) = dCAdtau
    dxdtau(2) = dCBdtau
    dxdtau(3) = dCCdtau
    dxdtau(4) = dTdtau
endfunction

// CONSTANTES
L = 20; // dm
D = 2; // dm
V = %pi/4*D^2*L // dm3
F = 10; // L/h
TAU = V/F // h
H = -4.3E4 // cal/mol
TJ = 280; // K
T = 300; // K
RHO = 1000 ; // g/dm3   
// ENTRADA
CA0 = 1; // mol/L
CB0 = 1; // mol/L
T0 = 300; // K
x0 = [CA0;CB0;T0];

// TIEMPO DE RESIDENCIA
tau = 0:TAU/1000:TAU; // min
l = 0:L/1000:L; // m

// ENTRADA
T20 = 280:0.01:300; // 
for i = 1:length(T20)
T2i = T20(i); // K
x02 = [CA0;CB0;CC0;T2i];  // Enfriamiento: T20
// RESOLVER
x2 = ode(x02,0,tau2,f);
T2  = x2(4,:); T2s(i)  = T2($)
end
Tobj= 290; // K
indexT2 = find(abs((T2s-Tobj)/Tobj) < 1e-5); // tolerancia de 0.01%
if isempty(indexT2)
    disp("No se encontró una temperatura de 600 K en el reactor 2.");
else
    indexT2 = indexT2(1); // Tomar el primer índice encontrado
    Tobt = T20(indexT2); // Temperatura objetivo
end
//RESOLVER PARA Tobt
x02 = [CA0;CB0;CC0;Tobt];  // Enfriamiento: T20
x2 = ode(x02,0,tau2,f);
// RESULTADOS
// RESULTADOS
CA = x2(1,:); CAs = CA($);
CB = x2(2,:); CB2s = CB2($);
CC = x2(3,:); CC2s = CC2($);
T  = x2(4,:); Ts  = T($);
XA = 1 - CA2/CA0; XAs = XA($);
C2result  = '  Las concentraciones de salida del reactor 2 son:' +ascii(10) + '  CA = ' + string(CA2s) + ' mol/L y CB = ' + string(CB2s) + ' mol/L' + ascii(10) + '  CC = ' + string(CC2s) + ' mol/L';
Xa2result = '  La conversión del reactor 2 es: ' + string(XAs) + ' y la temperatura es: ' + string(T2s) + ' K';




Indexobj = XA>= 0.5 & XA <= 0.7;
CAgraph = CA(Indexobj);
CBgraph = CB(Indexobj);
CCgraph = CC(Indexobj);
Tgraph = T(Indexobj);
lgraph = l(Indexobj);
// GRÁFICAS
scf(1); clf(1); 
plot(l,XA);
xgrid; xlabel('l'); legend('XA',-2,%f);
//
scf(2); clf(2);
plot(l,CAgraph,'r',l,CBgraph,'b--',l,CCgraph,'g-.');
xgrid; xlabel('l'); legend('CA','CB','CC',-2,%f);

scf(3); clf(3);
plot(T,CAgraph,'r--',T,CBgraph,'b--',T,CCgraph,'g-.');
xgrid; xlabel('T'); legend('CA','CB','CC',-2,%f);


// Maximos y minimos  de T

MaxT = max(T);
MinT = min(T);
findexT = find(T == MaxT);
disp("Índice de la temperatura máxima: " + string(findexT) + ", Temperatura máxima: " + string(MaxT) + " K");
findexT = find(T == MinT);
disp("Índice de la temperatura mínima: " + string(findexT) + ", Temperatura mínima: " + string(MinT) + " K");

//puntos de inflexion 

diffT = diff(T,2);
inf  = find(diffT(1:end-1) .* diffT(2:end) < 0) + 1; // +1 para ajustar el índice

disp("Puntos de inflexión en T: " + string(inf));