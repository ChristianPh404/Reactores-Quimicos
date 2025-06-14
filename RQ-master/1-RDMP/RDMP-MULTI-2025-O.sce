clear; clc;
// 2024-06-10-A.sce

// RDMP-MULT
// 1)   A => B
// 2)   B => C   B> C
// No adiabático

// SISTEMA DE ECUACIONES ALGEBRAICAS
function dxdt = f(t,x)
    // Variables 
    CA = x(1);
    CB = x(2);
    CC = x(3);
    T  = x(4);
    TJ = x(5);

    // Ecuaciones de Arrhenius
    k1 = exp(-2500/T+8); // L/(h·mol)
    k2 = exp(-6000/T+17); // h-1

    // Velocidades de reacción
    r1 = k1*CA;
    r2 = k2*CB;

    // Balance de materia para A
    dCAdt = - r1;

    // Balance de materia para B
    dCBdt = r1 - r2;

    // Balance de materia para C
    dCCdt = r2;
    // Selección de condiciones de operación
    if t < tfin/3 then
        TJ0 = 300;
        Fj = 0.01;
    else
        TJ0 = 280;
        Fj = 0.05;
    end
    // Calor transferido del reactor a la camisa
    Q = UA*(T-TJ)

    // Balance de energía
    dTdt =  - (H1*r1 + H2*r2)/(RHO*CP) - Q/(V*RHO*CP);
    dTjdt= Fj*(TJ0-TJ)/vj+Q/(vj*RHOj*cpj)

    // Derivadas
    dxdt(1) = dCAdt;
    dxdt(2) = dCBdt;
    dxdt(3) = dCCdt;
    dxdt(4) = dTdt;
    dxdt(5) = dTjdt;
endfunction

// CONSTANTES
V    = 1; // m^3
H1   = -75; // kj/mol
H2   = -50; // Kj/mol
UA   = 200; // Kj/(h*K) 
vj   = 0.1; // m^3
CP   = 4.186; // Kj/(kg*K)
cpj  = 4.186;//kj/kg*k
RHO  = 1000; // kg/m^3
RHOj = 1000; // kg/m^3


// CONDICIONES INICIALES
CAini = 1000; // mol/m^3
CBini = 0; // mol/m^3
CCini = 0; // mol/m^3   
Tini  = 300; // K
TJini = 300; // K
xini  = [CAini;CBini;CCini;Tini;TJini];

// TIEMPO
tfin = 12; dt = 0.01; t = 0:dt:tfin; // h

// RESOLVER
x = ode(xini,0,t,f);
CA = x(1,:); CAfin = CA($);
CB = x(2,:); CBfin = CB($);
CC = x(3,:); CCfin = CC($);
T  = x(4,:); Tfin = T($);
TJ = x(5,:); TJfin = TJ($);
//! Apartado A
disp("CA final: " + string(CAfin) + " mol/m^3");
disp("CB final: " + string(CBfin) + " mol/m^3");
disp("CC final: " + string(CCfin) + " mol/m^3");
disp("T final: " + string(Tfin) + " K");
disp("TJ final: " + string(TJfin) + " K");

//! Apartado B 
[CBmax, indexCB] = max(CB);
disp("CB max: " + string(CBmax) + " mol/m^3");
// Graficar

scf(1); clf(1);
plot(t,CA,'r',t,CB,'g',t,CC,'b',t(indexCB),CB(indexCB),'k*'); // Graficar CA, CB, CC y el máximo de CB
xlabel('Tiempo (h)');
ylabel('Concentración (mol/m^3)');
legend('CA','CB','CC','CB máximo');
//!Apartado C
indices = find(T > 320); // Todos los índices donde T > 320
if indices <> [] then
    tiempo_T_superior_320 = length(indices) * dt;
    disp("Tiempo total en que T > 320 K: " + string(tiempo_T_superior_320) + " h");
end
//!Apartado D 
// seria de la siguiente manera

indextj = find(TJ > 310); // Todos los índices donde TJ > 310
if indextj <> [] then
    tiempo_T_superior_310 = length(indextj) * dt;
    disp("Tiempo total en que TJ > 310 K: " + string(tiempo_T_superior_310) + " h");
end
// Apartado E TJ máximo
indexTJmax = find(TJ == max(TJ)); // Índice donde TJ es máximo
disp("TJ máximo: " + string(TJ(indexTJmax)) + " K en el tiempo: " + string(t(indexTJmax)) + " h");
indexTmax = find(T == max(T)); // Índice donde T es máximo
disp("T máximo: " + string(T(indexTmax)) + " K en el tiempo: " + string(t(indexTmax)) + " h");
// Graficar Temp/tiempo
scf(2); clf(2);
plot(t, T, 'k', t, TJ, 'm',t(indices), T(indices), 'r.',t(indextj), TJ(indextj), 'b.',t(indexTJmax),TJ(indexTJmax),'r*',t(indexTmax),T(indexTmax),'b*'); 
xlabel('Tiempo (h)');
ylabel('Temperatura (K)');
title('Temperatura del reactor y de la camisa');
legend('T del reactor', 'T de la camisa', 'T > 320 K', 'TJ > 310 K', 'TJ máximo', 'T máximo');

