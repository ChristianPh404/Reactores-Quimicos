clear; clc;
// RCFP-E-2026.sce
// A + B => C
// 3 reactores adiabáticos en serie con enfriamiento intermedio

// FUNCION DE ECUACIONES DIFERENCIALES
function dxdtau = f(tau, x, H_val)
    CA = x(1);
    CB = x(2);
    CC = x(3);
    T  = x(4);
    
    // Cinética
    k = A_arr * exp(-E/(R*T)); // L/(mol*h)
    r = k * CA * CB;
    
    // Balances de materia
    dCAdtau = -r;
    dCBdtau = -r;
    dCCdtau = r;
    
    // Balance de energía adiabático
    // dTdtau = (-H * r) / (rhoCp)
    dTdtau = (-H_val) * r / rhoCp;
    
    dxdtau(1) = dCAdtau;
    dxdtau(2) = dCBdtau;
    dxdtau(3) = dCCdtau;
    dxdtau(4) = dTdtau;
endfunction

// PARÁMETROS GENERALES
A_arr = 1E10; // L/(mol*h)
E = 60000; // J/mol
R = 8.314; // J/(mol*K)
// Como es en medio acuoso, asumimos propiedades del agua:
rhoCp = 4180; // J/(L*K) 
F = 50; // L/h

// REACTOR 1
L1 = 20; // dm
D1 = 1.5; // dm
V1 = (%pi/4) * D1^2 * L1; // dm3 = L
TAU1 = V1 / F; // h

// REACTOR 2
L2 = 25; // dm
D2 = 2; // dm
V2 = (%pi/4) * D2^2 * L2; // L
TAU2 = V2 / F; // h

// REACTOR 3
L3 = 30; // dm
D3 = 2.5; // dm
V3 = (%pi/4) * D3^2 * L3; // L
TAU3 = V3 / F; // h

// CONDICIONES INICIALES R1
CA0 = 1; // mol/L
CB0 = 1; // mol/L
CC0 = 0; // mol/L
T0 = 300; // K

// VECTORES DE TIEMPO DE RESIDENCIA
N = 1000;
tau1 = 0:TAU1/N:TAU1;
tau2 = 0:TAU2/N:TAU2;
tau3 = 0:TAU3/N:TAU3;

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
disp("Iterando para hallar el calor de reacción (H) que cumple con XA = 80%...");
H_rango = -200000:100:-10000; // J/mol (búsqueda entre -200 y -10 kJ/mol)
XA3_finales = zeros(1, length(H_rango));

for i = 1:length(H_rango)
    H_test = H_rango(i);
    
    // R1
    x01 = [CA0; CB0; CC0; T0];
    x1 = ode(x01, 0, tau1, list(f, H_test));
    
    // R2 (Enfriamiento de 10 K)
    T2_in = x1(4, $) - 10;
    x02 = [x1(1, $); x1(2, $); x1(3, $); T2_in];
    x2 = ode(x02, 0, tau2, list(f, H_test));
    
    // R3 (Enfriamiento de 15 K)
    T3_in = x2(4, $) - 15;
    x03 = [x2(1, $); x2(2, $); x2(3, $); T3_in];
    x3 = ode(x03, 0, tau3, list(f, H_test));
    
    XA3_finales(i) = 1 - x3(1, $) / CA0;
end

[min_err, idx_opt] = min(abs(XA3_finales - 0.80));
H_opt = H_rango(idx_opt);

disp("La entalpía de reacción óptima estimada es H = " + string(H_opt) + " J/mol");

// -------------------------------------------------------------------------
// RESOLUCIÓN FINAL CON EL H ENCONTRADO
// -------------------------------------------------------------------------

// REACTOR 1
x01 = [CA0; CB0; CC0; T0];
x1 = ode(x01, 0, tau1, list(f, H_opt));
CA1_out = x1(1, $);
T1_out = x1(4, $);
XA1 = 1 - CA1_out / CA0;

// REACTOR 2
T2_in = T1_out - 10;
x02 = [x1(1, $); x1(2, $); x1(3, $); T2_in];
x2 = ode(x02, 0, tau2, list(f, H_opt));
CA2_out = x2(1, $);
T2_out = x2(4, $);
XA2 = 1 - CA2_out / CA0;

// REACTOR 3
T3_in = T2_out - 15;
x03 = [x2(1, $); x2(2, $); x2(3, $); T3_in];
x3 = ode(x03, 0, tau3, list(f, H_opt));
CA3_out = x3(1, $);
T3_out = x3(4, $);
XA3 = 1 - CA3_out / CA0;

// a) Mostrar resultados en consola
disp("--- RESULTADOS (ESTADO ESTACIONARIO) ---");
disp("REACTOR 1:");
disp("Conversión salida R1 = " + string(XA1*100) + " %");
disp("Temperatura salida R1 = " + string(T1_out) + " K");
disp("-------------------------");
disp("REACTOR 2:");
disp("Conversión salida R2 = " + string(XA2*100) + " %");
disp("Temperatura salida R2 = " + string(T2_out) + " K");
disp("-------------------------");
disp("REACTOR 3:");
disp("Conversión salida R3 = " + string(XA3*100) + " %");
disp("Temperatura salida R3 = " + string(T3_out) + " K");

// GRÁFICAS
// Calculamos la coordenada espacial z (longitud) para cada reactor
z1 = tau1 * F / ((%pi/4)*D1^2);
z2 = L1 + tau2 * F / ((%pi/4)*D2^2);
z3 = L1 + L2 + tau3 * F / ((%pi/4)*D3^2);

scf(1); clf(1);
plot(z1, 1 - x1(1,:)/CA0, 'r-', z2, 1 - x2(1,:)/CA0, 'g-', z3, 1 - x3(1,:)/CA0, 'b-');
xgrid; xlabel("Longitud del reactor z (dm)"); ylabel("Conversión X_A");
legend("Reactor 1", "Reactor 2", "Reactor 3", 4);
title("Evolución de la conversión a lo largo de los reactores");

scf(2); clf(2);
plot(z1, x1(4,:), 'r-', z2, x2(4,:), 'g-', z3, x3(4,:), 'b-');
plot([L1, L1], [x1(4,$), x1(4,$)-10], 'k--');
plot([L1+L2, L1+L2], [x2(4,$), x2(4,$)-15], 'k--');
xgrid; xlabel("Longitud del reactor z (dm)"); ylabel("Temperatura (K)");
legend("Reactor 1", "Reactor 2", "Reactor 3", 4);
title("Evolución de la temperatura con enfriamientos intermedios");

scf(3); clf(3);
plot(H_rango, XA3_finales, 'b-', 'LineWidth', 2);
plot(H_opt, XA3_finales(idx_opt), 'ro', 'MarkerSize', 8);
xgrid; xlabel("Entalpía supuesta H (J/mol)"); ylabel("Conversión final R3");
legend("XA vs H", "Óptimo para 80%");
title("Búsqueda del calor de reacción");
