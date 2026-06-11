clear; clc; 
// RDMP-E-2020.sce

function dxdt = f(t, x)
    CA = x(1);
    CB = x(2);
    CC = x(3);
    T = x(4);
    TJ = x(5);
    
    // Cinética
    k1 = A1 * exp(-Ea1/(R*T));
    k2 = A2 * exp(-Ea2/(R*T));
    
    r1 = k1 * CA;
    r2 = k2 * CB;
    
    // Balances de materia
    dCAdt = -r1;
    dCBdt = r1 - r2;
    dCCdt = r2;
    
    // Calor transferido
    Q = UA * (T - TJ);
    
    // Balance de energía del reactor
    // calor generado por unidad de volumen
    q_gen = (-H1)*r1 + (-H2)*r2; 
    dTdt = ( q_gen - Q/V ) / rhoCp;
    
    // Balance de energía de la camisa
    dTJdt = ( Fc_val * rhoCp_c * (Tce - TJ) + Q ) / (Vc * rhoCp_c);
    
    dxdt(1) = dCAdt;
    dxdt(2) = dCBdt;
    dxdt(3) = dCCdt;
    dxdt(4) = dTdt;
    dxdt(5) = dTJdt;
endfunction

// Constantes
A1 = 6.72; // min-1
Ea1 = 2980; // cal/mol
H1 = -1e4; // cal/mol

A2 = 367; // min-1
Ea2 = 5960; // cal/mol
H2 = -2e4; // cal/mol

V = 1000; // L
UA = 1e4; // cal/(min K)
Vc = 100; // L

// Asumimos propiedades del agua
rhoCp = 1000; // cal/(L K) (densidad 1 kg/L, Cp 1000 cal/(kg K))
rhoCp_c = 1000; // cal/(L K)

Tce = 280; // K
R = 1.987; // cal/(mol K)

// Condiciones Iniciales
CA0 = 1; // mol/L
CB0 = 0;
CC0 = 0;
T0 = 300; // K
TJ0 = 280; // K
x0 = [CA0; CB0; CC0; T0; TJ0];

tfin = 120; // min
dt = 0.5;
t = 0:dt:tfin;

Fc_val = 25; // L/min (Caudal nominal)

// --- Apartado A ---
x = ode(x0,0,t,f);

CA_nom = x(1,:);
CB_nom = x(2,:);
CC_nom = x(3,:);
T_nom  = x(4,:);
TJ_nom = x(5,:);
[Tmax,IndTmax] = max(T_nom);
tTmax = t(IndTmax);

disp("=== Apartado A ===");
disp("Tmax = " + string(Tmax) + " K en t = " + string(tTmax) + " min");

// --- Apartado B ---
disp("=== Apartado B: Cruces de concentraciones ===");
[min_AB, idx_AB] = min(abs(CA_nom - CB_nom));
t_cruce_AB = t(idx_AB);
y_AB = CA_nom(idx_AB);
disp("Cruce CA = CB aprox en t = " + string(t_cruce_AB) + " min");

[min_BC, idx_BC] = min(abs(CB_nom - CC_nom));
t_cruce_BC = t(idx_BC);
y_BC = CB_nom(idx_BC);
disp("Cruce CB = CC aprox en t = " + string(t_cruce_BC) + " min");

[min_AC, idx_AC] = min(abs(CA_nom - CC_nom));
t_cruce_AC = t(idx_AC);
y_AC = CA_nom(idx_AC);
disp("Cruce CA = CC aprox en t = " + string(t_cruce_AC) + " min");


// --- Apartado C ---
disp("=== Apartado C: Máxima transmisión de calor ===");
Q_nom = UA * (T_nom - TJ_nom);
[Q_max, idx_Qmax] = max(Q_nom);
t_Qmax = t(idx_Qmax);
disp("Máxima transmisión de calor Q = " + string(Q_max) + " cal/min en t = " + string(t_Qmax) + " min");

// --- Apartado D ---
disp("=== Apartado D: Intervalo de caudal ===");
Fcs = 0:0.1:200; // barrido fino de caudales
Fcs_validos = [];

vec_CC_end = zeros(1, length(Fcs));
vec_T_max = zeros(1, length(Fcs));

for i = 1:length(Fcs)
    Fc_val = Fcs(i); // Cambiamos la variable global
    sol = ode(x0, 0, t, f);
    
    CC_end_val = sol(3,$);
    T_max_val = max(sol(4,:));
    
    vec_CC_end(i) = CC_end_val;
    vec_T_max(i) = T_max_val;
    
    if (CC_end_val > 0.85) & (T_max_val < 310) then
        Fcs_validos = [Fcs_validos, Fc_val];
    end
end

if isempty(Fcs_validos) then
    disp("No hay ningún caudal que cumpla las condiciones.");
else
    disp("Rango de caudales válido: [" + string(min(Fcs_validos)) + " , " + string(max(Fcs_validos)) + "] L/min");
end

// Restauramos el caudal nominal
Fc_val = 25;

// Gráficas
scf(1); clf(1);
plot(t, CA_nom, 'b-', t, CB_nom, 'm-', t, CC_nom, 'g-');
plot(t_cruce_AB, y_AB, 'ko');
plot(t_cruce_BC, y_BC, 'ko');
plot(t_cruce_AC, y_AC, 'ko');
xgrid; xlabel("Tiempo (min)"); ylabel("Concentración (mol/L)");
legend("CA", "CB", "CC", "Cruces");

scf(2); clf(2);
plot(t, T_nom, 'r-', t, TJ_nom, 'b--');
plot(t_Qmax, T_nom(idx_Qmax), 'ro');
xgrid; xlabel("Tiempo (min)"); ylabel("Temperatura (K)");
legend("T Reactor", "T Camisa", "Punto Q_{max}", 4);

scf(3); clf(3);
subplot(2,1,1);
plot(Fcs, vec_CC_end, 'r-');
plot(Fcs, 0.85*ones(1,length(Fcs)), 'g--');
xgrid; xlabel("Caudal (L/min)"); ylabel("C_C (mol/L)");
legend("Concentración final", "Mínimo exigido");

subplot(2,1,2);
plot(Fcs, vec_T_max, 'b-');
plot(Fcs, 310*ones(1,length(Fcs)), 'g--');
xgrid; xlabel("Caudal (L/min)"); ylabel("Temperatura (K)");
legend("Temperatura máxima", "Máximo permitido");