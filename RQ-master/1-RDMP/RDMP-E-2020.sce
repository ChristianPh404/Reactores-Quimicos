clear; clc; 
// RDMP-E-2020.sce

function t_cruce = encontrar_cruces(C1, C2, t)
    t_cruce = [];
    diffC = C1 - C2;
    for i = 1:(length(t)-1)
        if diffC(i)*diffC(i+1) < 0 then
            // interpolación lineal cruzando cero
            t_int = t(i) - diffC(i) * (t(i+1) - t(i)) / (diffC(i+1) - diffC(i));
            t_cruce = [t_cruce, t_int];
        elseif diffC(i) == 0 then
            // Evitar duplicados (ya que el if capturaría el i, el sgt i+1 no dispara porque <0 no se cumple y ==0 no ocurrirá salvo que se mantenga)
            if i == 1 | diffC(i-1) <> 0 then
                t_cruce = [t_cruce, t(i)];
            end
        end
    end
endfunction

function dxdt = f(t, x, Fc_val)
    CA = x(1);
    CB = x(2);
    CC = x(3);
    T = x(4);
    TJ = x(5);
    
    R = 1.987; // cal/(mol K)
    
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

Fc_nominal = 25; // L/min

// --- Apartado A ---
x_nom = ode(x0, 0, t, list(f, Fc_nominal));

CA_nom = x_nom(1,:);
CB_nom = x_nom(2,:);
CC_nom = x_nom(3,:);
T_nom  = x_nom(4,:);
TJ_nom = x_nom(5,:);

// --- Apartado B ---
disp("=== Apartado B: Cruces de concentraciones ===");
cruce_AB = encontrar_cruces(CA_nom, CB_nom, t);
disp("Cruce CA = CB aprox en t = " + string(cruce_AB) + " min");
y_AB = interp1(t, CA_nom, cruce_AB, 'linear');

cruce_BC = encontrar_cruces(CB_nom, CC_nom, t);
disp("Cruce CB = CC aprox en t = " + string(cruce_BC) + " min");
y_BC = interp1(t, CB_nom, cruce_BC, 'linear');

cruce_AC = encontrar_cruces(CA_nom, CC_nom, t);
disp("Cruce CA = CC aprox en t = " + string(cruce_AC) + " min");
y_AC = interp1(t, CA_nom, cruce_AC, 'linear');

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
    Fc_test = Fcs(i);
    sol = ode(x0, 0, t, list(f, Fc_test));
    
    CC_end_val = sol(3,$);
    T_max_val = max(sol(4,:));
    
    vec_CC_end(i) = CC_end_val;
    vec_T_max(i) = T_max_val;
    
    if (CC_end_val > 0.85) & (T_max_val < 310) then
        Fcs_validos = [Fcs_validos, Fc_test];
    end
end

if isempty(Fcs_validos) then
    disp("No hay ningún caudal que cumpla las condiciones.");
else
    disp("Rango de caudales válido: [" + string(min(Fcs_validos)) + " , " + string(max(Fcs_validos)) + "] L/min");
end

// Gráficas
scf(1); clf(1);
plot(t, CA_nom, 'b-', t, CB_nom, 'm-', t, CC_nom, 'g-');
plot(cruce_AB, y_AB, 'ko', 'MarkerSize', 5);
plot(cruce_BC, y_BC, 'ko', 'MarkerSize', 5);
plot(cruce_AC, y_AC, 'ko', 'MarkerSize', 5);
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