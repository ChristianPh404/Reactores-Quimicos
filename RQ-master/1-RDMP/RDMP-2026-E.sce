clear; clc; 
// RDMP-2026-E.sce
// Reacciones en paralelo en reactor discontinuo con camisa
// A => B (deseada)
// A => C (no deseada)

// FUNCION SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = f(t,x)
    CA = x(1);
    CB = x(2);
    CC = x(3);
    T  = x(4);
    TJ = x(5);
    
    // Cinética
    k1 = 1E1 * exp(-2600/T); // min-1
    k2 = 5E4 * exp(-5500/T); // min-1
    
    r1 = k1 * CA;
    r2 = k2 * CA;
    
    // Operación de la camisa con salto en escalón
    if t <= 600 then
        FJ = 0.5; // L/min
        TJin = 280; // K 
    else
        FJ = 0.6; // L/min
        TJin = 285; // K
    end
    
    // Calor transferido
    Q = UA * (T - TJ);
    
    // Balances de materia
    dCAdt = -r1 - r2;
    dCBdt = r1;
    dCCdt = r2;
    
    // Balances de energía 
    // Asumiendo disolución acuosa, rhoCp = 4.18 kJ/(L*K)
    dTdt = ((-H1)*r1 + (-H2)*r2) / rhoCp - Q / (V * rhoCp);
    
    // Balance en la camisa
    dTJdt = FJ * (TJin - TJ) / VJ + Q / (VJ * rhoCp);
    
    dxdt(1) = dCAdt;
    dxdt(2) = dCBdt;
    dxdt(3) = dCCdt;
    dxdt(4) = dTdt;
    dxdt(5) = dTJdt;
endfunction

// PARÁMETROS
H1 = -160; // kJ/mol
H2 = -80;  // kJ/mol
V = 100;   // L
VJ = 20;   // L
UA = 1;    // kJ/(min*K)
rhoCp = 4.18; // kJ/(L*K) (calor específico volumétrico del agua)

// CONDICIONES INICIALES 
CA0 = 1; // mol/L
CB0 = 0;
CC0 = 0;
TJ0 = 280; // K (asumimos que la camisa empieza a la temperatura de entrada)

// TIEMPO
tfin = 1000; // min
dt = 1;
t = 0:dt:tfin;

// BÚSQUEDA DE LA TEMPERATURA INICIAL ÓPTIMA
// Se iterará Tini para encontrar la que maximice la concentración final de B
Tini_rango = 270:0.5:320; // K
CB_finales = [];

for i = 1:length(Tini_rango)
    x0 = [CA0; CB0; CC0; Tini_rango(i); TJ0];
    sol_temp = ode(x0, 0, t, f);
    CB_finales(i) = sol_temp(2, $);
end

[max_CB, idx_opt] = max(CB_finales);
Tini_optima = Tini_rango(idx_opt);

disp("La temperatura inicial óptima es: " + string(Tini_optima) + " K")
disp("Concentración final de B alcanzada: " + string(max_CB) + " mol/L")

// RESOLVER CON TINI ÓPTIMA
x0_opt = [CA0; CB0; CC0; Tini_optima; TJ0];
x = ode(x0_opt, 0, t, f);

CA = x(1,:); CAfin = CA($);
CB = x(2,:); CBfin = CB($);
CC = x(3,:); CCfin = CC($);
T  = x(4,:); Tfin  = T($);
TJ = x(5,:); TJfin = TJ($);

// RESULTADOS FINALES
disp("--- RESULTADOS FINALES ---")
disp("CA final = " + string(CAfin) + " mol/L")
disp("CB final = " + string(CBfin) + " mol/L")
disp("CC final = " + string(CCfin) + " mol/L")
disp("T final reactor = " + string(Tfin) + " K")
disp("T final camisa = " + string(TJfin) + " K")

// --- APARTADOS ADICIONALES ---
// a) La concentración de A está comprendida entre las de B y C
index_A_mid = find((CA > CB & CA < CC) | (CA < CB & CA > CC));
if isempty(index_A_mid) then
    disp("a) La concentración de A nunca está comprendida entre B y C.")
else
    tiempo_a = length(index_A_mid) * dt;
    disp("a) El tiempo en el que CA está comprendida entre CB y CC es de: " + string(tiempo_a) + " min.")
end

// b) La temperatura del reactor está comprendida entre 298 y 300 K
index_T_mid = find(T > 298 & T < 300);
if isempty(index_T_mid) then
    disp("b) La temperatura del reactor nunca está entre 298 y 300 K.")
else
    tiempo_b = length(index_T_mid) * dt;
    disp("b) El tiempo en el que T del reactor está entre 298 y 300 K es de: " + string(tiempo_b) + " min.")
end

// c) La temperatura de la camisa está comprendida entre 286 y 288 K
index_TJ_mid = find(TJ > 286 & TJ < 288);
if isempty(index_TJ_mid) then
    disp("c) La temperatura de la camisa nunca está entre 286 y 288 K.")
else
    tiempo_c = length(index_TJ_mid) * dt;
    disp("c) El tiempo en el que TJ de la camisa está entre 286 y 288 K es de: " + string(tiempo_c) + " min.")
end

// GRÁFICAS
scf(1); clf(1);
if isempty(index_A_mid) then
    plot(t, CA, 'b-', t, CB, 'g-', t, CC, 'r-');
    legend("CA", "CB", "CC");
else
    plot(t, CA, 'b-', t, CB, 'g-', t, CC, 'r-', t(index_A_mid), CA(index_A_mid), 'k*');
    legend("CA", "CB", "CC", "a) CA en [CB, CC]");
end
xgrid; xlabel("Tiempo (min)"); ylabel("Concentración (mol/L)");

scf(2); clf(2);
if isempty(index_T_mid) & isempty(index_TJ_mid) then
    plot(t, T, 'r-', t, TJ, 'm-');
    legend("T Reactor", "T Camisa", 4);
elseif ~isempty(index_T_mid) & isempty(index_TJ_mid) then
    plot(t, T, 'r-', t, TJ, 'm-', t(index_T_mid), T(index_T_mid), 'k*');
    legend("T Reactor", "T Camisa", "b) T en [298, 300]", 4);
elseif isempty(index_T_mid) & ~isempty(index_TJ_mid) then
    plot(t, T, 'r-', t, TJ, 'm-', t(index_TJ_mid), TJ(index_TJ_mid), 'c*');
    legend("T Reactor", "T Camisa", "c) TJ en [286, 288]", 4);
else
    plot(t, T, 'r-', t, TJ, 'm-', t(index_T_mid), T(index_T_mid), 'k*', t(index_TJ_mid), TJ(index_TJ_mid), 'c*');
    legend("T Reactor", "T Camisa", "b) T en [298, 300]", "c) TJ en [286, 288]", 4);
end
xgrid; xlabel("Tiempo (min)"); ylabel("Temperatura (K)");

scf(3); clf(3);
plot(Tini_rango, CB_finales, 'b-', 'LineWidth', 2);
plot(Tini_optima, max_CB, 'ro', 'MarkerSize', 8);
xgrid; xlabel("Tini (K)"); ylabel("CB final (mol/L)");
legend("CB vs Tini", "Óptimo");
