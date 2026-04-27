clear; clc; 
// RDMP-O-Unknown.sce

function dxdt = f(t, x, TJ_val)
    CA = x(1);
    CB = x(2);
    CC = x(3);
    T = x(4);  
    // Cinética
    kd = A_d * exp(-Ea/(R*T));
    Keq = A_vh * exp(-H/(R*T)); // H es negativo
    
    r = kd * (CA * CB - CC / Keq);
    
    // Balances
    dCAdt = -r;
    dCBdt = -r;
    dCCdt = r;
    
    // Calor extraído
    Q = U * S * (T - TJ_val);
    
    // Balance de energía
    dTdt = ( (-H)*r*V - Q ) / (V * rho * Cp);
    
    dxdt(1) = dCAdt;
    dxdt(2) = dCBdt;
    dxdt(3) = dCCdt;
    dxdt(4) = dTdt;
endfunction

// Constantes
A_d = 1.8E8; // m3/(kmol h)
A_vh = 2.5E-22; // m3/kmol
Ea = 1.5E4; // kcal/kmol
H = -3.5E4; // kcal/kmol
R = 1.987; // kcal/(kmol K)
V = 1; // m3
S = 5; // m2
U = 1.5; // asumimos kcal/(h m2 K) - en el enunciado pone kg pero es un typo común

rho = 1200; // kg/m3
Cp = 0.9; // kcal/(kg K)

// Condiciones iniciales
CA0 = 1; // kmol/m3
CB0 = 1.5;
CC0 = 0;
T0 = 320; // K

TJ_ab = 300; // K

tfin = 120; // h
dt = 0.1;
t = 0:dt:tfin;
x0 = [CA0; CB0; CC0; T0];

// -- Apartados A y B --
x = ode(x0, 0, t, list(f, TJ_ab));

CA = x(1,:);
CB = x(2,:);
CC = x(3,:);
T = x(4,:);

disp("=== Apartado A ===");
disp("Concentraciones finales (t=120h):");
disp("CA = " + string(CA($)) + " kmol/m3");
disp("CB = " + string(CB($)) + " kmol/m3");
disp("CC = " + string(CC($)) + " kmol/m3");

disp("=== Apartado B ===");
index_330 = find(T > 330);
if isempty(index_330) then
    disp("La temperatura nunca supera los 330 K.");
else
    tiempo_mayor_330 = length(index_330) * dt;
    disp("El reactor está a más de 330 K por aproximadamente " + string(tiempo_mayor_330) + " horas.");
end
disp("=== Apartado C ===");

Tj_rango = 312:0.01:314; // Primero probe de 280 a 380 y salio 313
XA_finales = [];

for i = 1:length(Tj_rango)
    sol_temp = ode(x0, 0, t, list(f, Tj_rango(i)));
    CA_f = sol_temp(1, $);
    XA_f = (CA0 - CA_f) / CA0;
    XA_finales(i) = XA_f;
end

[max_XA, idx_opt] = max(XA_finales);
Tj_optima = Tj_rango(idx_opt);

disp("La temperatura de camisa óptima es: " + string(Tj_optima) + " K" + " que da una conversion de A de: " + string(max_XA*100) + " %");
// Gráficas
scf(1); clf(1);
plot(t, CA, 'b-', t, CB, 'm-', t, CC, 'g-');
xgrid; xlabel("Tiempo (h)"); ylabel("Concentración (kmol/m3)");
legend("CA", "CB", "CC");


scf(2); clf(2);
plot(t, T, 'r-', t, 330*ones(1, length(t)), 'k--');
xgrid; xlabel("Tiempo (h)"); ylabel("Temperatura (K)");
legend("T Reactor", "T = 330 K", 4);


scf(3); clf(3);
plot(Tj_rango, XA_finales, 'LineWidth', 2);
plot(Tj_optima, max_XA, 'ro', 'MarkerSize', 8);
xgrid; xlabel("Tj (K)"); ylabel("XA en t = 120h");
legend("XA vs Tj", "Tj Óptima", 4);