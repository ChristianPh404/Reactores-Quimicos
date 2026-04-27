clear; clc; 
// RDMP-Unknown-2.sce

function dxdt = f(t, x)
    CA = x(1);
    CB = x(2);
    T = x(3);
    TJ = x(4);
    
    // Cinética
    // En el enunciado original pone:
    // Factor preexponencial = 2.5E4 s^{-1}
    // k = 1E8*exp(-8000/T) s^{-1}
    // Entendemos que la fórmula explícita prevalece y ya incorpora preexponencial o Ea/R.
    k = A_pre * exp(-E_R / T);
    r = k * CA;
    
    // Balances de materia
    dCAdt = -r;
    dCBdt = r;
    
    // Calor transferido a la camisa
    Q = U * Area * (T - TJ);
    
    // Generación de calor por reacción termodinámica
    q_gen = r * V * (-H);
    
    // Balance de energía del reactor (J/s = W)
    dTdt = (q_gen - Q) / (V * rhoCp);
    
    // Balance de energía de la camisa (J/s = W)
    dTJdt = (Fc * rhoCp_c * (Tce - TJ) + Q) / (Vc * rhoCp_c);
    
    // Vector de derivadas
    dxdt(1) = dCAdt;
    dxdt(2) = dCBdt;
    dxdt(3) = dTdt;
    dxdt(4) = dTJdt;
endfunction

// Parámetros proporcionados (S.I.)
A_pre = 1E8; // s-1
E_R = 8000;  // K

H = -2E5; // J/mol

V = 10; // m^3
Vc = 1; // m^3
Area = 10; // m^2
U = 400; // J/(m^2 s K)
Fc = 0.01; // m^3/s
Tce = 280; // K

// Propiedades termodinámicas asumidas para "fase acuosa"
// Densidad 1000 kg/m3 y Cp 4184 J/(kg K)
rhoCp = 1000 * 4184; // J/(m^3 K)
rhoCp_c = 1000 * 4184; // J/(m^3 K)

// Condiciones iniciales
CA0 = 1000; // mol/m^3
CB0 = 0;
T0 = 280; // K
TJ0 = 280; // K
x0 = [CA0; CB0; T0; TJ0];

tfin = 10000; // s
dt = 10;
t = 0:dt:tfin;

// Integración ODE
disp("Calculando la dinámica del reactor...");
x = ode(x0, 0, t, f);

CA = x(1,:);
CB = x(2,:);
T = x(3,:);
TJ = x(4,:);

// Análisis final
disp("=== Resultados ===");
disp("Conversión final de A: " + string((CA0 - CA($))/CA0 * 100) + " %");
disp("Temperatura máxima del reactor: " + string(max(T)) + " K");
disp("Temperatura máxima de la camisa: " + string(max(TJ)) + " K");

// Gráficas
scf(1); clf(1);
plot(t, CA, 'b-', t, CB, 'm-', 'LineWidth', 2);
xgrid; xlabel("Tiempo (s)"); ylabel("Concentración (mol/m^3)");
legend("CA", "CB");


scf(2); clf(2);
plot(t, T, 'r-', t, TJ, 'b--', 'LineWidth', 2);
xgrid; xlabel("Tiempo (s)"); ylabel("Temperatura (K)");
legend("T Reactor", "T Camisa", 4);