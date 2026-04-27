clear; clc;
// RDMP-O-2018.sce
// A => B
// Reacción química exotérmica en fase líquida

function dxdt = f(t,x)
    // Variables diferenciales
    CA = x(1);
    T = x(2);
    TJ = x(3);
    // Ecuación de Arrhenius
    k = A0 * exp(-Ea/(R*T));
    // Velocidad de reacción
    r = k * CA;
    // Calor transferido del reactor a la camisa
    Q = U * A_area * (T - TJ);
    // Balance de materia para A
    dCAdt = -r;
    // Balance de energía del reactor
    // d(V*rho*Cp*T)/dt = -H*r*V - Q
    dTdt = (-H * r) / (rho * Cp) - Q / (V * rho * Cp);
    // Balance de energía de la camisa
    // d(Vj*rho*Cp*TJ)/dt = Fj*rho*Cp*(Tjin - TJ) + Q
    dTJdt = (Fj / Vj) * (Tjin - TJ) + Q / (Vj * rho * Cp);
    // Derivadas
    dxdt(1) = dCAdt;
    dxdt(2) = dTdt;
    dxdt(3) = dTJdt;
endfunction

// CONSTANTES CINETICAS Y TERMODINAMICAS
A0 = 2.5E4; // s-1 (Factor preexponencial)
Ea = 42000; // J/mol (Energía de activación)
H = -5E5; // J/mol (Entalpía de reacción)
R = 8.314; // J/(mol*K) (Constante de los gases)

// CONSTANTES DEL REACTOR Y CAMISA
V = 5; // m3 (Volumen del reactor)
Vj = 1; // m3 (Volumen de la camisa)
A_area = 10; // m2 (Área de intercambio de calor)
Fj = 0.01; // m3/s (Caudal de agua de refrigeración)
Tjin = 280; // K (Temperatura de entrada del agua)
U = 500; // J/(m2*s*K) (Coeficiente de transmisión de calor)

// PROPIEDADES FISICAS DEL AGUA
rho = 1000; // kg/m3 (Densidad)
Cp = 4180; // J/(kg*K) (Capacidad calorífica)

// CONDICIONES INICIALES
CAini = 500; // mol/m3
Tini = 280; // K
TJini = 280; // K
xini = [CAini; Tini; TJini];

// TIEMPO
tfin = 1200; // s
dt = 1; 
t = 0:dt:tfin; 

// RESOLVER
x = ode(xini, 0, t, f);
CA = x(1,:); CAfin = CA($);
T = x(2,:); Tfin = T($);
TJ = x(3,:); TJfin = TJ($);

// A - Calcular la dinámica del proceso
disp("=== Apartado A: Dinámica del proceso ===");
disp("Concentración final de A = " + string(CAfin) + " mol/m3");
disp("Temperatura final del reactor = " + string(Tfin) + " K");
disp("Temperatura final de la camisa = " + string(TJfin) + " K");

// B - Localizar el 50% de conversión
disp("=== Apartado B: Localizar el 50% de conversión ===");
XA = 1 - CA/CAini;
XAobj = 0.5;
indexXA = find(XA > XAobj, 1);
tXA = t(indexXA);
disp("El 50% de conversión se alcanza en el tiempo t = " + string(tXA) + " s");

// C - Determinar durante cuánto tiempo la temperatura está comprendida entre 300 y 310 K
disp("=== Apartado C: Temperatura comprendida entre 300 y 310 K ===");
indexT = find(T > 300 & T < 310);
tiempoT = length(indexT) * dt;
disp("Tiempo con la temperatura del reactor entre 300 y 310 K = " + string(tiempoT) + " s");

// GRÁFICAS
scf(1); clf(1);
plot(t, CA, 'b-');
xgrid; xlabel('Tiempo (s)'); ylabel('Concentración de A (mol/m3)');
legend('CA');

scf(2); clf(2);
plot(t, T, 'r-', t, TJ, 'm-',t(indexT), T(indexT), 'go');
xgrid; xlabel('Tiempo (s)'); ylabel('Temperatura (K)');
legend('T Reactor', 'T Camisa', 2);

scf(3); clf(3);
plot(t, XA, 'g-', tXA, XAobj, 'ro');
xgrid; xlabel('Tiempo (s)'); ylabel('Conversión XA');
legend('XA', '50%', 4);