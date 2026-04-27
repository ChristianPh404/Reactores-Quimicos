clear; clc;
// RDMP-O-2022.sce
// A <=> B 
// Reacción reversible exotérmica en reactor discontinuo mezcla perfecta

function dxdt = f(t,x)
    // Variables de estado
    CA = x(1);
    CB = x(2);
    T = x(3);
    
    // Ecuaciones de la cinética (Arrhenius y Van't Hoff)
    k = A0 * exp(-Ea/(R*T));
    Keq = Keq0 * exp(-H/(R*T));
    // Ecuación de velocidad reacción elemental
    r = k * (CA - CB/Keq);
    
    // Calor extraído (flujo de energía térmica a la camisa)
    Q = UA * (T - Tj);
    
    // Balances
    dCAdt = -r;
    dCBdt = r;
    // d(V*rho*Cp*T)/dt = -H*r*V - Q
    dTdt = (-H * r) / (rho * Cp) - Q / (V * rho * Cp);
    
    // Vector de derivadas
    dxdt(1) = dCAdt;
    dxdt(2) = dCBdt;
    dxdt(3) = dTdt;
endfunction

// Constantes dadas por el problema
A0 = 9.1E14; // min-1
Keq0 = 4.21E-5; 
Ea = 1.12E5; // J/mol
H = -1.46E5; // J/mol
R = 8.314; // J/(mol K)

V = 2000; // L
Tj = 280; // K
UA = 1500; // J/(min K)
rho = 0.9; // kg/L
Cp = 4500; // J/(kg K)

// Condiciones iniciales
CAini = 1; // mol/L
CBini = 0; // mol/L
Tini = 360; // K
xini = [CAini; CBini; Tini];

// Configuración del tiempo (lo ponemos alto para asegurar que llegue a la conversión)
tfin = 5000; 
dt = 0.1;
t = 0:dt:tfin;

// Resolución numérica de la ODE
x = ode(xini, 0, t, f);
CA = x(1,:);
CB = x(2,:);
T = x(3,:);
XA = 1 - CA/CAini;

// RESULTADOS
// === Apartado A ===
XAobj = 0.95;
index_95 = find(XA >= XAobj, 1);

if isempty(index_95) then
    disp("No se alcanza el 95% de conversión en " + string(tfin) + " min");
    t_95 = t($);
    index_95 = length(t);
else
    t_95 = t(index_95);
    T_95 = T(index_95);
    disp("=== Apartado A ===");
    disp("El 95% de conversión se alcanza en el tiempo t = " + string(t_95) + " min");
    disp("Temperatura en este punto = " + string(T_95) + " K");
end

// Acotamos los datos hasta el 95% para unicamente mostrar hasta el 95%
t_graf = t(1:index_95);
CA_graf = CA(1:index_95);
CB_graf = CB(1:index_95);
T_graf = T(1:index_95);
XA_graf = XA(1:index_95);

// === Apartado B ===
// Punto de inflexión en la evolución temporal de la temperatura (donde la derivada cambia)
dTdt = diff(T_graf)/dt;
[max_dTdt, index_inflexion] = max(dTdt);
t_inflexion = t_graf(index_inflexion);
T_inflexion = T_graf(index_inflexion);

disp("=== Apartado B ===");
disp("El punto de inflexión en temperatura ocurre en t = " + string(t_inflexion) + " min");
disp("Temperatura en el punto de inflexión = " + string(T_inflexion) + " K");
//Bonus donde CA = CB
index_CA_CB = find(CA < CB, 1);
if isempty(index_CA_CB) then
    disp("No se alcanza la igualdad de concentraciones en " + string(tfin) + " min");
    t_CA_CB = t($);
    index_CA_CB = length(t);
else
    t_CA_CB = t(index_CA_CB);
    T_CA_CB = T(index_CA_CB);
    disp("=== Apartado Bonus ===");
    disp("La igualdad de concentraciones se alcanza en el tiempo t = " + string(t_CA_CB) + " min" + " y la concentracion es = " + string(CA(index_CA_CB)));
    disp("Temperatura en este punto = " + string(T_CA_CB) + " K");
end
// Gráficas
scf(1); clf(1);
plot(t_graf, CA_graf, 'b-', t_graf, CB_graf, 'm-', t_CA_CB, CA(index_CA_CB), 'ko');
xgrid; xlabel('Tiempo (min)'); ylabel('Concentración (mol/L)');
legend('CA', 'CB','CA=CB',-1);

scf(2); clf(2);
plot(t_graf, T_graf, 'r-', t_inflexion, T_inflexion, 'ko');
xgrid; xlabel('Tiempo (min)'); ylabel('Temperatura (K)');
legend('T', 'Punto Inflexión',2);

scf(3); clf(3);
plot(t_graf, XA_graf, 'g-', t_graf($), XA_graf($), 'ro');
xgrid; xlabel('Tiempo (min)'); ylabel('Conversión XA');
legend('XA', '95%',2);