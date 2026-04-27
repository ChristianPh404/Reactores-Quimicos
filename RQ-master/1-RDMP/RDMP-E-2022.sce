clear; clc; 
// RDMP-E-2022.sce
// Reacciones en paralelo 
// 1) 2A -> B
// 2) A -> C 

Qc = 1
TJ0 = 300
error_tolerado = 0.01;
max_iter = 1000; // Para evitar bucles infinitos

for i = 1:max_iter
function dxdt = f(t,x)
    CA = x(1);
    CB = x(2);
    CC = x(3);
    T = x(4);
    TJ = x(5);
    // Arrhenius
    R = 1.987e-3; // cal/(mmol K)
    k1 = A1 * exp(-E1/(R*T));
    k2 = A2 * exp(-E2/(R*T));
    
    r1 = k1 * CA^2;
    r2 = k2 * CA;
    
    // Balances de materia
    dCAdt = -2*r1 - r2;
    dCBdt = r1;
    dCCdt = r2;
    
    // Balance de energía del reactor
    Q = U * Vj * (T - TJ);
    dTdt = ( (-H1)*r1 + (-H2)*r2 - (Q/V) ) / rhoCp;
    
    // Balance de energía de la camisa
   
    dTJdt = (Q / (Vj * rhoCp)) + Qc * (TJ0 - TJ);
    
    dxdt(1) = dCAdt;
    dxdt(2) = dCBdt;
    dxdt(3) = dCCdt;
    dxdt(4) = dTdt;
    dxdt(5) = dTJdt;
endfunction

// Constantes
A1 = 6.1E17; // cm3/(h*mmol)
E1 = 28; // cal/mmol
H1 = -25; // cal/mmol

A2 = 5.7E14; // h-1
E2 = 21; // cal/mmol
H2 = -20; // cal/mmol

U = 10; // cal/(cm3*h*K)

V = 125; // cm3
Vj = 75; // cm3
rhoCp = 1; // cal/(cm3*K)
// Condiciones iniciales
CAini = 1; // mmol/cm3
CBini = 0; // mmol/cm3
CCini = 0; // mmol/cm3
Tini = 310; // K
TJini = 310; // K

xini = [CAini; CBini; CCini; Tini; TJini];
tfin = 2.5; // h
dt = 0.001;
t = 0:dt:tfin;
// Resolvemos la ODE completa para este QC
    sol = ode(xini, 0, t, f);
    TJ_sim = sol(5,:);
    T_pico = max(TJ_sim);
    err = T_pico - 312;
    Kp = 0.5
    // Lógica de ajuste que propusiste
    if abs(err) < error_tolerado then
        disp("Caudal (Qc) final: " + string(Qc) + ". Pico de TJ: " + string(T_pico) + " K" + " exito en la iteracion " + string(i));
        break; // Salimos del bucle si ya es 312
    else
        Qc = Qc +(err * Kp); // Si calienta mucho, quitamos caudal
    end
end
disp("El valor de Qc es: " + string(Qc) + " cm3/min");
// Integración
x = ode(xini, 0, t, f);
CA = x(1,:);
CB = x(2,:);
CC = x(3,:);
T = x(4,:);
TJ = x(5,:);

// Apartado A: Tiempo con CA < 0.1 mmol/cm3
index_CA = find(CA < 0.1, 1);
if isempty(index_CA) then
    disp("Apartado A: La concentración de A no baja de 0.1 mmol/cm3 en el tiempo simulado.");
else
    t_obj = t(index_CA);
    tCalc = tfin - t_obj;
    disp("Apartado A: La concentración de A fue inferior a 0.1 durante " + string(tCalc) + " horas.");
end

// Apartado B: Extremos de temperatura
// Punto de temperatura máxima  
[Tmax, index_Tmax] = max(T);
t_Tmax = t(index_Tmax);

[Tmin, index_Tmin] = min(T);
t_Tmin = t(index_Tmin);

disp("La temperatura máxima alcanzada es " + string(Tmax) + " K en t = " + string(t_Tmax) + " h.");
disp("La temperatura final del reactor es " + string(T($)) + " K." + "Temperatura final de la camisa es " + string(TJ($)) + " K.");

// Gráficas
scf(1); clf(1);
plot(t, CA, 'b-', t, CB, 'g-', t, CC, 'r-');
xgrid; xlabel('Tiempo (h)'); ylabel('Concentración (mmol/cm^3)');
legend('CA', 'CB', 'CC',-1);


scf(2); clf(2);
plot(t, T, 'r-', t, TJ, 'b--');
plot(t_Tmax, Tmax, 'ko');
xgrid; xlabel('Tiempo (h)'); ylabel('Temperatura (K)');
legend('T Reactor', 'T Camisa', 'T Max',3);