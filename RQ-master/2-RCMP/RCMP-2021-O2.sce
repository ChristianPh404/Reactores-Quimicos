clear; clc; 
// 2021-04-23-A.sce

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt=f(t, x)
    // Variables
    CA = x(1)
    T = x(2)
    TJ = x(3)
    // Ecuación de Arrhenius
    k = k0*exp(-E/(R*T))
    // Velocidad de reacción
    r = k*CA
    // Calor transferido del reactor a la camisa
    Q = UA*(T-TJ)
    // Balance de materia para A
    // d(V*CA)dt = F*CA0 - F*CA - r*V
    dCAdt = F*(CA0-CA)/V - r
    // Balance de energía en el reactor
    // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - F*RHO*CP*T -H*r*V - Q
    dTdt = F*(T0-T)/V - H*r/(RHO*CP) - Q/(V*RHO*CP)
    // Balance de energía en la camisa
    // d(VJ*RHOJ*CPJ*TJ)dt = FJ*RHOJ*CPJ*TJ0 - FJ*RHOJ*CPJ*TJ + Q
    dTJdt = FJ*(TJ0-TJ)/VJ + Q/(VJ*RHOJ*CPJ)
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dTdt
    dxdt(3) = dTJdt
endfunction

// CONSTANTES
k0 = 1.5E12; // min-1
E = 8.5E4; // J/mol
R = 8.314 // J/(mol*K)
H = -3.1E5; // J/mol
V = 1000; // L
F = 1200; // L/min
CA0 = 2; // mol/L
T0 = 290; // K
VJ = 100; // L
UA = 2.5E6; // J/(min*K)
FJ = 600; // L/min
RHOJ = 1000; // g/L
CPJ = 4.18; // J/(g*K)
TJ0 = 280; //// K
RHO = 800; // g/L
CP = 3.5; // J/(g*K)

// CONDICIONES INICIALES
CAini = 0; // mol/L
Tini = 350; // K
TJini = TJ0; // K
xini = [CAini; Tini; TJini];

// TIEMPO
tfin = 10; dt = 0.01; t = 0:dt:tfin; // min

// RESOLVER
x = ode(xini,0,t,f);
CA = x(1,:); CAee = CA($)
T = x(2,:); Tee = T($)
TJ = x(3,:); TJee = TJ($)

// GRÁFICAS
scf(1); clf(1);
plot(t,CA,'b-');
xgrid; xtitle('A.sce','t','CA');

scf(2); clf(2);
plot(t,T,'r-',t,TJ,'g-');
xgrid; xtitle('A.sce','t','T(rojo), TJ(verde)');

// Resultados del estado estacionario
disp("=== Estado Estacionario (Tini = 350 K) ===");
disp("CAee = " + string(CAee) + " mol/L");
disp("Tee = " + string(Tee) + " K");
disp("TJee = " + string(TJee) + " K");

DT = T-TJ;

// (b)
DTobj = 5; // K
indexDTobj = find(DT<DTobj,1);
tDTobj = t(indexDTobj);
TDTobj = T(indexDTobj);
TJDTobj = TJ(indexDTobj);
CAindex = CA(indexDTobj);

disp('La diferencia de temperatura entre el reactor y la camisa es inferior a '+string(DTobj)+' K a partir del instante '+string(tDTobj)+' min')


// (c)
disp("=== Apartado C ===");
Tini_c = 420; // K
xini_c = [CAini; Tini_c; TJini];
x_c = ode(xini_c, 0, t, f);

T_c = x_c(2,:);
TJ_c = x_c(3,:);
DT_c = T_c - TJ_c;

[DTmin, indexDTmin] = min(DT_c);
tDTmin = t(indexDTmin);

disp("Para Tini = 420 K, la diferencia de temperatura mínima es " + string(DTmin) + " K en t = " + string(tDTmin) + " min");

scf(3); clf(3);
plot(t, DT, 'g-'); // Tini = 350 K
plot(t, DT_c, 'b-'); // Tini = 420 K
plot(tDTobj, DTobj, 'mo'); // Punto apartado b
plot(tDTmin, DTmin, 'ro'); // Punto apartado c
xgrid; xtitle('Evolución de DT (T - TJ)', 't (min)', 'DT (K)');
legend("Tini = 350 K", "Tini = 420 K", "Objetivo (b)", "Mínimo (c)", 4);
//estudiando la influencia de la temperatura inicial
disp("=== Influencia de Tini en la concentración y conversión final ===");
Tini_inf = 200:5:400;
CAee_inf = zeros(1, length(Tini_inf));
XAee_inf = zeros(1, length(Tini_inf));

for i = 1:length(Tini_inf)
    xini_inf = [CAini; Tini_inf(i); TJini];
    x_inf = ode(xini_inf, 0, t, f);
    
    CAee_val = x_inf(1, $); // Solo nos interesa CA
    CAee_inf(i) = CAee_val;
    XAee_inf(i) = (CA0 - CAee_val) / CA0; // Conversión de A
end

scf(6); clf(6);
plot(Tini_inf, CAee_inf, 'b-');
xgrid; xtitle('Influencia de Tini en la Concentración Final', 'Tini (K)', 'CAee (mol/L)');

scf(7); clf(7);
plot(Tini_inf, XAee_inf, 'r-');
xgrid; xtitle('Influencia de Tini en la Conversión', 'Tini (K)', 'XA');