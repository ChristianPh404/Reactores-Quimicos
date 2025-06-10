clear; clc; 
// 2024-07-04-A.sce

// RCMP
// 1) A => B
// 2) B => C

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = f(t,x)
    // Variables
    CA = x(1)
    CB = x(2)
    CC = x(3)
    T  = x(4)
    // Ecuación de Arrhenius    
    k1 = exp(6-2500/T)  // h-1
    k2 = exp(15-6000/T) // h-1
    // Velocidad de reacción
    r1 = k1*CA
    r2 = k2*CB
    // Balance de materia para A
    // d(V*CA)dt = F*CA0 - F*CA - r1*V
    dCAdt = F*(CA0-CA)/V - r1
    // Balance de materia para B
    // d(V*CB)dt = F*CB0 - F*CB + r1*V - r2*V
    dCBdt = F*(CB0-CB)/V + r1 - r2
    // Balance de materia para C
    // d(V*CC)dt = F*CC0 - F*CC + r2*V
    dCCdt = F*(CC0-CC)/V + r2
    // Temperatura de la camisa
    if t < tee then TJ = 300;
    else TJ = 300 + 10*sin(2*%pi*t/10);
    end
    // Calor transferido del reactor a la camisa
    Q = UA*(T-TJ)
    // Balance de energía
    // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - F*RHO*CP*T - H1*r1*V - H2*r2*V - Q
    dTdt  = F*(T0-T)/V - (H1*r1+H2*r2)/(RHO*CP) - Q/(V*RHO*CP)
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dCBdt
    dxdt(3) = dCCdt
    dxdt(4) = dTdt
endfunction

// CONSTANTES
H1 = -35;  // kcal/mol
H2 = -18;  // kcal/mol
RHO = 1; // kg/L
CP = 1; // kcal/(kg*K)
V = 147; //L  (Tantear para CBee = 0.5)
F = 15; // L/h
T0 = 300; // K
CA0 = 1; CB0 = 0; CC0 = 0; // mol/L
UA = 6; // kcal/(h*K)

// CONDICIONES INICIALES
CAini = 1; CBini = 0; CCini = 0; // mol/L 
Tini = 300; // K
xini = [CAini; CBini; CCini; Tini];

// TIEMPO
tee = 200; tfin = 400; dt = 0.1; t = 0:dt:tfin; // h

// RESOLVER
x = ode(xini,0,t,f);
CA = x(1,:); CAee = CA(t==tee)
CB = x(2,:); CBee = CB(t==tee)
CC = x(3,:); CCee = CC(t==tee)
T  = x(4,:); Tee  = T(t==tee)

// GRÁFICAS
scf(1); clf(1);
plot(t,CA,t,CB,t,CC);
xgrid; xlabel('t'); legend('CA','CB','CC',-2,%f);

scf(2); clf(2);
plot(t,T,);
xgrid; xlabel('t'); legend('T',-2,%f);

// Efecto de la perturbación
tp = tee:dt:tfin;
TJp = 300 + 10*sin(2*%pi*tp/10);

scf(3); clf(3);
subplot(511); plot(tp,TJp,'g'); xgrid; ylabel('TJ');
subplot(512); plot(tp,CA(t>=tee)); xgrid; ylabel('CA');
subplot(513); plot(tp,CB(t>=tee)); xgrid; ylabel('CB');
subplot(514); plot(tp,CC(t>=tee)); xgrid; ylabel('CC');
subplot(515); plot(tp, T(t>=tee)); xgrid; xlabel('t'); ylabel('T');
