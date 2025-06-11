clear; clc;
// 2024-06-10-A.sce

// RDMP-MULT
// 1) 2 A => B
// 2)   A => C   B> C
// No adiabático

// SISTEMA DE ECUACIONES ALGEBRAICAS
function dxdt = f(t,x)
    // Variables 
    CA = x(1);
    CB = x(2);
    CC = x(3);
    T  = x(4);

    // Selección de condiciones de operación
    if t < tfin/3 then
        TJ = TJ1;
        m = m1;
    else
        TJ = TJ2;
        m = m2;
    end

    // Ecuaciones de Arrhenius
    k1 = exp(-2500/(T+8)); // L/(h·mol)
    k2 = exp(-6000/(T+17)); // h-1

    // Velocidades de reacción
    r1 = k1*CA;
    r2 = k2*CB;

    // Balance de materia para A
    dCAdt = - r1;

    // Balance de materia para B
    dCBdt = r1 - r2;

    // Balance de materia para C
    dCCdt = r2;

    // Calor transferido del reactor a la camisa
    Q = m*CP*(TJ-T);

    // Balance de energía
    dTdt =  - (H1*r1 + H2*r2)/(RHO*CP) - Q/(V*RHO*CP);

    // Derivadas
    dxdt = [dCAdt; dCBdt; dCCdt; dTdt];
endfunction

// CONSTANTES
H1  = -75; // kj/mol
H2  = -50; // Kj/mol
CP  =  4.180; // Kj/(kg*K)
RHO =  1000; // g/L
V   =  0.1 // m^3
UA  =  200; // Kj/(h*K) 

// CONDICIONES INICIALES
CAini = 1000; // mol/m^3
CBini = 0; // mol/m^3
CCini = 0; // mol/m^3   
Tini  = 320; // K
xini  = [CAini;CBini;CCini;Tini];

// TIEMPO
tfin = 12; dt = 0.01; t = 0:dt:tfin; // h
TJ1 = 300; // K
TJ2 = 280; // K
m1 = 0.01*RHO //m^3/h * RHO // kg/h
m2 = 0.05*RHO // m^3/h * RHO // kg/h


// RESOLVER
x = ode(xini,0,t,f);
CA = x(1,:); CAfin = CA($);
CB = x(2,:); CBfin = CB($);
CC = x(3,:); CCfin = CC($);
T  = x(4,:); Tfin = T($);
disp("CA final: " + string(CAfin) + " mol/m^3");
disp("CB final: " + string(CBfin) + " mol/m^3");
disp("CC final: " + string(CCfin) + " mol/m^3");
disp("T final: " + string(Tfin) + " K");
// Graficar

scf(1); clf(1);
plot(t,CA,'r',t,CB,'g',t,CC,'b');
xlabel('Tiempo (h)');
ylabel('Concentración (mol/m^3)');
legend('CA','CB','CC');

indexCB = findmax(CB);
CBmax = CB(indexCB);
disp("CB max: " + string(CBmax) + " mol/m^3");
indexT = findT(T>320,1);
T_superior_a_320 = T-(indexT*dt);
disp("Tiempo en que T supera 320 K: " + string(T_superior_a_320) + " h");

//!Apartado C aunque esta mal ya que no estoy calculando bien la camisa
// seria de la siguiente manera

// indextj = find(T > 310, 1);
// Ttj = t*(T-(TJ(indextj)*dt));
//  disp("Tiempo en que T supera 310 K: " + string(Ttj) + " h");

//indexTJmax = findmax(TJ);
//disp("Tiempo en que TJ es máximo: " + string(t(indexTJmax)) + " h");
// trep = indextj:dt:T($);
// Tiempo= t(indextj:dt:tfin);

// Graficar temperatura
//scf(2); clf(2);
//plot(t,T,'k',Tiempo,trep,'r');
//xlabel('Tiempo (h)');
//ylabel('Temperatura (K)');
//legend('T,k','TJ> 310');
