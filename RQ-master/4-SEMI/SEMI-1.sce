clear; clc; 
// SEMI-1.sce 
// 1) A + B => P*
// 2) A + B => Q
// Isotermo

// (a)

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = f(t,x)
    // Variables diferenciales
    V  = x(1)
    NA = x(2)
    NB = x(3)
    NP = x(4)
    NQ = x(5)
    // Concentraciones
    CA = NA/V
    CB = NB/V
    CP = NP/V
    CQ = NQ/V
    // Velocidades de reacción    
    r1 = k1*CA*CB
    r2 = k2*CA*CB^2
    // Caudal de alimentación
    if t < tB then 
        F = FB;  // Semicontinuo
        // ejemplo caudal de forma lineal F = 2*VB /tfin -2 *VB/tfin^2*t
    else 
        F = 0;   // Discontinuo
    end
    // Balance de materia global
    // d(V*RHO)dt = F*RHO
    dVdt = F
    // Balance de materia para A    
    dNAdt = (-r1-r2)*V
    // Balance de materia para B    
    dNBdt= F*CB0 + (-r1-r2)*V
    // Balance de materia para P    
    dNPdt = r1*V
    // Balance de materia para Q    
    dNQdt = r2*V
    // Derivadas
    dxdt(1) = dVdt
    dxdt(2) = dNAdt
    dxdt(3) = dNBdt
    dxdt(4) = dNPdt
    dxdt(5) = dNQdt
endfunction

// CONSTANTES
k1 = 1; // L/(mol*min)
k2 = 2; // L2/(mol2*min)
VB = 1; // L
CB0 = 1; // mol/L

// CONDICIONES INICIALES
Vini = 0.5; // L
CAini = 2; // mol/L
NAini = Vini*CAini; NBini = 0; NPini = 0; NQini = 0; // mol
xini = [Vini;NAini; NBini; NPini; NQini];

// TIEMPO
tfin = 50; dt = 0.01; t = 0:dt:tfin; // min

// OPTIMIZAR tB para MAXIMIZAR NPfin
tBinterval = 0.5:0.5:50; // min

for i = 1:length(tBinterval)
    tB = tBinterval(i);
    FB = VB/tB; // L/min
    // RESOLVER
    x = ode(xini,0,t,f);
    NPfin(i) = x(4,$);
end

[NPfinmax,indexNPfinmax] = max(NPfin)
tBopt = tBinterval(indexNPfinmax)


// GRÁFICAS
scf(1); clf(1);
plot(tBinterval,NPfin,'ro');
plot(tBopt,NPfinmax,'x');
xgrid; xlabel('tB'); ylabel('NPfin');

// (b)

tB = tBopt;
FB = VB/tB;
disp("tiempo optimo tB: " + string(tBopt) + " min");
disp("Caudal optimo FB: " + string(FB) + " L/min");
// RESOLVER
x = ode(xini,0,t,f);
V =  x(1,:); Vfin = V($);
NA = x(2,:); NAfin = NA($);
NB = x(3,:); NBfin = NB($);
NP = x(4,:); NPfin = NP($);
NQ = x(5,:); NQfin = NQ($);
disp("Vfinal: " + string(Vfin) + " L");
disp("NAfinal: " + string(NAfin) + " mol" + " (CAfinal: " + string(NAfin/Vfin) + " mol/L)");
disp("NBfinal: " + string(NBfin) + " mol" + " (CBfinal: " + string(NBfin/Vfin) + " mol/L)");
disp("NPfinal: " + string(NPfin) + " mol" + " (CPfinal: " + string(NPfin/Vfin) + " mol/L)");
disp("NQfinal: " + string(NQfin) + " mol" + " (CQfinal: " + string(NQfin/Vfin) + " mol/L)");
// GRÁFICAS
scf(2); clf(2);
plot(t,V);
xgrid; xlabel('t'); legend('V',-2,%f);

scf(3); clf(3);
plot(t,NA,t,NB,t,NP,t,NQ);
xgrid; xlabel('t'); legend('NA','NB','NP','NQ',-2,%f);

Sg = x(4,$) / (x(5,$)+x(4,$));              // NPfin / NQfin
XA = (NAini - x(2,$)) / NAini;    // Conversión de A
s = x(4,$) / x(5,$);
disp("Selectividad global final: " + string(Sg));
printf("Conversión final de A: %.2f%%\n", XA);
XA_t = (NAini - NA) / NAini;
scf(4); clf(4);
plot(t, XA_t);
xgrid; xlabel('t'); ylabel('XA');
scf(6); clf(6);
plot(t, s);
xgrid; xlabel('t'); ylabel('selectividad P/Q');
R = XA*Sg;
printf("El rendimiento es: %.2f%%\n", R);

//Representacion de las concentraciones
CA = NA./V;
CB = NB./V;
CP = NP./V;
CQ = NQ./V; 
//TODO retocar para que  sea normalizado 
tol = 1e-5;
dif = abs(CA - CB);
indexCA = find(dif <= tol, 1); 

//* Grafica 

scf(5); clf(5);
plot(t,CA,'k-.',t,CB,t,CP,t,CQ,t(indexCA),CA(indexCA),'ro');
xgrid; xlabel('t'); legend('CA','CB','CP','CQ','CA=CB',-2,%f);
title ('concentraciones frente al tiempo');

//if length(indexCA) == 0
//    printf("No se iguala la concentracion de A y B.\n");
//else
//    printf("A y B se igualan en t = %.2f min\n", t(indexCA));
//end
