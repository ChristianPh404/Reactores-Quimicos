clear; clc; 
// RDMP-4.sce
// 2 A <=> B
// Isotermo

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = f(t,x)
    // Variables diferenciales
    CA = x(1)
    CB = x(2)
    // Velocidad de reacción
    // r = rd - ri = kd*CA^2 - ki*CB = kd*CA^2 - kd*CB/Keq
    r = kd*(CA^2 - CB/Keq)
    // Balance de materia para A
    // d(V*CA)dt = -r*V
    dCAdt = -r 
    // Balance de materia para B
    // d(V*CB)dt = r/2*V
    dCBdt = r/2
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dCBdt
endfunction

// CONSTANTES
T = 420; // K
kd0 = 1.4E12; // L/(mol*h)
E = 105000; // J/mol
Keq0 = 6.9E8; // L/mol
H = 63000; // J/mol
R = 8.314; // J/(mol*K)
kd = kd0*exp(-E/(R*T))   // Ecuación de Arrhenius
Keq = Keq0*exp(-H/(R*T)) // Ecuación de Van't Hoff

// CONDICIONES INICIALES
CAini = 5; CBini = 0; // mol/L
xini = [CAini; CBini];

// TIEMPO
tfin = 100; dt = 0.1; t = 0:dt:tfin; // h

// RESOLVER
x = ode(xini,0,t,f);
CA = x(1,:); CAfin = CA($);
CB = x(2,:); CBfin = CB($);
disp('Concentracion final de A: '+string(CAfin)+' mol/L');
disp('Concentracion final de B: '+string(CBfin)+' mol/L');

dCAdt = diff(CA)/dt; dCBdt = diff(CB)/dt; 
dCAdteq = 1E-5; dCBdteq = 1E-5; // mol/(L*h)
indexCAeq = find(abs(dCAdt)<dCAdteq,1); indexCBeq = find(abs(dCBdt)<dCBdteq,1);
tCAeq = t(indexCAeq), tCBeq = t(indexCBeq);
CAeq = CA(indexCAeq), CBeq = CB(indexCBeq);
disp('Tiempo de equilibrio para CA: '+string(tCAeq)+' h');
disp('Concentracion de equilibrio para CA: '+string(CAeq)+' mol/L');
disp('Tiempo de equilibrio para CB: '+string(tCBeq)+' h');
disp('Concentracion de equilibrio para CB: '+string(CBeq)+' mol/L');
teq = max(tCAeq,tCBeq);
disp('Tiempo de equilibrio: '+string(teq)+' h');

indexCACB = find(CA<CB,1);
tCACB = t(indexCACB);
CACB = CA(indexCACB);
disp('Tiempo en que CA = CB: '+string(tCACB)+' h');

// GRÁFICAS
scf(1); clf(1); 
plot(t,CA,t,CB);
plot(tCAeq,CAeq,'b.',tCBeq,CBeq,'g.');
plot(tCACB,CACB,'bo');
xgrid; xlabel('t'); legend('CA','CB',-2,%f);

scf(2); clf(2); 
plot(t(1:$-1),abs(dCAdt),t(1:$-1),abs(dCBdt));
plot(tCAeq,dCAdteq,'b.',tCBeq,dCBdteq,'g.'); 
xgrid; xlabel('t'); legend('|dCAdt|','|dCBdt|',-2,%f);
a2 = gca; a2.log_flags = "nl";
//! gca da el control de los ejes del grafico
//! log_flags = "nl" es para que el eje y sea logaritmico
