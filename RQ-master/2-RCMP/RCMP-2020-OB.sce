clc;
clear all;
// A => B
// B => C
// B => D
// SISTEMA DE ECUACIONES ALGEBRAICAS
function dxdt=f(t, x)
  // Variables
  CA = x(1)
  CB = x(2)
  CC = x(3)
  CD = x(4)
  T  = x(5)
  // Ecuaciones de Arrhenius
  k1 = k01*exp(-E1/(R*T))
  k2 = k02*exp(-E2/(R*T))
  k3 = k03*exp(-E3/(R*T))
  // Velocidades de reacción
  r1 = k1*CA
  r2 = k2*CB
  r3 = k3*CB
  // Balance de materia para A
  // d(V*CA)dt = F*CA0 - F*CA - r1*V
  dCAdt = F*(CA0-CA)/V - r1
  // Balance de materia para B
  // d(V*CB)dt = F*CB0 - F*CB - r2*V - r3*V
  dCBdt = F*(CB0-CB)/V + r1- r2 - r3
  // Balance de materia para C
  dCCdt = F*(CC0-CC)/V + r2
  // Balance de materia para D
  dCDdt = F*(CD0-CD)/V + r3
  // Calor transferido del reactor a la camisa
  Q = U*(T-TJ)
  // Balance de energía
  // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - F*RHO*CP*T - H1*r1*V - H2*r2*V - H3*r3*VQ
  dTdt = F*(T0-T)/V - (H1*r1 + H2*r2 + H3*r3)/(RHO*CP) - Q/(V*RHO*CP)
  // Derivadas
  dxdt(1) = dCAdt
  dxdt(2) = dCBdt
  dxdt(3) = dCCdt
  dxdt(4) = dCDdt
  dxdt(5) = dTdt
endfunction

// CONSTANTES
V = 150; // L
F = 5; // L/min
CA0 = 1; // mol/L
CB0 = 0; // mol/L
CC0 = 0; // mol/L
CD0 = 0; // mol/L
T0 = 300; // K
TJ = 350; // K
U = 700; // cal/(min*K)
CP = 1; // cal/(g*K)
RHO = 1000; // g/L
k01 = 4.03E2; // 1/min
E1 = 5166; // cal/mol
H1 = -1000; // cal/mol
k02 = 2.41E7; // 1/min
E2 = 12319; // cal/mol
H2 = -500; // J/mol
R = 1.9872; // cal/(mol*K)
k03 = 1.09E3; // 1/min
H3 = -2000 // cal/mol
E3 = 4371 // cal/mol

// CONDICIONES INICIALES
CAini = 1; // mol/L
CBini = 0; // mol/L
CCini = 0; // mol/L
CDini = 0; // mol/L
Tini = 300; // K
xini = [CAini;CBini;CCini;CDini;Tini];

//TIEMPO
tfin=1000; dt=10; t=0:dt:tfin;

// RESOLVER
x = ode(xini,0,t,f);
xfin = x(:,$)
dxdT= f(tfin,xfin)
Estacionario = abs(dxdT ./ xfin) < 0.001

CA = x(1,:); CAee = CA($)
CB = x(2,:); CBee = CB($)
CC = x(3,:); CCee = CC($)
CD = x(4,:); CDee = CD($)
T = x(5,:); Tee = T($)

scf(1); clf(1);
plot(t,CA,t,CB,t,CC,t,CD);
xgrid; xtitle('t','CA(azul), CB(verde), CC(rojo), CD(cian)');

// DURANTE CUANTO TIEMPO PREDOMINA A
indexA= find(CA == max(CA,CB,CC,CD));
tD = dt*length(indexA)
plot(t(indexA),CD(indexA),'go');

//EL MOMENTO DE MÁX CONCENTRACIÓN
