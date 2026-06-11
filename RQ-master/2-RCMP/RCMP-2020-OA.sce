clc;
clear all;
// A => B
// B => C
// B => D
// SISTEMA DE ECUACIONES ALGEBRAICAS
function dxdt=f(x)
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
  // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - F*RHO*CP*T - H1*r1*V - H2*r2*V - H3*r3*V - Q
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
H2 = -500 ; // J/mol  //! Revisar unidades
R = 1.9872; // cal/(mol*K)
k03 = 1.09E3; // 1/min
H3 = -2000 // cal/mol
E3 = 4371 // cal/mol

// SOLUCIÓN SUPUESTA
CAeeguess = 1; // mol/L
CBeeguess = 0; // mol/L
CCeeguess = 0; // mol/L
CDeeguess = 0; // mol/L
Teeguess = 300; // K
xeeguess = [CAeeguess;CBeeguess;CCeeguess;CDeeguess;Teeguess];

// RESOLVER
[xee,fxee,info] = fsolve(xeeguess,f)
CAee = xee(1)
CBee = xee(2)
CCee = xee(3)
CDee = xee(4)
Tee = xee(5)

// OPTIMIZAR LA TEMPERATURA

// GRÁFICAS
scf(1);
plot(CAee,Tee,'x');

// LINEALIZACIÓN
J = numderivative(f,xee); // Jacobiano

// ESTABILIDAD DE UN SISTEMA LINEAL DE ECUACIONES DIFERENCIALES
// Valores propios
lambda = spec(J)

// Critero de estabilidad
Estable = and(real(lambda) < 0)
