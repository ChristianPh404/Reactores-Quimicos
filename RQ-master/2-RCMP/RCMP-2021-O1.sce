clear; clc; 
// 2021-06-04-A.sce
// https://en.wikipedia.org/wiki/Brusselator

// A => X        r1 = k1*CA
// 2X + Y => 3X  r2 = k2*CX^2*CY
// B + X => Y + D  r3 = k3*CB*CX
// X => E        r4 = k4*CX

// CONSTANTES
k1 = 1; // min-1
k2 = 1; // L2/(mol2*min)
k3 = 1; // L/(mol*min)
k4 = 1; // min-1
CA = 1.5; // mol/L
CB = 3; // mol/L

// a) BALANCES GRÁFICOS

scf(1); clf(1);
xgrid; xtitle('A.sce','CX','CY');

CXmin = 0; dCX = 0.01; CXmax = 1; // mol/L
CYmin = 0; dCY = 0.01; CYmax = 1; // mol/L

CX = CXmin:dCX:CXmax;

// Balance de materia para X
// dCXdt = r1 + r2 - r3 - r4 = 
//       = k1*CA + k2*CX^2*CY - k3*CB*CX - k4*CX = 0
CY = (-k1*CA + k3*CB*CX + k4*CX)./(k2*CX^2);
plot(CX,CY,'r-');

// Balance de materia para Y
// dCYdt = - r2 + r3 =
//       = - k2*CX^2*CY + k3*CB*CX = 0
CY = k3*CB*CX./(k2*CX.^2);
plot(CX,CY,'r--');

// b) ESTADO ESTACIONARIO

// Sistema de ecuaciones algebraicas
function dxdt=f(x)
    // Variables diferenciales
    CX = x(1)
    CY = x(2)
    // Balance de materia para X
    dCXdt = k1*CA + k2*CX^2*CY - k3*CB*CX - k4*CX
    // Balance de materia para Y
    dCYdt = - k2*CX^2*CY + k3*CB*CX
    // Derivadas
    dxdt(1) = dCXdt
    dxdt(2) = dCYdt
endfunction

// Solución supuesta
CXeeguess = 2; CYeeguess = 2; // mol/L
xeeguess = [CXeeguess;CYeeguess];

// Resolver
[xee,fxee,info] = fsolve(xeeguess,f)
CXee = xee(1)
CYee = xee(2)
plot(CXee,CYee,'x');

J = numderivative(f,xee) // Jacobiano

lambda = spec(J) // Valores propios
Estable = and(real(lambda) < 0)

// c) CAMPO VECTORIAL

// Sistema de ecuaciones diferenciales
function dxdt=g(t, x)
    dxdt = f(x)
endfunction

fchamp(g,0,CXmin:25*dCX:CXmax,CYmin:25*dCY:CYmax);

// d) DINÁMICA

// Condiciones iniciales
CXini = 0; CYini = 0; // mol/L
xini = [CXini;CYini];

// Tiempo
tfin = 60; dt = 0.01; t = 0:dt:tfin; // min

// Resolver
x = ode(xini,0,t,g);
CX = x(1,:); CXee = CX($)
CY = x(2,:); CYee = CY($)
plot(CX,CY,'o-');

CXobj = 2; CYobj = 2; // mol/L
indexobj = find(CX>CXobj & CY>CYobj);
plot(CX(indexobj),CY(indexobj),'g.');
tobj = dt*length(indexobj)

a1 = gca;
a1.data_bounds = [CXmin,CYmin;CXmax,CYmax];
