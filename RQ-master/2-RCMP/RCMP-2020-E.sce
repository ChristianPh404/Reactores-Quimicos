clear; clc;
// A + B <=> C
// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt=f(t, x)
  // Variables
  CA = x(1)
  CB = x(2)
  CC = x(3)
  T  = x(4)
  // Ecuación de Arrhenius
  k = k0*exp(-E/(R*T))
  // Ecuación de Van't Hoff
  Keq = Keq0*exp(-H/(R*T))
  // Velocidad de reacción
  r = k*(CA*CB-CC/Keq)
  // Calor transferido del reactor a la camisa
  Q = UA*(T-TJ)
  // Balance de materia para A
  // d(V*CA)dt = F*CA0 - F*CA - r*V
  dCAdt = F*(CA0-CA)/V - r
  // Balance de materia para B
  // d(V*CB)dt = F*CB0 - F*CB - r*V
  dCBdt = F*(CB0-CB)/V - r
  // Balance de materia para C
  // d(V*CC)dt = F*CC0 - F*CC + r*V
  dCCdt = F*(CC0-CC)/V + r
  // Balance de energía
  // d(V*RHO*CP*T)dt = F*RHO*CP*T0 - F*RHO*CP*T -H*r*V - Q
  dTdt = F*(T0-T)/V - H*r/(RHO*CP) - Q/(V*RHO*CP)
  // Derivadas
  dxdt(1) = dCAdt
  dxdt(2) = dCBdt
  dxdt(3) = dCCdt
  dxdt(4) = dTdt
endfunction

// CONSTANTES
F = 10; // L/min
CA0 = 1; CB0 = 1.5; CC0 = 0; // mol/L
V = 1500; // L
k0 = 1.6E4; // L/(mol*min)
Keq0 = 0.11; // L/mol
E = 3.2E4; // J/mol
H = -1.2E4; // J/mol
TJ = 272; // K
T0 = 320; // K
UA = 1.5E4; // J/(min*K)
R = 8.314; // J/(mol*K)
CP = 3.8; // J/(g*K)
RHO = 1000; // g/L

// CONDICIONES INICIALES
CAini = 1; CBini = 1.5; CCini = 0; // mol/L
Tini = 320;
xini = [CAini; CBini; CCini; Tini];

// TIEMPO
tfin = 500; dt = 0.01; t = 0:dt:tfin; // min

// RESOLVER
x = ode(xini,0,t,f);
xfin = x(:,$)
dxdtfin = f(tfin,xfin)
Estacionario = abs(dxdtfin ./ xfin) < 1E-4

CA = x(1,:); CAee = CA($)
CB = x(2,:); CBee = CB($)
CC = x(3,:); CCee = CC($)
T = x(4,:); Tfin = T($)

// APARTADO B
// Puntos de corte
indexCBCC = find(CB<CC,1);
tCBCC = t(indexCBCC)
CBCC = CC(indexCBCC)

indexCACC = find(CA<CC,1);
tCACC = t(indexCACC)
CACC = CC(indexCACC)

// Gráfica de las concentraciones
scf(1);clf(1);
plot(t,CA,t,CB,t,CC);
plot(tCBCC,CBCC,'ro');
plot(tCACC,CACC,'ro');
xgrid;xtitle('t','CA(azul), CB(verde), CC(rojo)');

// APARTADO C
indexT=find(T<Tini,1);tT=t(indexT)
scf(2);clf(2);
plot(t,T,t(indexT),T(indexT),'ro');
xgrid; xtitle('t','T');

// APARTADO D
// Máximo global
[Tmax,indexTmax] = max(T)
tTmax = t(indexTmax)
plot(tTmax,Tmax,'ro');

// Mínimo global
[Tmin,indexTmin] = min(T)
tTmin = t(indexTmin)
plot(tTmin,Tmin,'ro');

Tmaxtest = Tmax < 325
Tmintest = Tmin > 315
scf(3);
if Tmaxtest & Tmintest then
    plot([TJ,TJ],[Tmin,Tmax],'go-');
else
    plot([TJ,TJ],[Tmin,Tmax],'ro-');
end
