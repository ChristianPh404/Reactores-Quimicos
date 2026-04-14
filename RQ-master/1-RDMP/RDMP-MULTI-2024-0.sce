clear; clc;
// 2024-06-10-A.sce

// RDMP-MULT
// 1) 2 A => B
// 2)   A => C 
// No adiabático

function dxdt=f(t, x) 
//Variables diferenciales 
CA = x(1) 
CB = x(2) 
CC = x(3) 
T = x(4) 
//Ecuación de Arrhenius 
k1 = 1.2E16*exp(-14500/T) 
k2 = 9.4E12*exp(-10600/T) 
//Coeficiente transmisión de calor  
//1/UA = 1/hi  
//1/UA = 1/(m*CP*(T-TJ)) 
//1/UA = 1/(RHO*F*CP*(T-TJ)) 
//1/UA = 1/(RHO*V*CA0*CP*(T-TJ)) 
UA = RHO*CA0*CP*(T-TJ); 
//Ecuación de velocidad de reacción 
r1 = k1*CA^2 
r2 = k2*CA 
//Calor transferido 
Q= UA*(T-TJ) 
//Balance de materia para A  
//dCAdt = -rA*V 
dCAdt = -r1-r2 
//Balance de materia para B  
//dCBdt = rB*V 
dCBdt = r1 
//Balance de materia para C  
//dCCdt = rC*V 
dCCdt = r2 
//Balance de energía 
//d(M*CP*T)dt = -H1*r1*V - H2*r2*V -Q/(RHO*CP*V) 
dTdt = -H1*r1/(RHO*CP) - H2*r2/(RHO*CP) -Q/(RHO*CP*V) 
//Derivadas 
dxdt (1) = dCAdt  
dxdt (2) = dCBdt  
dxdt (3) = dCCdt  
dxdt (4) = dTdt     
endfunction 
//Constantes 
H1 = -29000; // cal/mol 
H2 = -23000; // cal/mol 
RHO = 1000; //g/L 
CP = 1 ; //cal/(g*K) 
V = 0.1; // L  
TJ = 320; //K 
CA0 = 5; //mol/L 
//Tiempo 
tfin = 50; dt = 0.1; t= 0:dt:tfin; //h 
//Condiciones iniciales 
CAini = 5; CBini = 0; CCini = 0; //mol/L 
Tini = 320; //K 
xini = [CAini; CBini; CCini; Tini]; 
//Resolver 
x = ode(xini, 0, t, f); 
CA = x(1,:); CAfin = CA ($); 
disp('La concentracion final de A es: ' + string(CAfin) + ' mol/L');
CB = x(2,:); CBfin = CB($); 
disp('La concentracion final de B es: ' + string(CBfin) + ' mol/L');
CC = x(3,:); CCfin = CC ($); 
disp('La concentracion final de C es: ' + string(CCfin) + ' mol/L');
T = x(4,:); Tfin = T ($); 
disp('La temperatura final es: ' + string(Tfin) + ' K');
XA= 1-CA/CAini; XAfin = XA($); 
disp('La conversion final de A es: ' + string(XAfin) + ' mol/L');
XAobj = 0.80; 
indexXAobj = find(XA>XAobj,1); 
tXAobj = t(indexXAobj); 

TXAobj = T(indexXAobj);
disp('La conversion objetivo de A es: ' + string(XAobj) + ' mol/L' +'que se consigue a los ' + string(tXAobj) + ' h y a una temperatura de ' + string(TXAobj) + ' K');
//Apartado B)  
[Tmax, indexTmax] = max(T) ;
tTmax = t(indexTmax) ;
disp('La temperatura maxima es: ' + string(Tmax) + ' K' +' que se alcanza a los ' + string(tTmax) + ' h');
//Gráficas 
scf(1); clf(1) 
plot(t,CA,t,CB,t,CC) 
xgrid; xlabel('t');legend('CA','CB','CC',-2,%f); 
scf(2); clf(2) 
plot(t,T,'r-',tTmax,Tmax,'ro') 
xgrid; xlabel('t');legend('T',-2,%f); 
scf(3);clf(3) 
plot(T,XA,'m-',TXAobj,XAobj,'mo')
xgrid; xlabel('T');legend('XA',-2,%f); 
scf(4);clf(4) 
plot(t,XA,tXAobj,XAobj,'o') 
xgrid; xlabel('t');legend('XA',-2,%f);