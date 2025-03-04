clear; clc; 
// RDMP-2.sce
// A => B
// Adiabático

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = f(t,x)
    CA = x(1)
    T  = x(2)
    // Ecuación de Arrhenius
    k = k0*exp(-E/(R*T))
    // Velocidad de reacción
    r = k*CA
    // Balance de materia para A
    // d(V*CA)dt = -r*V consume a por eso es -r y el reacto de a es r lo consumido *v
    //v es constante por lo que puede salir de la derivada
    dCAdt = -r  
    // Variables diferenciales como ya tenemos CA en una ec definida la ponemos como constante 
       //CA = x(1) que pasa si ponemos la ec. aqui comprobar debido a que tenemos ya
       //la variable aunque no conozcamos su valor similar abajo para T         
    // Balance de energía
    // d(V*RHO*CP*T)dt  //<····(m·cp·t)//= -H*r*V ; el signo menos de la entalpia es por que al ser 
    //exotermica por convenio la h es negativa y como tiene que aumentar 
    //la temperatura es necesario para -·-=*
    dTdt  = -H*r/(RHO*CP)   
   
    // Derivadas
    dxdt(1) = dCAdt
    dxdt(2) = dTdt
    //si la reaccion es irreversible no se suele poner los balances de materia del producto 
    //se puede hacer aunque no sea necesario ya que la potencia de calculo no la complica
endfunction

// CONSTANTES
//en los ejercicios nunca vamos  a tener que realizar cambios de unidades pero si poner las unidades
//constante que buscar siempre la r es decir tenerlas tabuladas
CP = 0.9; // cal/(g*K)
RHO = 1070; // g/L
H = -50400; // cal/mol  //* Reacción exotérmica ya que la entalpia de reaccion tiene signo negativo
k0 = 4.15E5; // s-1
E = 11200; // cal/mol
R = 1.987; // cal/(mol*K)

// CONDICIONES INICIALES
CAini = 0.5; // mol/L
Tini = 285; // K
xini = [CAini; Tini];

// TIEMPO
tfin = 1000; dt = 1; t = 0:dt:tfin;// s
//el tiempo no nos lo da el ejercicio ponemos un numero provisional de prueba y error

// RESOLVER
x = ode(xini,0,t,f);
CA = x(1,:); CAfin = CA($)
T = x(2,:); Tfin = T($)

XA = 1 - CA/CAini; XAfin = XA($) //vector de conversion por definicion
XAobj = 0.90;
indexXAobj = find(XA>XAobj,1); //primera vez que se cumple la condición
tXAobj = t(indexXAobj)//que temperatura se alcanza el 90% de conversion
TXAobj = T(indexXAobj)

// GRÁFICAS
scf(1); clf(1); 
plot(t,CA);
xgrid; xlabel('t'); legend('CA',-2,%f);

scf(2); clf(2); 
plot(t,XA,'m-',tXAobj,XAobj,'mo');
xgrid; xlabel('t'); legend('XA',-2,%f);

scf(3); clf(3); 
plot(t,T,'r-',tXAobj,TXAobj,'ro');
xgrid; xlabel('t'); legend('T',-2,%f);// el ,-2 es la posicion arriba izquiera fuera
// el -3 es abajo izquierda fuera del grafico , el %f es para eliminar la caja 
//de la leyenda
disp('Concentracion final de A: '+string(CAfin)+' mol/L');
disp('Temperatura final: '+string(Tfin)+' K');
disp('Conversión final de A: '+string(XAfin)+' mol/L');
disp('Tiempo para alcanzar el 90% de conversión: '+string(tXAobj)+' s');