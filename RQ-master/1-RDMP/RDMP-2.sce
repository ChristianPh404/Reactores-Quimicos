// clear; clc; 
// // RDMP-2.sce
// // A => B
// // Adiabático

// // SISTEMA DE ECUACIONES DIFERENCIALES
// function dxdt = f(t,x)
//     CA = x(1)
//     T  = x(2)
//     // Ecuación de Arrhenius
//     k = k0*exp(-E/(R*T))
//     // Velocidad de reacción
//     r = k*CA
//     // Balance de materia para A
//     // d(V*CA)dt = -r*V consume a por eso es -r y el reacto de a es r lo consumido *v
//     //v es constante por lo que puede salir de la derivada
//     dCAdt = -r  
//     // Variables diferenciales como ya tenemos CA en una ec definida la ponemos como constante 
//        //CA = x(1) que pasa si ponemos la ec. aqui comprobar debido a que tenemos ya
//        //la variable aunque no conozcamos su valor similar abajo para T         
//     // Balance de energía
//     // d(V*RHO*CP*T)dt  //<····(m·cp·t)//= -H*r*V ; el signo menos de la entalpia es por que al ser 
//     //exotermica por convenio la h es negativa y como tiene que aumentar 
//     //la temperatura es necesario para -·-=*
//     dTdt  = -H*r/(RHO*CP)   
   
//     // Derivadas
//     dxdt(1) = dCAdt
//     dxdt(2) = dTdt
//     //si la reaccion es irreversible no se suele poner los balances de materia del producto 
//     //se puede hacer aunque no sea necesario ya que la potencia de calculo no la complica
// endfunction

// // CONSTANTES
// //en los ejercicios nunca vamos  a tener que realizar cambios de unidades pero si poner las unidades
// //constante que buscar siempre la r es decir tenerlas tabuladas
// CP = 0.9; // cal/(g*K)
// RHO = 1070; // g/L
// H = -50400; // cal/mol
// k0 = 4.15E5; // s-1
// E = 11200; // cal/mol
// R = 1.987; // cal/(mol*K)

// // CONDICIONES INICIALES
// CAini = 0.5; // mol/L
// Tini = 285; // K
// xini = [CAini; Tini];

// // TIEMPO
// tfin = 1000; dt = 1; t = 0:dt:tfin;// s
// //el tiempo no nos lo da el ejercicio ponemos un numero provisional de prueba y error

// // RESOLVER
// x = ode(xini,0,t,f);
// CA = x(1,:); CAfin = CA($)
// T = x(2,:); Tfin = T($)

// XA = 1 - CA/CAini; XAfin = XA($) //vector de conversion por definicion
// XAobj = 0.90;
// indexXAobj = find(XA>XAobj,1); //index para encontrar la posicion es decir a 
// tXAobj = t(indexXAobj)//que temperatura se alcanza el 90% de conversion
// TXAobj = T(indexXAobj)

// // GRÁFICAS
// scf(1); clf(1); 
// plot(t,CA);
// xgrid; xlabel('t'); legend('CA',-2,%f);

// scf(2); clf(2); 
// plot(t,XA,'m-',tXAobj,XAobj,'mo');
// xgrid; xlabel('t'); legend('XA',-2,%f);

// scf(3); clf(3); 
// plot(t,T,'r-',tXAobj,TXAobj,'ro');
// xgrid; xlabel('t'); legend('T',-2,%f);// el ,-2 es la posicion arriba izquiera fuera
// // el -3 es abajo izquierda fuera del grafico , el %f es para eliminar la caja 
// //de la leyenda

clear; clc;
// RDMP-2.sce
// A => B
// Adiabático

// SISTEMA DE ECUACIONES DIFERENCIALES
function dxdt = f(t, x)
    CA = x(1);
    T  = x(2);
    
    // Ecuación de Arrhenius
    k = k0 * exp(-E / (R * T));
    
    // Velocidad de reacción
    r = k * CA;
    
    // Balance de materia para A
    dCAdt = -r;
    
    // Balance de energía
    dTdt = -H * r / (RHO * CP);
    
    // Retornar como columna
    dxdt = [dCAdt; dTdt]; 
endfunction

// CONSTANTES
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

// TIEMPO INICIAL
tfin = 50; dt = 1;
XAobj = 0.90; // Conversión objetivo
XAalcanzado = 0; // Variable para comprobar si se alcanza XAobj
tfinal_max = 10000; // Tiempo máximo permitido
t = 0:dt:tfin; 

// SOLVER
while (XAalcanzado < XAobj & tfin < tfinal_max)
    x = ode(xini, 0, t, f);
    
    // Extraer valores de concentración y temperatura
    CA = x(1, :);
    T = x(2, :);
    
    // Calcular conversión
    XA = 1 - CA / CAini;
    
    // Revisar si se alcanzó el valor objetivo
    if max(XA) >= XAobj then
        XAalcanzado = max(XA);
        break; // Termina el bucle si se alcanza la conversión objetivo
    else
        // Aumentar el tiempo de simulación si no se alcanza XAobj
        tfin = tfin + 50;
        t = 0:dt:tfin;
    end
end

// Encontrar el tiempo en que se alcanza XAobj
indexXAobj = find(XA >= XAobj, 1);
if indexXAobj <> [] then
    tXAobj = t(indexXAobj);
    TXAobj = T(indexXAobj);
else
    tXAobj = NaN;
    TXAobj = NaN;
end

// GRÁFICAS
scf(1); clf(1); 
plot(t, CA, 'b');
xgrid; xlabel('t'); legend('CA', -2, %f);

scf(2); clf(2); 
plot(t, XA, 'm-', tXAobj, XAobj, 'mo');
xgrid; xlabel('t'); legend('XA', -2, %f);

scf(3); clf(3); 
plot(t, T, 'r-', tXAobj, TXAobj, 'ro');
xgrid; xlabel('t'); legend('T', -2, %f);
disp('la conversion de A es de '+string(XAalcanzado*100)+'%');
disp('la temperatura en la que se alcanza el 90% de conversion es de '+string(TXAobj)+' K');
disp('el tiempo en el que se alcanza el 90% de conversion es de '+string(tXAobj)+' s');
disp('el tiempo maximo permitido es de '+string(tfinal_max)+' s');
