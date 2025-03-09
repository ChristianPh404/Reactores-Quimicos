// Parámetros adicionales
C_A0 = 1; // Concentración inicial del reactivo [mol/L]
k0 = 0.01; // Constante de velocidad a referencia [1/s]
Ea = 5000; // Energía de activación [J/mol]
R = 8.314; // Constante de los gases [J/mol*K]

// Ecuación diferencial de balance de materia + energía
function dXdt = reactorSistema(X, t, Tj)
    T = X(1);
    C_A = X(2);
    
    // Ecuación de Arrhenius para k(T)
    k = k0 * exp(-Ea / (R * T));
    
    // Balance de energía
    dTdt = UA * (Tj - T) / (rho * Cp * V);
    
    // Balance de materia
    dCAdt = -k * C_A;
    
    dXdt = [dTdt; dCAdt];
endfunction

// Simulación en los 4 casos
C_A_finales = [];

for i = 1:4
    select i
        case 1 then
            Tj = Tj1; // Siempre a 293K
        case 2 then
            Tj = Tj2; // Siempre a 353K
        case 3 then
            Tj = ifelse(t < t_total/2, Tj1, Tj2); // Cambio 293K -> 353K
        case 4 then
            Tj = ifelse(t < t_total/2, Tj2, Tj1); // Cambio 353K -> 293K
    end
    
    // Condiciones iniciales
    X0 = [T0; C_A0];
    
    // Resolver la EDO
    X_sol = ode(X0, t, list(reactorSistema, Tj));
    
    // Extraer la concentración de A y la temperatura
    T_sol = X_sol(:,1);
    C_A_sol = X_sol(:,2);
    
    // Guardar el valor final de la concentración de A
    C_A_finales = [C_A_finales; C_A_sol($)];
end

// Encontrar el caso con menor concentración final de A
[min_CA, idx] = min(C_A_finales);

// Mostrar resultado
disp("El caso con menor concentración final de A es el caso: " + string(idx));
disp("Concentración final mínima de A: " + string(min_CA));
