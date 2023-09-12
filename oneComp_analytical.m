% defines a symbolic function 
syms C(t)
k = 0.1;

% defines the ODE
ode = diff(C, t) == -k*C;

% solves the ODE with initial condition
cond = C(0) == 100;
Csolve(t) = dsolve(ode, cond);

fplot(Csolve)
xlabel('Time (hours)');
ylabel('Drug Concentration');
title('One-Compartment Pharmacokinetic Model');
grid on;
