% Define the one-compartment model as an inline function

% Define parameters
k = 0.1; % Elimination rate constant (you can adjust this value)

% Define the time span for the simulation
tspan = [0 10]; % Start at time 0 and end at time 10 (you can adjust this)

% Initial drug concentration
C0 = 100; % Initial concentration (you can adjust this value)

% Solve the differential equation using ode45 and the oneCompModel function
[t,C] = ode45(@(t,C) (-k)*C,  tspan, C0);

% Plot the drug concentration vs. time
plot(t, C, '-o');
xlabel('Time (hours)');
ylabel('Drug Concentration');
title('One-Compartment Pharmacokinetic Model');
grid on;

