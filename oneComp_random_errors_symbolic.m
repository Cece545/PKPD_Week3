% Define symbolic function
syms C(t)

%Random Errors from a Normal Distribution, adapted from oneComp_analytical.m
% parameters
k = 0.1;

% state the ODE
ode = diff(C, t) == -k * C;

% solve the ODE with initial condition
cond = C(0) == 100;
Csolve(t) = dsolve(ode, cond);

% create a time vector
time = 0.5:2:20.5;

% solve the symbolic solution at the time points
concentration = double(Csolve(time));

% create random errors
error_proportion = 0.05; % 5% proportionality
error_additive = 0.05 * mean(concentration); % 5% of average concentration for additive error

proportional_errors = error_proportion * randn(size(concentration));
additive_errors = error_additive * randn(size(concentration));

% add errors to the concentration data
noisy_concentration_proportional = concentration .* (1 + proportional_errors);
noisy_concentration_additive = concentration + additive_errors;

% Define initial guess values for C0 and ke
beta0 = [50, 1];

% fit the equations to noisy data (proportional error)
[parameters_proportional] = nlinfit(time, noisy_concentration_proportional, @conc, beta0);
disp('Parameter Estimates of C0 and ke (Proportional Error, 20.5h): ')
disp(parameters_proportional)

% fit the equations to noisy data (additive error)
[parameters_additive] = nlinfit(time, noisy_concentration_additive, @conc, beta0);
disp('Parameter Estimates of C0 and ke (Additive Error, 20.5h): ')
disp(parameters_additive)

% Plot noisy data and fits
figure;
plot(time, concentration, 'b-', 'LineWidth', 2, 'DisplayName', 'Symbolic Solution');
hold on;
plot(time, noisy_concentration_proportional, 'r--', 'LineWidth', 2, 'DisplayName', 'Proportional Error');
plot(time, noisy_concentration_additive, 'g-.', 'LineWidth', 2, 'DisplayName', 'Additive Error');
xlabel('Time (t)');
ylabel('Concentration (C)');
title('Noisy Data with Fitted Curves (Symbolic Solution)');
legend('Location', 'best');
grid on;
hold off;

% Function for modeling the fit of the data
function output = conc(c, t)
    C0 = c(1); % Finds value of C0
    k = c(2);  % Finds the rate constant of elimination
    output = C0 * exp(-k * t);
end
