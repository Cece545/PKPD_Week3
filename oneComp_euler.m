% number of steps n
n = 84;

% time start t0
t0 = 0;

% time end tf
tf = 20.5;

% initial condition C0
C0 = 100;

% elimination constant k
k=0.1;

% defines a symbolic function 
syms y(t)
k = 0.1;

% defines the ODE
ode = diff(y, t) == -k*y;

% solves the ODE with initial condition
cond = y(0) == 100;
ysolve(t) = dsolve(ode, cond);

% function dC/dt = f(t, C(t))
fun = @(t, C) -k*C;

% calling the euler function
C = euler(fun, t0, tf, C0, n);

% vector C is n x 2
% first column: linspace for time
% second column: values after each step

% plotting
plot(C(:,1),C(:,2), 'linewidth', 2)
title(['Euler Method With ',num2str(n),' steps'])
xlabel('Time (Hours)')
ylabel('Concentration')

errors = zeros(n,1);
for i=1:n
    % subtract euler value from analytical solution value
    errors(i) = C(i,2) - ysolve(C(i,1));
end

errors = errors'
percent_error = transpose((errors./C(i,2)).*100)
percent_error = percent_error'

% generate samples for fitting
time = 0.5:2:20.5;
samples = zeros(11,1);
for i=1:11
    samples(i) = C(8*i-5,2);
end

euler_samples = [time(:), samples(:)]
time = transpose(time)
samples = transpose(samples)

function output = euler(fun, t0, tf, C0, n)
    % define the size of the step
    h = (tf-t0)/n;

    % defining time vector and empty y-vector
    t = (t0:h:tf)';
    C = zeros(n+1, 1);

    C(1) = C0;

    for i=1:n
       C(i+1)=C(i)+h*fun(t(i),C(i));
    end
    output = [t,C]
end