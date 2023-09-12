% number of steps n
n = 10;

% time start t0
t0 = 0.5;

% time end tf
tf = 20.5;

% initial condition C0
C0 = 100;

% elimination constant k
k=0.1;

fun = @(t, C) -k*C;

euler_data = euler(fun, t0, tf, C0, n)

time = euler_data(:,1);
concentration = euler_data(:,2);



function euler_output = euler(fun, t0, tf, C0, n)
    % define the size of the step
    h = (tf-t0)/n;

    % defining time vector and empty y-vector
    t = (t0:h:tf)';
    C = zeros(n+1, 1);

    C(1) = C0;

    for i=1:n
       C(i+1)=C(i)+h*fun(t(i),C(i));
    end
    euler_output = [t,C];
end