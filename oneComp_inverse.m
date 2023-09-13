function data_entry = Inverse_Problem()

clear all;
close all;
data_entry = 0;

% sampled euler data
time = 0.5:2:20.5;
concentration = [95.1786, 78.1082, 64.0994, 52.6031, 43.1687, 35.4263, 29.0726, 23.8584, 19.5794, 16.0678, 13.1860];

%beta0 is the guess vector with the initial guess values of C0 and ke
%we'll just guess 50 and 1 for the starting values

beta0 = [50 , 1];

% Generate a vector of the coefficient estimates using nlinfit
                                 
[parameters] = nlinfit(time, concentration, @conc, beta0);    
disp('Parameter Estimates of C0 and ke (20.5h): ')
disp(parameters)

data_entry = 1
return;

% function for modeling the fit of the data
function output = conc(c, t) 
      C0   = c(1); %finds value of C0
      k   = c(2); %finds the rate constant of elimination
      output  =  C0*exp(-k*t); 
return;