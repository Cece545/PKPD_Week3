%Random Errors from a Normal Distribution, adapted from oneComp_inverse.m
function data_entry = Inverse_Problem()

    clear all;
    close all;
    data_entry = 0;

    % Sampled euler data
    time = 0.5:2:20.5;
    concentration = [95.1786, 78.1082, 64.0994, 52.6031, 43.1687, 35.4263, 29.0726, 23.8584, 19.5794, 16.0678, 13.1860];

    % random proportional and additive errors
    % 5% proportionality
    error_proportion = 0.05; 
    % 5% of average concentration for additive error
    error_additive = 0.05 * mean(concentration); 
    % create random errors  
    proportional_errors = error_proportion * randn(size(concentration));
    additive_errors = error_additive * randn(size(concentration));

    % add errors to the concentration data
    noisy_concentration_proportional = concentration .* (1 + proportional_errors);
    noisy_concentration_additive = concentration + additive_errors;

    % Define initial guess values for C0 and ke
    beta0 = [50, 1];

    % fit equations to noisy data (proportional error)
    [parameters_proportional] = nlinfit(time, noisy_concentration_proportional, @conc, beta0);
    disp('Parameter Estimates of C0 and ke (Proportional Error, 20.5h): ')
    disp(parameters_proportional)

    % fit equations to noisy data (additive error)
    [parameters_additive] = nlinfit(time, noisy_concentration_additive, @conc, beta0);
    disp('Parameter Estimates of C0 and ke (Additive Error, 20.5h): ')
    disp(parameters_additive)

    data_entry = 1;

    % plot noisy data and fits
    figure;
    plot(time, concentration, 'b-', 'LineWidth', 2, 'DisplayName', 'Original Data');
    hold on;
    plot(time, noisy_concentration_proportional, 'r--', 'LineWidth', 2, 'DisplayName', 'Proportional Error');
    plot(time, noisy_concentration_additive, 'g-.', 'LineWidth', 2, 'DisplayName', 'Additive Error');
    xlabel('Time (t)');
    ylabel('Concentration (C)');
    title('Noisy Data with Fitted Curves');
    legend('Location', 'best');
    grid on;
    hold off;
    
    return;

    % Ffnction for modeling the fit of the data
    function output = conc(c, t)
        C0 = c(1); % finds value of C0
        k = c(2);  % finds rate constant of elimination
        output = C0 * exp(-k * t);
    end
end