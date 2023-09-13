function data_entry = Inverse_Problem()

    clear all;
    close all;
    data_entry = 0;

    % Sampled time points
    time = 0.5:2:20.5;
    
    % Original concentration data (decreasing)
    concentration = [95.1786, 78.1082, 64.0994, 52.6031, 43.1687, 35.4263, 29.0726, 23.8584, 19.5794, 16.0678, 13.1860];

    % Generate random errors
    error_proportion = 0.05; % 5% proportionality
    error_additive = 0.05 * mean(concentration); % 5% of average concentration for additive error

    proportional_errors = error_proportion * randn(size(concentration));
    additive_errors = error_additive * randn(size(concentration));

    % Add errors to the concentration data
    noisy_concentration_proportional = concentration .* (1 + proportional_errors);
    noisy_concentration_additive = concentration + additive_errors;

    % Define initial guess values for C0 and ke
    beta0 = [50, 1];

    % Fit equations to noisy data (proportional error)
    [parameters_proportional] = nlinfit(time, noisy_concentration_proportional, @conc, beta0);
    disp('Parameter Estimates of C0 and ke (Proportional Error, 20.5h): ')
    disp(parameters_proportional)

    % Fit equations to noisy data (additive error)
    [parameters_additive] = nlinfit(time, noisy_concentration_additive, @conc, beta0);
    disp('Parameter Estimates of C0 and ke (Additive Error, 20.5h): ')
    disp(parameters_additive)

    data_entry = 1;

    % Plot noisy data and fits
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

    % show the difference between IV and oral administration
    % ex 1: IV Administration (no absorption)
    C_IV = concentration; % IV administration data (no errors)
    
    % ex 2: Oral Administration with Absorption (ka = 0.2)
    Dose = 200; % Initial dose for oral administration
    ka = 0.2; % Absorption rate constant for oral administration
    C_oral_1 = OralAdministration(time, Dose, ka);
    
    % ex 3: Oral Administration with Slower Absorption (ka = 0.05)
    ka = 0.05; % Slower absorption rate constant for oral administration
    C_oral_2 = OralAdministration(time, Dose, ka);

    % Plot the examples
    figure;
    plot(time, C_IV, 'b-', 'LineWidth', 2, 'DisplayName', 'IV Administration');
    hold on;
    plot(time, C_oral_1, 'r--', 'LineWidth', 2, 'DisplayName', 'Oral Administration (ka = 0.2)');
    plot(time, C_oral_2, 'g-.', 'LineWidth', 2, 'DisplayName', 'Oral Administration (ka = 0.05)');
    xlabel('Time (t)');
    ylabel('Concentration (C)');
    title('Comparison of IV and Oral Administration');
    legend('Location', 'best');
    grid on;
    hold off;
    
    return;

    % Function for modeling the fit of the data
    function output = conc(c, t)
        C0 = c(1); % finds value of C0
        k = c(2);  % finds rate constant of elimination
        output = C0 * exp(-k * t);
    end

    % Function to calculate drug concentration after oral administration
    function C_oral = OralAdministration(time, Dose, ka)
        k = 0.1; % Elimination rate constant
        Vd = 100; % Volume of distribution
        C_oral = zeros(size(time)); % Initialize concentration vector
        
        for i = 1:length(time)
            t = time(i);
            if t <= 0
                C_oral(i) = 0; % No drug before administration
            else
                % Modified oral administration model to decrease concentration over time
                C_oral = (((Dose/Vd) * ka) / ((k-ka))) * (1 - exp(-k * time)) + exp(k * time);
            end
        end
    end
end

