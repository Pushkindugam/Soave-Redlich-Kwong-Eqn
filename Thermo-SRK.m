% Given data
R=0.08314; % (we could also use 0.0821)
T = 400; % Temperature in Kelvin
P = 8.234; % Pressure in bar

% SRK constants for propane
a = 9.385; % Atm·L^2/mol^2
b = 0.09044; % L/mol
alpha = 0.944; % dimensionless

% Initial guesses for molar volume of vapor and liquid
V_vapor = 1 * R * T / P; % Initial guess for molar volume of vapor
V_liquid = 0.5 * R * T / P; % Initial guess for molar volume of liquid

% Tolerance for convergence
tolerance = 1e-6;

% Maximum number of iterations
maxIterations = 100;

% SRK equation and its derivative with respect to molar volume
SRK = @(V) P - R * T / (V - b) + a * alpha / (V * (V + b));
dSRK_dV = @(V) R * T / (V - b)^2 - a * alpha / (V + b)^2 - 2 * a * alpha / (V * (V + b)^2);

% Newton-Raphson iteration
for iteration = 1:maxIterations
    % Calculate the function value and its derivative at the current guess
    f = SRK(V_vapor);
    df_dV = dSRK_dV(V_vapor);
    
    % Update the guess for the molar volume using Newton-Raphson formula
    V_vapor = V_vapor - f / df_dV;
    
    fprintf('Molar volume of saturated-vapor propane at 400 K after %d iteration: %.4f L/mol\n', iteration, V_vapor);

    % Check for convergence
    if abs(f) < tolerance
        break;
    end
end

% Display results
fprintf('Molar volume of saturated-vapor propane at 400 K after %d iterations: %.4f L/mol\n', iteration, V_vapor);

% Repeat the process for saturated-liquid propane
for iteration = 1:maxIterations
    % Calculate the function value and its derivative at the current guess
    f = SRK(V_liquid);
    df_dV = dSRK_dV(V_liquid);
    
    % Update the guess for the molar volume using Newton-Raphson formula
    V_liquid = V_liquid - f / df_dV;
    fprintf('Molar volume of saturated-Liquid propane at 400 K after %d iteration: %.4f L/mol\n', iteration, V_liquid);

    % Check for convergence
    if abs(f) < tolerance
        break;
    end
end

% Display results
fprintf('Molar volume of saturated-liquid propane at 400 K after %d iteration: %.4f L/mol\n',iteration, V_liquid)