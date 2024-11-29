clear;
a = -20;                       % Left end point
b = +20;                       % Right end point
L = b-a;                        % Width of the space
N = 512;                       % No. of cells
X = a+L*(0:N-1)/N;                % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1]; % Dimensionless momentum
T = 10*pi;                     % Time duration of the evolution
M = 40^2; % Reduced number of steps in the evolution                      % Total No. of steps in the evolution
dt = T/M;                      % Time step
A = 0.1;                       % Constant A value
lambda = 0.05;                 % Constant lambda value
omega1_vals = linspace(0.5, 2, 20); % Increased omega1 values to enhance resolution
omega2_vals = linspace(0.01, 0.3, 20); % Increased omega2 values to enhance resolution

% Transition probability storage for maximum probabilities over time
max_transition_probabilities = zeros(length(omega1_vals), length(omega2_vals));
max_prob = 0;
max_theta1_idx = 0;
max_theta2_idx = 0;

% Loop over different parameter values to create the 3D torus plot data
for omega1_idx = 1:length(omega1_vals)
    omega1 = omega1_vals(omega1_idx);
    for omega2_idx = 1:length(omega2_vals)
        omega2 = omega2_vals(omega2_idx);
        
        % Initialize wavefunction and parameters
        n = 0;
        Hn = hermiteH(n, X);
        psi_1 = (1 / sqrt(2^n * factorial(n) * sqrt(pi))) * Hn .* exp(-X.^2 / 2);
        psi = psi_1 / sqrt(sum(abs(psi_1).^2)); % normalized state

        n2 = 1;
        Hn2 = hermiteH(n2, X);
        psiH_2 = (1 / sqrt(2^n2 * factorial(n2) * sqrt(pi))) * Hn2 .* exp(-X.^2 / 2);
        psiH_2 = psiH_2 / sqrt(sum(abs(psiH_2).^2));

        natural_frequency = (n2 - n) * 1;

        % Split-step propagator definitions
        UV = exp(-1i * (X.^2 / 2 + lambda * X.^4) * dt / 2);
        UT = exp(-1i * (P.^2 / 2) * dt);
        psi_0 = psi; % Initialize state for evolution

        % Initialize transition probability tracking for each time step
        current_max_prob = 0;

        for m = 1:M
            t = m * dt;
            V_t = A .* sin(X) * (cos(omega1 * t) + cos(omega2 * t)) + lambda * X.^4;
            UV_t = exp(-1i * (X.^2 / 2 + V_t) * dt / 2);

            % Split-operator propagation
            psi_1 = UV_t .* psi_0;
            phi_2 = fft(psi_1);
            phi_3 = UT .* phi_2;
            psi_3 = ifft(phi_3);
            psi_4 = UV_t .* psi_3;
            psi_0 = psi_4;

            % Calculate the transition probability at this time step
            transition_prob = abs(sum(conj(psiH_2) .* psi_0))^2;
            if transition_prob > current_max_prob
                current_max_prob = transition_prob;
            end
        end

        % Store the maximum transition probability for this combination of omega1 and omega2
        max_transition_probabilities(omega1_idx, omega2_idx) = current_max_prob;

        % Update the global maximum probability and corresponding indices
        if current_max_prob > max_prob
            max_prob = current_max_prob;
            max_theta1_idx = omega1_idx;
            max_theta2_idx = omega2_idx;
        end
    end
end

% Retrieve the corresponding Theta1 and Theta2 values for the maximum transition probability
max_theta1 = omega1_vals(max_theta1_idx);
max_theta2 = omega2_vals(max_theta2_idx);

fprintf('Maximum transition probability is %f at:\n', max_prob);
fprintf('Theta1: %f radians\n', max_theta1);
fprintf('Theta2: %f radians\n', max_theta2);

% Plot the maximum transition probability as a color map on the torus
[Theta1_grid, Theta2_grid] = ndgrid(linspace(0.5, 2, length(omega1_vals)), linspace(0, 2*pi, length(omega2_vals)));

% Define torus coordinates
R = lambda;
r = A;
X_torus = (R + r .* cos(Theta2_grid)) .* cos(Theta1_grid);
Y_torus = (R + r .* cos(Theta2_grid)) .* sin(Theta1_grid);
Z_torus = r .* sin(Theta2_grid);

% Plot the torus with the colormap representing the maximum transition probability
figure;
surf(X_torus, Y_torus, Z_torus, max_transition_probabilities, 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
zlabel('Z');
title({'Maximum Transition Probability on 3D Torus', 'Maximum at:', sprintf('Theta1: %f, Theta2: %f', max_theta1, max_theta2)});
colorbar;
caxis([0, 1]);
view(3);
axis equal;
grid on;