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
omega1_vals = linspace(0.5, 2.0, 100); % Fewer omega1 values to reduce computation % Array of omega1 values
omega2_vals = linspace(0.05, 0.2, 100); % Fewer omega2 values to reduce computation % Array of omega2 values

% Transition probability storage
transition_probabilities_summary = zeros(length(omega1_vals), length(omega2_vals));

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

        % Initialize transition probability tracking
        transition_prob_final = 0;

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
        end

        % Calculate the final transition probability
        transition_prob_final = abs(sum(conj(psiH_2) .* psi_0))^2;
        transition_probabilities_summary(omega1_idx, omega2_idx) = transition_prob_final;
    end
end

% Plot the results as a 3D torus
[Theta1_grid, Theta2_grid] = ndgrid(linspace(0, 2*pi, length(omega1_vals)), linspace(0, 2*pi, length(omega2_vals)));

% Define torus coordinates with finer resolution
[Theta1_grid, Theta2_grid] = ndgrid(linspace(0, 2*pi, 50), linspace(0, 2*pi, 50)); % Increase resolution to get a full torus

% Define torus coordinates
R = lambda;
r = A;
X_torus = (R + r .* cos(Theta2_grid)) .* cos(Theta1_grid);
Y_torus = (R + r .* cos(Theta2_grid)) .* sin(Theta1_grid);
Z_torus = r .* sin(Theta2_grid);

% Interpolate the transition probabilities for a smoother torus surface
transition_prob_interp = interp2(...
    linspace(0, 2*pi, length(omega1_vals)), linspace(0, 2*pi, length(omega2_vals)), ...
    transition_probabilities_summary, ...
    Theta1_grid, Theta2_grid, 'linear');

% Plot the torus with a surface plot
figure;
surf(X_torus, Y_torus, Z_torus, transition_prob_interp, 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Transition Probability on 3D Torus');
colorbar;
caxis([0, 1]);
view(3);
axis equal;

