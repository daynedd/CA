%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bifurcation Diagram with Edge of Chaos Analysis for a Wavepacket in a 1D 
% Harmonic Trap with a Time-Dependent Potential, with Frequency ω as a Control Parameter
% Over Multiple Loops to Capture Long-Term Behavior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for space and simulation
a = -20;                       % Left end point 
b = +20;                       % Right end point 
L = b - a;                     % Width of the space
N = 512;                       % No. of cells
X = a + L * (0:N-1) / N;       % Dimensionless coordinates
P = (2 * pi / L) * [0:N/2-1, -N/2:-1]; % Dimensionless momentum
T = 5 * pi;                    % Time duration of one evolution loop
M = 10^3;                      % Total No. of steps in the evolution for each loop
dt = T / M;                    % Time step
A = 0.5;                       % Amplitude for the sin term
B = 0.3;                       % Amplitude for the cos term in torus-like potential

% Number of loops (repetitions) to capture long-term behavior
num_loops = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the initial state
% Initial state: Gaussian wavepacket located at X0
X0 = 4.0;                             % Initial position of the wave packet
sigma = 1.0;                          % Width of the initial wavepacket
psiprep = exp(-(X - X0).^2 / (2 * sigma^2));  % Gaussian state
psi = psiprep / sqrt(sum(abs(psiprep).^2));   % Normalized state

% Slight perturbation for Lyapunov analysis
initial_psi_perturbed = psi * (1 + 1e-5); % Slight perturbation to the initial state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bifurcation Analysis with Frequency as Parameter
% Parameter range for bifurcation analysis (ω1)
omega1_params = linspace(0.01, 0.2, 50);  % Frequency ω1 range
bifurcation_results = [];                 % Store bifurcation results
lyapunov_exponents = [];                  % Store Lyapunov exponents

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Evolution Loop for Bifurcation Analysis
for param_idx = 1:length(omega1_params)
    % Update control parameter for bifurcation analysis
    omega1 = omega1_params(param_idx);

    % Reset wavefunction for each parameter
    psi_0 = psi;
    psi_perturbed = initial_psi_perturbed;
    initial_distance = norm(psi_perturbed - psi_0);

    % Evolution with current frequency ω1 over multiple loops
    lyapunov_sum = 0;  % To accumulate distance growth over all loops
    for loop = 1:num_loops
        for m = 1:M
            t = (m + (loop - 1) * M) * dt;  % Current time across loops
            
            % Define the torus-like parameters theta1 and theta2
            theta1 = omega1 * t;  
            theta2 = 0.02 * t;  % Keep omega2 constant
            
            % Update the time-dependent potential with both periodic components
            V_t = A * sin(X) .* cos(theta1) + B * cos(X) .* sin(theta2);
            UV_t = exp(-1i * (X.^2 / 2 + V_t) * dt / 2);  % Time-dependent potential propagator
            
            % Split-step evolution for the original wavepacket
            psi_1 = UV_t .* psi_0;          
            phi_2 = fft(psi_1);              
            phi_3 = exp(-1i * (P.^2 / 2) * dt) .* phi_2;  
            psi_3 = ifft(phi_3);             
            psi_4 = UV_t .* psi_3;           
            psi_0 = psi_4;  % Update for next cycle

            % Split-step evolution for the perturbed wavepacket
            psi_1_pert = UV_t .* psi_perturbed;
            phi_2_pert = fft(psi_1_pert);
            phi_3_pert = exp(-1i * (P.^2 / 2) * dt) .* phi_2_pert;
            psi_3_pert = ifft(phi_3_pert);
            psi_4_pert = UV_t .* psi_3_pert;
            psi_perturbed = psi_4_pert;

            % Calculate and accumulate Lyapunov distance growth
            current_distance = norm(psi_perturbed - psi_0);
            lyapunov_sum = lyapunov_sum + log(current_distance / initial_distance);
        end

        % Store bifurcation results (use last state of each loop)
        max_prob_density = max(abs(psi_0).^2);
        bifurcation_results = [bifurcation_results; omega1, max_prob_density];
    end

    % Calculate the average Lyapunov exponent over all loops
    lyapunov_exponent = lyapunov_sum / (num_loops * M * dt);
    lyapunov_exponents = [lyapunov_exponents; omega1, lyapunov_exponent];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify Edge of Chaos
% Define a threshold for identifying the edge of chaos
lyapunov_threshold = 0.02; % Close to zero to identify edge of chaos
edge_of_chaos_params = lyapunov_exponents(abs(lyapunov_exponents(:,2)) < lyapunov_threshold, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Bifurcation Diagram
figure;
plot(bifurcation_results(:,1), bifurcation_results(:,2), 'k.');
hold on;
% Highlight edge of chaos regions
if ~isempty(edge_of_chaos_params)
    scatter(edge_of_chaos_params(:,1), ones(size(edge_of_chaos_params, 1), 1) * max(bifurcation_results(:,2)), ...
        50, 'r', 'filled', 'DisplayName', 'Edge of Chaos');
end
xlabel('Frequency \omega_1');
ylabel('Maximum Probability Density');
title('Bifurcation Diagram with Frequency \omega_1 as Parameter (Multiple Loops)');
legend('show');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Lyapunov Exponents
figure;
plot(lyapunov_exponents(:,1), lyapunov_exponents(:,2), 'b-o');
hold on;
% Highlight edge of chaos regions
if ~isempty(edge_of_chaos_params)
    scatter(edge_of_chaos_params(:,1), edge_of_chaos_params(:,2), 50, 'r', 'filled', 'DisplayName', 'Edge of Chaos');
end
xlabel('Frequency \omega_1');
ylabel('Lyapunov Exponent');
title('Lyapunov Exponent Analysis with Frequency \omega_1 as Parameter');
legend('show');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
