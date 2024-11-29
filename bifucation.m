clear;
a = -20;                       % Left end point
b = +20;                       % Right end point
L = b - a;                     % Width of the space
N = 512;                       % No. of cells
X = a + L * (0:N-1) / N;       % Dimensionless coordinates
P = (2 * pi / L) * [0:N/2-1, -N/2:-1]; % Dimensionless momentum
T = 10 * pi;                   % Time duration of the evolution
M = 40^2;                      % Reduced number of steps in the evolution
dt = T / M;                    % Time step
A = 0.1;                       % Constant A value
omega1 = 1.0526;               % Fixed value from previous analysis (resonance)
omega2 = 0.2084;               % Fixed value from previous analysis (resonance)

% Range of lambda values for bifurcation analysis
lambdas = linspace(0, 1, 100);  % Lambda values to be explored
max_transition_probabilities = zeros(length(lambdas), 1);  % To store max transition probability

% Loop over different lambda values to create the bifurcation plot data
for lambda_idx = 1:length(lambdas)
    lambda = lambdas(lambda_idx);
    
    % Initialize wavefunction and parameters
    n = 0;
    Hn = hermiteH(n, X);
    psi_1 = (1 / sqrt(2^n * factorial(n) * sqrt(pi))) * Hn .* exp(-X.^2 / 2);
    psi = psi_1 / sqrt(sum(abs(psi_1).^2)); % normalized state

    n2 = 1;
    Hn2 = hermiteH(n2, X);
    psiH_2 = (1 / sqrt(2^n2 * factorial(n2) * sqrt(pi))) * Hn2 .* exp(-X.^2 / 2);
    psiH_2 = psiH_2 / sqrt(sum(abs(psiH_2).^2));

    % Split-step propagator definitions
    UV = exp(-1i * (X.^2 / 2 + lambda * X.^4) * dt / 2);
    UT = exp(-1i * (P.^2 / 2) * dt);
    psi_0 = psi; % Initialize state for evolution

    % Initialize transition probability tracking for each time step
    current_max_prob = 0;

    for m = 1:M
        t = m * dt;
        V_t = A * sin(X) * (cos(omega1 * t) + cos(omega2 * t)) + lambda * X.^4;
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

    % Store the maximum transition probability for this lambda
    max_transition_probabilities(lambda_idx) = current_max_prob;
end

% Plotting the bifurcation diagram for lambda vs maximum transition probability
figure;
scatter(lambdas, max_transition_probabilities, '.');
xlabel('\lambda');
ylabel('Maximum Transition Probability');
title('Bifurcation Diagram for Transition Probability vs Nonlinearity (\lambda)');
grid on;
