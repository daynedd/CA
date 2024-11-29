clear;
a = -20;                       % Left end point
b = +20;                       % Right end point
L = b-a;                        % Width of the space
N = 512;                       % No. of cells
X = a+L*(0:N-1)/N;                % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1]; % Dimensionless momentum
T = 20*pi;                      % Time duration of the evolution
M = 40^3;                       % Total No. of steps in the evolution
dt = T/M;                       % Time step
A = 0.5;
omega = 0.5;
epsilon = 1e-10;                % Small value to avoid division by zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UV = exp(-1i*(X.^2/2)*dt/2);    % One-step propagator in position space
UT = exp(-1i*(P.^2/2)*dt);      % One-step propagator in momentum space
% note, hbar=1 in our dimensionless units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state
%   As a typical example, we consider the initial state to be a Gaussian
%   wavepacket located at X0
figure;
set(gcf, 'Position', [100, 100, 1200, 600]);
n = 0;
Hn = hermiteH(n, X);
psi_1 = (1/sqrt(2^n * factorial(n) * sqrt(pi))) * Hn .* exp(-X.^2 / 2);
psi = psi_1 / sqrt(sum(abs(psi_1).^2)); % normalized state
plot(X(1:N), abs(psi(1:N)).^2, 'r');   % plotting initial state
hold on

% Define target state (n=1)
n2 = 1;
Hn2 = hermiteH(n2, X);
psiH_2 = (1 / sqrt(2^n2 * factorial(n2) * sqrt(pi))) * Hn2 .* exp(-X.^2 / 2);
psiH_2 = psiH_2 / sqrt(sum(abs(psiH_2).^2));
plot(X, abs(psiH_2).^2, 'g');

natural_frequency = (n2 - n) * 1;

% Define transition probability storage vectors
transition_probabilities = zeros(1, M);
transition_probabilities_Pert = zeros(1, M);
transition_probabilities_SecondPert = zeros(1, M); % Adding second-order perturbation storage
transition_probabilities_Total = zeros(1, M); % Total transition probability

% Calculate prefactor for first-order perturbation theory
prefactor = abs((A/2) * sum(conj(psiH_2) .* sin(X) .* psi))^2;

% Calculate prefactor for second-order perturbation theory
second_order_prefactor = 0;
max_intermediate_state = 5; % Maximum intermediate state to consider for calculation
for n1 = 1:max_intermediate_state
    if n1 ~= n && n1 ~= n2 % Ensure intermediate state is not the same as initial or final state
        Hn1 = hermiteH(n1, X);
        psiH_1 = (1 / sqrt(2^n1 * factorial(n1) * sqrt(pi))) * Hn1 .* exp(-X.^2 / 2);
        psiH_1 = psiH_1 / sqrt(sum(abs(psiH_1).^2));

        % Contribution from intermediate state n1
        matrix_element_1 = sum(conj(psiH_1) .* sin(X) .* psi);
        matrix_element_2 = sum(conj(psiH_2) .* sin(X) .* psiH_1);
        energy_diff = natural_frequency - (n1 - n);
        contribution = abs(A * matrix_element_1 * matrix_element_2 / (energy_diff + epsilon))^2;
        second_order_prefactor = second_order_prefactor + contribution;
    end
end

fprintf('Second-order prefactor: %.5f\n', second_order_prefactor);

% Initial wave function
psi_0 = psi;

for m = 1:M
    t = m * dt;
    V_t = A * sin(X) * cos(omega * t);
    UV_t = exp(-1i * (X.^2 / 2 + V_t) * dt / 2);
    
    % Numerical evolution
    psi_1 = UV_t .* psi_0;
    phi_2 = fft(psi_1);
    phi_3 = UT .* phi_2;
    psi_3 = ifft(phi_3);
    psi_4 = UV_t .* psi_3;
    psi_0 = psi_4;

    % Store transition probabilities
    transition_probabilities(m) = abs(sum(conj(psiH_2) .* psi_0))^2;

    % First-order perturbation theory probability
    transition_probabilities_Pert(m) = prefactor * abs(((exp(1i * (omega + natural_frequency) * t) - 1) / (omega + natural_frequency + epsilon)) + ((exp(1i * (natural_frequency - omega) * t) - 1) / (natural_frequency - omega + epsilon))).^2;
    
    % Second-order perturbation theory probability using RWA (no intermediate states considered explicitly)
    % Using the rotating wave approximation (RWA) to focus on resonant terms
    transition_probabilities_SecondPert(m) = second_order_prefactor *abs((exp(1i * (2 * omega) * t) - 1) / (2 * omega + epsilon)).^2;
    
    % Total transition probability (sum of first and second order, capped at 1)
    transition_probabilities_Total(m) = abs(min(transition_probabilities_Pert(m) + transition_probabilities_SecondPert(m), 1));
end

psi = psi_0; % Final state update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting initial, target, and final states
plot(X(1:N), abs(psi(1:N)).^2, 'b');  % Plotting the final state profile
xlabel('X');
ylabel('Probability Density');
legend('Initial State(n=0)', 'Target State(n=1)', 'Final State');
figTitle1 = sprintf('QHO 1st eigenstate Evolution with Time-Periodic Perturbation, A=%.2f, omega=%.2f', A, omega);
title(figTitle1);
legend('Initial State(n=0)', 'Target State(n=1)', 'Final State');
hold off;

% Plotting the transition probability over time
figure;
set(gcf, 'Position', [100, 100, 1200, 600]);
plot((1:M) * dt, transition_probabilities, 'r'); % Numerical result
hold on;
plot((1:M) * dt, transition_probabilities_Pert, 'g'); % First-order perturbation
plot((1:M) * dt, transition_probabilities_SecondPert, 'b'); % Second-order perturbation (RWA)
xlabel('Time');
ylabel('Transition Probability to First Excited State');
legend('Numerical', 'First-order Perturbation Theory', 'Second-order Perturbation Theory (RWA)');
figTitle2 = sprintf('Transition Probability from Ground State to First Excited State Over Time, A=%.2f, omega=%.2f', A, omega);
title(figTitle2);
hold off;

toc;