%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Numerical simulation of the evolution of a wavepacket in a 1D harmonic
%   trap with an added time-dependent potential using fast Fourier transform (FFT) method
%   Potential: V(t) = A * sin(x) * cos(ω * t)
%   Unit of energy: hbar * omega, where hbar is the Planck constant and
%   omega is the frequency of the trap
%   Unit of length: l = sqrt(hbar / (m * omega)), where sqrt(...) is the square
%   root function and m is the mass of the particle
%--------------------------------------------------------------------------
a = -20;                       % Left end point 
b = +20;                       % Right end point 
L = b - a;                     % Width of the space
N = 512;                       % No. of cells
X = a + L * (0:N-1) / N;       % Dimensionless coordinates
P = (2 * pi / L) * [0:N/2-1, -N/2:-1]; % Dimensionless momentum
T = 10 * pi;                   % Time duration of the evolution
M = 10^3;                      % Total No. of steps in the evolution
dt = T / M;                    % Time step
A = 3;                         % Amplitude of the time-dependent potential
omega = 5.0;                   % Frequency of the time-dependent potential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split-step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UV = exp(-1i * (X.^2 / 2) * dt / 2);    % One-step propagator in position space (harmonic potential)
UT = exp(-1i * (P.^2 / 2) * dt);        % One-step propagator in momentum space
% note, hbar = 1 in our dimensionless units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state using Hermite polynomial (Analytical Eigenstate)
%   We consider the initial state to be the ground state of the harmonic oscillator
n = 0;  % Quantum number for ground state
psi_hermite = hermiteH(n, X) .* exp(-X.^2 / 2);  % Analytical form of eigenstate
psi = psi_hermite / sqrt(sum(abs(psi_hermite).^2));   % Normalized state

plot(X, abs(psi).^2, 'DisplayName', 'Initial State');   % Plotting initial state
legend('IC','Location','southwest')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full Evolution with Harmonic Oscillator and Time-Dependent Potential
psi_0 = psi;
P_transition_time = zeros(1, M);  % To store transition probabilities over time

% Define the first excited state using Hermite polynomial
n_1 = 1;  % Quantum number for first excited state
psi_hermite_1 = hermiteH(n_1, X) .* exp(-X.^2 / 2);  % Analytical form of first excited state
psi_hermite_1 = psi_hermite_1 / sqrt(sum(abs(psi_hermite_1).^2));   % Normalized first excited state

for m = 1:M
    t = m * dt;  % Current time
    % Update the potential propagator with the time-dependent potential
    V_t = A * sin(X) * cos(omega * t);
    UV_t = exp(-1i * (X.^2 / 2 + V_t) * dt / 2);  % Time-dependent potential propagator
    % Split-step evolution
    psi_1 = UV_t .* psi_0;           % Apply potential part (first half step)
    phi_2 = fft(psi_1);              % Transform to momentum space
    phi_3 = UT .* phi_2;             % Apply kinetic part
    psi_3 = ifft(phi_3);             % Transform back to position space
    psi_4 = UV_t .* psi_3;           % Apply potential part (second half step)    
    psi_0 = psi_4;  % Prepare for the next cycle    
    
    % Calculate the transition probability to the first excited state
    P_transition_time(m) = abs(sum(conj(psi_hermite_1) .* psi_0)).^2;
end

psi = psi_0;  % Final state updated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the final state profiles
plot(X, abs(psi).^2, 'DisplayName', 'Full Evolution');

xlabel('Position X')
ylabel('Probability Density')
title('Wavepacket Evolution in 1D Harmonic Trap with Time-Dependent Potential')
legend('show', 'Location', 'southwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the transition probability as a function of time
figure;
plot((1:M) * dt, P_transition_time, 'DisplayName', 'Transition Probability');
xlabel('Time')
ylabel('Transition Probability')
title('Transition Probability from Ground State to First Excited State')
legend('show')
