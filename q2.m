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
T = 5 * pi;                    % Time duration of the evolution
M = 10^3;                      % Total No. of steps in the evolution
dt = T / M;                    % Time step
A = 0.1;                       % Amplitude of the time-dependent potential
omega = 1.0;                   % Frequency of the time-dependent potential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split-step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UV = exp(-1i * (X.^2 / 2) * dt / 2);    % One-step propagator in position space (harmonic potential)
UT = exp(-1i * (P.^2 / 2) * dt);        % One-step propagator in momentum space
% note, hbar = 1 in our dimensionless units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state as an excited state of the harmonic oscillator
n_excited = 1;  % Set the quantum number for the desired excited state (e.g., n = 1 for the first excited state)
psi_ho_excited = hermiteH(n_excited, X) .* exp(-X.^2 / 2);  % Harmonic oscillator eigenstate (Hermite polynomial)
psi_ho_excited = psi_ho_excited / sqrt(sum(abs(psi_ho_excited).^2));  % Normalized initial state

plot(X, abs(psi_ho_excited).^2);   % Plotting initial state
hold on

psi_0 = psi_ho_excited;
for m = 1:M
    t = m * dt;  % Current time
    % Update the potential propagator with the time-dependent potential
    V_t = A * sin(X) * cos(omega * t);
    UV_t = exp(-1i * V_t * dt / 2);  % Time-dependent potential propagator
    
    % Split-step evolution
    psi_1 = UV_t .* psi_0;           % Apply potential part (first half step)
    phi_2 = fft(psi_1);              % Transform to momentum space
    phi_3 = UT .* phi_2;             % Apply kinetic part
    psi_3 = ifft(phi_3);             % Transform back to position space
    psi_4 = UV_t .* psi_3;           % Apply potential part (second half step)
    
    psi_0 = psi_4;  % Prepare for the next cycle    
end

psi = psi_0;  % Final state updated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Projection onto eigenstates of the harmonic oscillator
n_max = 5;  % Maximum number of eigenstates to consider
coefficients = zeros(1, n_max);

for n = 0:n_max-1
    % Calculate the nth eigenstate of the harmonic oscillator
    psi_n = hermiteH(n, X) .* exp(-X.^2 / 2);
    psi_n = psi_n / sqrt(sum(abs(psi_n).^2));  % Normalize the eigenstate
    
    % Project the final wavefunction onto the nth eigenstate
    coefficients(n+1) = sum(conj(psi_n) .* psi) * (X(2) - X(1));
end

% Plotting the projection coefficients
figure;
stem(0:n_max-1, abs(coefficients).^2);
xlabel('Eigenstate Index n');
ylabel('Projection Probability |<ψ_n|ψ>|^2');
title('Projection of Final State onto Harmonic Oscillator Eigenstates');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
