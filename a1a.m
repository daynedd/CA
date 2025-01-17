%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Numerical simulation of the evolution of a wavepacket in a 1D harmonic
%   trap with an added time-dependent potential using fast Fourier transform (FFT) method
%   Potential: V(t) = A * sin(x) * cos(omega * t)
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
omega = 0.01;                  % Frequency of the time-dependent potential

% Set of A values from 0.1 to 10 with 5 elements
A_values = linspace(0.1, 10, 500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split-step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UV = exp(-1i * (X.^2 / 2) * dt / 2);    % One-step propagator in position space (harmonic potential)
UT = exp(-1i * (P.^2 / 2) * dt);        % One-step propagator in momentum space
% note, hbar = 1 in our dimensionless units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state
%   As a typical example, we consider the initial state to be a Gaussian
%   wavepacket located at X0
X0 = 4.0;
sigma = 1.0;  % sigma is the width of the initial wavepacket
psiprep = exp(-(X - X0).^2 / (2 * sigma^2));  % Gaussian state
psi = psiprep / sqrt(sum(abs(psiprep).^2));   % Normalized state

% Define meshgrid for 3D plot
[X_grid, A_grid] = meshgrid(X, A_values);
probability_density = zeros(length(A_values), length(X));

for idx = 1:length(A_values)
    A = A_values(idx);
    
    % Full Evolution with Harmonic Oscillator and Time-Dependent Potential
    psi_0 = psi;
    for m = 1:M
        t = m * dt;  % Current time
        % Update the potential propagator with the time-dependent potential
        V_t = A * sin(X) .* cos(omega * t);
        UV_t = exp(-1i * (X.^2 / 2 + V_t) * dt / 2);  % Time-dependent potential propagator
        
        % Split-step evolution
        psi_1 = UV_t .* psi_0;           % Apply potential part (first half step)
        phi_2 = fft(psi_1);              % Transform to momentum space
        phi_3 = UT .* phi_2;             % Apply kinetic part
        psi_3 = ifft(phi_3);             % Transform back to position space
        psi_4 = UV_t .* psi_3;           % Apply potential part (second half step)
        
        psi_0 = psi_4;  % Prepare for the next cycle    
    end
    
    psii = psi_0;  % Final state updated
    probability_density(idx, :) = abs(psii).^2;
end

% Plotting the 3D plot where A is an independent axis
figure;
surf(X_grid, A_grid, probability_density);
xlabel('Position X')
ylabel('Amplitude A')
zlabel('Probability Density')
title('3D Plot of Wavepacket Evolution in 1D Harmonic Trap with Time-Dependent Potential')
shading interp
colorbar
