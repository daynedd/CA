%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical simulation of the evolution of a wavepacket in a 1D harmonic
% trap with an added torus-like time-dependent potential using fast Fourier transform (FFT) method
% Potential: V(t) = A * sin(x) * cos(θ1) + B * cos(x) * sin(θ2)
% Then, visualize the result on a torus with the full surface represented by colormap.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for space and simulation
a = -20;                       % Left end point 
b = +20;                       % Right end point 
L = b - a;                     % Width of the space
N = 512;                       % No. of cells
X = a + L * (0:N-1) / N;       % Dimensionless coordinates
P = (2 * pi / L) * [0:N/2-1, -N/2:-1]; % Dimensionless momentum
T = 5 * pi;                    % Time duration of the evolution
M = 10^3;                      % Total No. of steps in the evolution
dt = T / M;                    % Time step
A = 0.5;                       % Amplitude for the sin term
B = 0.3;                       % Amplitude for the cos term in torus-like potential
omega1 = 0.01;                 % Frequency for θ1
omega2 = 0.02;                 % Frequency for θ2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define vectors to store split-step propagators in position and momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UV = exp(-1i * (X.^2 / 2) * dt / 2);    % One-step propagator in position space (harmonic potential)
UT = exp(-1i * (P.^2 / 2) * dt);        % One-step propagator in momentum space
% note, hbar = 1 in our dimensionless units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the initial state
% As a typical example, we consider the initial state to be a Gaussian wavepacket located at X0
X0 = 4.0;                             % Initial position of the wave packet
sigma = 1.0;                          % Width of the initial wavepacket
psiprep = exp(-(X - X0).^2 / (2 * sigma^2));  % Gaussian state
psi = psiprep / sqrt(sum(abs(psiprep).^2));   % Normalized state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full Evolution with Harmonic Oscillator and Torus-like Time-Dependent Potential
psi_0 = psi;
for m = 1:M
    t = m * dt;  % Current time
    
    % Define the torus-like parameters theta1 and theta2
    theta1 = omega1 * t;  
    theta2 = omega2 * t;  
    
    % Update the time-dependent potential with both periodic components
    V_t = A * sin(X) .* cos(theta1) + B * cos(X) .* sin(theta2);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the torus using the full surface with colormap representing the maximum probability
theta1_vals = linspace(0, 2 * pi, 100);
theta2_vals = linspace(0, 2 * pi, 100);
[Theta1, Theta2] = meshgrid(theta1_vals, theta2_vals);

% Constructing the torus
R = 3; % Major radius
r = 1; % Minor radius

X_torus = (R + r * cos(Theta2)) .* cos(Theta1);
Y_torus = (R + r * cos(Theta2)) .* sin(Theta1);
Z_torus = r * sin(Theta2);

% Map the final probability density onto the torus
% Using the absolute square of the wavefunction to represent probability density
prob_density = abs(psii).^2;
prob_density_interp = interp1(X, prob_density, r * cos(Theta2), 'linear', 0); % Interpolate prob_density onto torus surface

% Find the maximum value in each inner circle (fixed theta1) for color emphasis
max_prob_indices = zeros(size(theta1_vals));
max_prob_values = zeros(size(theta1_vals));

for i = 1:length(theta1_vals)
    [max_val, max_idx] = max(prob_density_interp(i, :));
    max_prob_indices(i) = max_idx;
    max_prob_values(i) = max_val;
end

% Construct a colormap that emphasizes the maximum values
% Amplify the color of the maximum values to make them stand out
highlighted_prob_density = prob_density_interp;
for i = 1:length(theta1_vals)
    highlighted_prob_density(i, max_prob_indices(i)) = max_prob_values(i) * 1.5; % Emphasize maximum value
end

% Plot the torus with the highlighted colormap
figure;
surf(X_torus, Y_torus, Z_torus, highlighted_prob_density);
shading interp
colormap jet
colorbar
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Torus Surface with Maximum Probability Emphasized by Colormap')
axis equal
hold off
