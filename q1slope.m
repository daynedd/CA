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

% Set of A values from 0.1 to 10 with 100 elements
A_values = linspace(0.1, 10, 100);
% Set of omega values from 0.01 to 1 with 10 elements
omega_values = linspace(0.01, 1, 10);

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

% Initialize matrix to store maximum probability density
max_probability_density = zeros(length(A_values), length(omega_values));

for idx_a = 1:length(A_values)
    A = A_values(idx_a);
    
    for idx_omega = 1:length(omega_values)
        omega = omega_values(idx_omega);
        
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
        max_probability_density(idx_a, idx_omega) = max(abs(psii).^2);
    end
end

% Toroidal Projection with Inner Circle (A * X) and Outer Circle (omega * t)
figure;
[A_mesh, omega_mesh] = meshgrid(A_values, omega_values);
[X_mesh, t_mesh] = meshgrid(X, linspace(0, T, length(X)));

% Adjust dimensions for compatibility
A_mesh_extended = repmat(A_values(:), 1, length(X));
X_mesh_extended = repmat(X, length(A_values), 1);
omega_mesh_extended = repmat(omega_values(:), 1, length(X));
t_mesh_extended = repmat(linspace(0, T, length(X)), length(omega_values), 1);

% Calculate toroidal coordinates
R = 15;  % Major radius of the torus
r = max(max_probability_density(:)) + 2;  % Minor radius based on max probability density

% Coordinates on torus with A * X as inner circle and omega * t as outer circle
inner_radius = A_mesh_extended .* X_mesh_extended;  % Inner circle defined by A * X
outer_angle = omega_mesh_extended .* t_mesh_extended;  % Outer circle defined by omega * t

% Calculate torus coordinates in 3D
x_torus = (R + r * cos(inner_radius)) .* cos(outer_angle);
y_torus = (R + r * cos(inner_radius)) .* sin(outer_angle);
z_torus = r * sin(inner_radius);

% Reshape max_probability_density to match torus coordinate sizes for plotting
color_data = repmat(max(max_probability_density, [], 2), 1, length(X));

% Plotting the toroidal projection
surf(x_torus, y_torus, z_torus, color_data);
xlabel('X (Torus Projection)');
ylabel('Y (Torus Projection)');
zlabel('Z (Torus Projection)');
colormap(jet);
colorbar;
title('Toroidal Representation with Inner Circle (A * X) and Outer Circle (omega * t)');
view(3);
shading interp;
