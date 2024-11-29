clear;
a = -20;                       % Left end point
b = +20;                       % Right end point
L = b - a;                     % Width of the space
N = 256;                       % No. of cells (reduced for better visualization performance)
X = a + L * (0:N-1) / N;       % Dimensionless coordinates
P = (2 * pi / L) * [0:N/2-1, -N/2:-1]; % Dimensionless momentum
dx = L / N;                    % Grid spacing in position space

% Range of lambda values for Wigner analysis
lambdas = [0.1, 0.5, 1.0];  % Different lambda values to explore nonlinearity impact

% Define phase space grid for Wigner distribution
x_range = linspace(a, b, N);
p_range = linspace(-3, 3, N);

figure;

for lambda_idx = 1:length(lambdas)
    lambda = lambdas(lambda_idx);
    
    % Initialize wavefunction as a Gaussian for given lambda
    n = 0;
    Hn = hermiteH(n, X);
    psi = (1 / sqrt(2^n * factorial(n) * sqrt(pi))) * Hn .* exp(-X.^2 / 2);
    psi = psi / sqrt(sum(abs(psi).^2)); % normalized state

    % Time evolution with nonlinear potential
    T = 5 * pi;                   % Time duration for the evolution
    M = 100;                      % Total number of steps in the evolution
    dt = T / M;                   % Time step
    A = 0.1;                      % Constant A value
    omega1 = 1.0526;              % Fixed value from previous analysis (resonance)
    omega2 = 0.2084;              % Fixed value from previous analysis (resonance)

    % Split-step propagator definitions
    UV = exp(-1i * (X.^2 / 2 + lambda * X.^4) * dt / 2);
    UT = exp(-1i * (P.^2 / 2) * dt);
    psi_0 = psi; % Initialize state for evolution

    % Time evolution
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
    end

    % Final wavefunction after evolution
    psi = psi_0;

    % Calculate the Wigner distribution for the final wavefunction
    W = zeros(length(x_range), length(p_range)); % Initialize Wigner function array

    for i = 1:length(x_range)
        x0 = x_range(i);
        for j = 1:length(p_range)
            p0 = p_range(j);
            % Calculate the Wigner function at (x0, p0)
            integral_sum = 0;
            for k = 1:length(X)
                x = X(k);
                y = x - x0;
                psi_star = conj(psi(k));
                psi_shift = interp1(X, psi, x + y / 2, 'linear', 0) .* conj(interp1(X, psi, x - y / 2, 'linear', 0));
                integral_sum = integral_sum + psi_shift * exp(1i * p0 * y);
            end
            W(i, j) = real(integral_sum * dx / pi);
        end
    end

    % Plotting the Wigner distribution for the current lambda
    subplot(1, length(lambdas), lambda_idx);
    imagesc(x_range, p_range, W');
    set(gca, 'YDir', 'normal');
    xlabel('Position X');
    ylabel('Momentum P');
    title(['\lambda = ', num2str(lambda)]);
    colormap jet;
    colorbar;
end

% Set a title for the entire figure
sgtitle('Wigner Distributions for Different \lambda Values');
