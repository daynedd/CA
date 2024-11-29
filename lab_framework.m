
a = -20;                       % Left end point
b = +20;                       % Right end point
L = b-a;                        % Width of the space
N = 512;                       % No. of cells
X = a+L*(0:N-1)/N;                % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1]; % Dimensionless momentum
T=10*pi;                         % Time duration of the evolution
M =40^3;                     % Total No. of steps in the evolution
dt = T/M;                     % Time step
A =0.1;
omega =0.95;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UV = exp(-1i*(X.^2/2)*dt/2);    % One-step propagator in position space, only taking diagonal form
UV = exp(-1i*(X.^2/2)*dt/2);
UT = exp(-1i*(P.^2/2)*dt);       % One-setp propagator in momentum space
% note, hbar=1 in our dimensionless units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state
%   As a typical example, we consider the initial state to be a Gaussian
%   wavepacket located at X0
figure;
set(gcf, 'Position', [100, 100, 1200, 600]);
n=0
Hn = hermiteH(n,X);
psi_1 =(1/sqrt(2^n*factorial(n)*sqrt(pi)))*Hn.*exp(-X.^2/2);
psi=psi_1/sqrt(sum(abs(psi_1).^2));%normalized state
plot (X(1:N),abs(psi(1:N)).^2,'r');   % plotting initial state
hold on
n2 = 1;
Hn2 = hermiteH(n2, X);
psiH_2 = (1 / sqrt(2^n2 * factorial(n2) * sqrt(pi))) * Hn2 .* exp(-X.^2 / 2);
psiH_2 = psiH_2 / sqrt(sum(abs(psiH_2).^2));
plot(X, abs(psiH_2).^2, 'g');
natural_frequency=(n2-n)*1
%plot (P(1:N),abs(phi(1:N)).^2)
transition_probabilities = zeros(1, M);
transition_probabilities_Pert = zeros(1, M);
prefactor = abs((A/2).*sum(conj(psiH_2) .*sin(X).* psi))^2 
psi_0=psi;
for m = 1:M
    t = m*dt;
    V_t = A.*sin(X)*cos(omega.*t);
    UV_t=exp(-1i*(X.^2/2+V_t)*dt/2);
    psi_1 = UV_t .* psi_0;
    phi_2 = fft(psi_1);
    phi_3 = UT .* phi_2;
    psi_3 = ifft(phi_3);
    psi_4 = UV_t .* psi_3;
    psi_0 = psi_4;
    transition_probabilities(m) = abs(sum(conj(psiH_2) .* psi_0) )^2;
    transition_probabilities_Pert(m)=prefactor*abs(((exp(1i * (omega + natural_frequency) * t) - 1) / (omega + natural_frequency))+((exp(1i * (natural_frequency - omega) * t) - 1) / (natural_frequency - omega))).^2;
end
psi=psi_0; %final state updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot (X(1:N),abs(psi(1:N)).^2,'b')  %plotting the final state profile





xlabel('X');
ylabel('Probability Density');
legend('Initial State(n=0)','target state(n=1)', 'Final State');
figTitle1 = sprintf('QHO 1st eigenstate Evolution with Time-Periodic Perturbation, A=%.2f, omega=%.2f', A, omega);
title(figTitle1);
legend( 'Initial State(n=0)','Target State(n=1)', 'Final State');
hold off;


% Plotting the transition probability over time
figure;
set(gcf, 'Position', [100, 100, 1200, 600]);
plot((1:M) * dt, transition_probabilities);
hold on
plot((1:M) * dt, transition_probabilities_Pert);
xlabel('Time');
ylabel('Transition Probability to First Excited State');
legend( 'numercial','perturbation Theory');
figTitle2 = sprintf('Transition Probability from Ground State to First Excited State Over Time, A=%.2f, omega=%.2f', A, omega);
title(figTitle2);
hold off;



=======
clear; 
a = -20;                       % Left end point
b = +20;                       % Right end point
L = b-a;                        % Width of the space
N = 512;                       % No. of cells
X = a+L*(0:N-1)/N;                % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1]; % Dimensionless momentum
T=5*pi;                         % Time duration of the evolution
M =40^3;                     % Total No. of steps in the evolution
dt = T/M;                     % Time step
A =0.1;
omega =0.95;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define vectors to store split step propagators in position and
%   momentum space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UV = exp(-1i*(X.^2/2)*dt/2);    % One-step propagator in position space, only taking diagonal form
UV = exp(-1i*(X.^2/2)*dt/2);
UT = exp(-1i*(P.^2/2)*dt);       % One-setp propagator in momentum space
% note, hbar=1 in our dimensionless units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Define the initial state
%   As a typical example, we consider the initial state to be a Gaussian
%   wavepacket located at X0
figure;
set(gcf, 'Position', [100, 100, 1200, 600]);
n=0
Hn = hermiteH(n,X);
psi_1 =(1/sqrt(2^n*factorial(n)*sqrt(pi)))*Hn.*exp(-X.^2/2);
psi=psi_1/sqrt(sum(abs(psi_1).^2));%normalized state
plot (X(1:N),abs(psi(1:N)).^2,'r');   % plotting initial state
hold on
n2 = 1;
Hn2 = hermiteH(n2, X);
psiH_2 = (1 / sqrt(2^n2 * factorial(n2) * sqrt(pi))) * Hn2 .* exp(-X.^2 / 2);
psiH_2 = psiH_2 / sqrt(sum(abs(psiH_2).^2));
plot(X, abs(psiH_2).^2, 'g');
natural_frequency=(n2-n)*1
%plot (P(1:N),abs(phi(1:N)).^2)
transition_probabilities = zeros(1, M);
transition_probabilities_Pert = zeros(1, M);
prefactor = abs((A/2).*sum(conj(psiH_2) .*sin(X).* psi))^2 
psi_0=psi;
for m = 1:M
    t = m*dt;
    V_t = A.*sin(X)*cos(omega.*t);
    UV_t=exp(-1i*(X.^2/2+V_t)*dt/2);
    psi_1 = UV_t .* psi_0;
    phi_2 = fft(psi_1);
    phi_3 = UT .* phi_2;
    psi_3 = ifft(phi_3);
    psi_4 = UV_t .* psi_3;
    psi_0 = psi_4;
    transition_probabilities(m) = abs(sum(conj(psiH_2) .* psi_0) )^2;
    transition_probabilities_Pert(m)=prefactor*abs(((exp(1i * (omega + natural_frequency) * t) - 1) / (omega + natural_frequency))+((exp(1i * (natural_frequency - omega) * t) - 1) / (natural_frequency - omega))).^2;
end
psi=psi_0; %final state updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot (X(1:N),abs(psi(1:N)).^2,'b')  %plotting the final state profile





xlabel('X');
ylabel('Probability Density');
legend('Initial State(n=0)','target state(n=1)', 'Final State');
figTitle1 = sprintf('QHO 1st eigenstate Evolution with Time-Periodic Perturbation, A=%.2f, omega=%.2f', A, omega);
title(figTitle1);
legend( 'Initial State(n=0)','Target State(n=1)', 'Final State');
hold off;

% Save the first figure
save_dir = 'C:\Users\ASUS\Desktop\projectdata_4230';
filename1 = sprintf('%s\\[%0.2f]_[%0.2f]_state_profile.jpg', save_dir, A, omega);
saveas(gcf, filename1);

% Plotting the transition probability over time
figure;
set(gcf, 'Position', [100, 100, 1200, 600]);
plot((1:M) * dt, transition_probabilities);
hold on
plot((1:M) * dt, transition_probabilities_Pert);
xlabel('Time');
ylabel('Transition Probability to First Excited State');
legend( 'numercial','perturbation Theory');
figTitle2 = sprintf('Transition Probability from Ground State to First Excited State Over Time, A=%.2f, omega=%.2f', A, omega);
title(figTitle2);
hold off;

% Save the second figure
filename2 = sprintf('%s\\[%0.2f]_[%0.2f]_transition_probability.jpg', save_dir, A, omega);
%saveas(gcf, filename2);

>>>>>>> 8bdc1254ca6006729352568a8feddee52b3c0c2e
toc;