clear; 
a = -20;                       % Left end point
b = +20;                       % Right end point
L = b-a;                        % Width of the space
N = 512;                       % No. of cells
X = a+L*(0:N-1)/N;                % Dimensionless coordinates
P = (2*pi/L)*[0:N/2-1,-N/2:-1]; % Dimensionless momentum
T=10*pi;                         % Time duration of the evolution
M =10^3;                     % Total No. of steps in the evolution
dt = T/M;                     % Time step
A =3;
omega =5;
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
n=0
Hn = hermiteH(n,X);
psi_1 =(1/sqrt(2^n*factorial(n)*sqrt(pi)))*Hn.*exp(-X.^2/2);
psi=psi_1/sqrt(sum(abs(psi_1).^2));%normalized state
plot (X(1:N),abs(psi(1:N)).^2,'r');   % plotting initial state
hold on
%plot (P(1:N),abs(phi(1:N)).^2)
psi_0=psi;
for m = 1:M
    t = (m-1)*dt;
    V_t = A.*sin(X)*cos(omega.*t);
    UV_t=exp(-1i*(X.^2/2+V_t)*dt/2);
    psi_1 = UV_t .* psi_0;
    phi_2 = fft(psi_1);
    phi_3 = UT .* phi_2;
    psi_3 = ifft(phi_3);
    psi_4 = UV_t .* psi_3;
    psi_0 = psi_4;
end
psi=psi_0; %final state updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot (X(1:N),abs(psi(1:N)).^2,'b')  %plotting the final state profile
n2 = 1;
Hn2 = hermiteH(n2, X);
psi_2 = (1 / sqrt(2^n2 * factorial(n2) * sqrt(pi))) * Hn2 .* exp(-X.^2 / 2);
psi_2 = psi_2 / sqrt(sum(abs(psi_2).^2));
plot(X, abs(psi_2).^2, 'g');
transition_prob = abs(sum(conj(psi_2) .* psi)).^2;
natural_frequency=(n2-n)*1

disp(['Transition Probability from n=0 to n=1: ', num2str(transition_prob)]);
xlabel('X');
ylabel('Probability Density');
figTitle = 'QHO 1st eigenstate Evolution with Time-Periodic Perturbation,A=1,omega=1';
title(figTitle);
legend('Initial State(n=0)', 'Final State','target state(n=1)');
annotation('textbox', [0.0,0.0, 0.3, 0.1], 'String', ...
    {['Transition Probability (n=0 -> n=1): ', num2str(transition_prob)], ...
     ['Transition Natural Frequency: ', num2str(natural_frequency)]}, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'none');
filename = strrep(figTitle, ' ', '_');
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
uniqueFilename = [filename, '_', timestamp, '.pdf'];
savePath = fullfile(getenv('USERPROFILE'), 'Downloads', uniqueFilename);
%saveas(gcf, savePath);
hold off