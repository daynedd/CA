%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 参数化设置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = 0.01;                   % 时间依赖频率
A = 4.5;                        % 位置依赖的振幅
a = -20;                        % 左端点
b = 20;                         % 右端点
L = b - a;                      % 宽度
N = 512;                        % 网格点数
X = a + L * (0:N-1) / N;        % 空间坐标
P = (2 * pi / L) * [0:N/2-1, -N/2:-1]; % 动量坐标
T = 5 * pi;                     % 总时间
M = 1000;                       % 时间步数
dt = T / M;                     % 时间步长

% 初始状态：高斯波包
X0 = 4.0;
sigma = 1.0;
psiprep = exp(-(X - X0).^2 / (2 * sigma^2)); % 高斯波包
psi = psiprep / sqrt(sum(abs(psiprep).^2));  % 归一化

% 定义位置和动量空间传播子
UV = exp(-1i * (X.^2 / 2) * dt / 2);    % 位置空间中的半步传播子
UT = exp(-1i * (P.^2 / 2) * dt);        % 动量空间传播子

% 创建二维网格参数化角度
num_points = 100;
t_vals = linspace(0, 2 * pi / omega, num_points);
x_vals = linspace(0, 2 * pi, num_points);
[T_vals, X_vals] = meshgrid(t_vals, x_vals);

% 准备存储不同时间和位置下的概率密度
prob_density_full = zeros(N, num_points);  % 使用 N 而不是 num_points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 波函数的时间演化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:M
    t = m * dt;
    % 更新时间依赖势能的传播子
    V_t = A * sin(X) * cos(omega * t);
    UV_t = exp(-1i * (X.^2 / 2 + V_t) * dt / 2);

    % 裂步法进行传播
    psi_1 = UV_t .* psi;
    phi_2 = fft(psi_1);
    phi_3 = UT .* phi_2;
    psi_3 = ifft(phi_3);
    psi_4 = UV_t .* psi_3;

    psi = psi_4;  % 更新波函数为下一次演化的初始状态

    % 计算并存储概率密度（只在特定时间步进行采样）
    if mod(m, M / num_points) == 0
        idx = m / (M / num_points);
        prob_density_full(:, idx) = abs(psi).^2; % 存储整个空间的概率密度
    end
end

% 将概率密度重采样为 100x100 的大小以匹配绘图网格
prob_density = imresize(prob_density_full, [num_points, num_points]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 环面的参数化和绘制
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 10; % 环的主半径
r = 3;  % 管的半径

% 计算环面坐标
theta = T_vals; % 时间作为内圆角度
phi = X_vals;   % 位置作为外圆角度

X_torus = (R + r * cos(phi)) .* cos(theta);
Y_torus = (R + r * cos(phi)) .* sin(theta);
Z_torus = r * sin(phi);

% 绘制环面并用颜色表示概率密度
figure;
surf(X_torus, Y_torus, Z_torus, prob_density, 'FaceColor', 'interp', 'EdgeColor', 'none');
colormap(viridis);  % 使用 viridis 颜色图
colorbar;           % 显示颜色条

% 设置图形参数
title('3D Torus (S1 x S1) with Probability Density');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
axis equal;
