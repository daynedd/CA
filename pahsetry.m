% 经典谐振子相空间轨迹
clear;

% 参数设置
k = 1;                % 弹性系数（设为1）
m = 1;                % 质量（设为1）
omega = sqrt(k/m);    % 谐振角频率
T = 2 * pi / omega;   % 谐振周期
dt = 0.01;            % 时间步长
t_max = 10 * T;       % 总时间长度
N = floor(t_max / dt); % 时间步数，使用 floor 确保为整数

% 初始化位置和动量
x = 2;                % 初始位置（设为2）
v = 0;                % 初始速度（设为0）
p = m * v;            % 初始动量

% 用于存储位置和动量的数组
x_vals = zeros(1, N);
p_vals = zeros(1, N);

% 时间演化：使用欧拉法求解经典谐振子的运动方程
for i = 1:N
    % 保存当前位置和动量
    x_vals(i) = x;
    p_vals(i) = p;
    
    % 更新位置和动量（欧拉法）
    x = x + dt * (p / m);         % 位置更新
    p = p - dt * (k * x);         % 动量更新
end

% 绘制相空间轨迹
figure;
plot(x_vals, p_vals, '-b');
xlabel('Position X');
ylabel('Momentum P');
title('Phase Space Trajectory of Classical Harmonic Oscillator');
grid on;
axis equal;


