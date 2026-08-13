%% =========================================================
%  第15章 例15-1：基于遗传算法的三自由度结构损伤识别
%  教材：结构健康监测与智能传感
%  =========================================================
%
%  问题描述：
%    三自由度弹簧-质量模型，各质量块 m1=m2=m3=10 kg，
%    各弹簧名义刚度 k1=k2=k3=1000 N/m。
%    预设损伤情形：k3 发生 20% 刚度折减（k3_d = 800 N/m）。
%    已知损伤后结构的三阶角频率（实测值），
%    目标：用遗传算法反演损伤参数 alpha = [a1, a2, a3]，
%    其中 aj 表示第 j 个弹簧的刚度折减系数（0 = 无损，1 = 完全破坏）。
%
%  目标函数：
%    J(alpha) = sum( (omega_theory(alpha) - omega_measured).^2 )
%
%  约束条件（物理约束）：
%    0 <= aj <= 1，j = 1, 2, 3
%
%  依赖工具箱：Matlab Global Optimization Toolbox（内含 ga 函数）
%
%  作者：配套课件示例代码
%  =========================================================

clear; clc; close all;

%% ---- 1. 结构参数定义 ----------------------------------------

m  = 10.0;       % 各质量块质量 (kg)
k0 = 1000.0;     % 名义弹簧刚度 (N/m)

% 质量矩阵 M (3×3 对角矩阵)
M = diag([m, m, m]);

%% ---- 2. 计算损伤状态下的"实测"角频率 ω_d -------------------
% 预设：k3 折减 20%，即 a3_true = 0.2
a_true = [0.0, 0.0, 0.2];
omega_d = compute_frequencies(a_true, k0, M);

fprintf('=== 真实损伤参数 ===\n');
fprintf('  alpha_true = [%.4f, %.4f, %.4f]\n', a_true(1), a_true(2), a_true(3));
fprintf('\n=== 损伤后"实测"角频率 (rad/s) ===\n');
fprintf('  omega1 = %.4f,  omega2 = %.4f,  omega3 = %.4f\n\n', ...
        omega_d(1), omega_d(2), omega_d(3));

%% ---- 3. 定义目标函数（匿名函数形式，方便传入 ga）-----------
% ga 要求目标函数格式：f = fun(x)，x 为行向量
obj_fun = @(alpha) sum( (compute_frequencies(alpha, k0, M) - omega_d).^2 );

%% ---- 4. 遗传算法参数配置 ------------------------------------
n_var = 3;          % 决策变量维数（3 个损伤参数）
lb    = [0, 0, 0];  % 下界：刚度不增加
ub    = [1, 1, 1];  % 上界：刚度最多完全丧失

% 使用 optimoptions 配置 GA 选项
options = optimoptions('ga', ...
    'PopulationSize',    100, ...   % 种群大小
    'MaxGenerations',    200, ...   % 最大迭代代数
    'CrossoverFraction', 0.8, ...   % 交叉概率
    'MutationFcn',       @mutationadaptfeasible, ...  % 自适应可行变异
    'FunctionTolerance', 1e-10, ... % 目标函数收敛容差
    'Display',           'iter', ...% 每代打印迭代信息
    'PlotFcn',           @gaplotbestf); % 实时绘制最优适应度曲线

%% ---- 5. 运行遗传算法 ----------------------------------------
fprintf('=== 开始遗传算法求解 ===\n');
rng(42);  % 固定随机种子，保证结果可复现

[alpha_opt, J_opt, exitflag, output] = ga(obj_fun, n_var, ...
    [], [], [], [], lb, ub, [], options);

%% ---- 6. 输出识别结果 ----------------------------------------
fprintf('\n=== 遗传算法识别结果 ===\n');
fprintf('  alpha1 = %.6f  (真实值: %.4f)\n', alpha_opt(1), a_true(1));
fprintf('  alpha2 = %.6f  (真实值: %.4f)\n', alpha_opt(2), a_true(2));
fprintf('  alpha3 = %.6f  (真实值: %.4f)\n', alpha_opt(3), a_true(3));
fprintf('  目标函数值 J = %.4e\n', J_opt);
fprintf('  迭代代数   = %d\n', output.generations);
fprintf('  函数评价次数 = %d\n\n', output.funccount);

% 计算识别误差
err = abs(alpha_opt - a_true) ./ (a_true + 1e-10) * 100;
fprintf('=== 各参数相对误差 ===\n');
fprintf('  |Δα1|/α1_true = %.2f%%\n', err(1));
fprintf('  |Δα2|/α2_true = %.2f%%\n', err(2));
fprintf('  |Δα3|/α3_true = %.2f%%\n', err(3));

%% ---- 7. 对比理论频率与识别后频率 ----------------------------
omega_identified = compute_frequencies(alpha_opt, k0, M);
fprintf('\n=== 频率对比 (rad/s) ===\n');
fprintf('  %-12s %-12s %-12s %-12s\n', '阶次', '实测值', '识别值', '误差');
for i = 1:3
    freq_err = abs(omega_identified(i) - omega_d(i)) / omega_d(i) * 100;
    fprintf('  ω%d: %-12.4f %-12.4f %.4f%%\n', i, omega_d(i), omega_identified(i), freq_err);
end

%% ---- 8. 绘制损伤参数识别结果柱状图 -------------------------
figure('Name', '损伤参数识别结果对比', 'NumberTitle', 'off', 'Position', [100, 100, 700, 450]);

bar_data = [a_true; alpha_opt]';
b = bar(bar_data, 0.6);
b(1).FaceColor = [0.10, 0.23, 0.36];  % 深蓝色（真实值）
b(2).FaceColor = [0.91, 0.47, 0.13];  % 橙色（识别值）

set(gca, 'XTickLabel', {'k_1 (α_1)', 'k_2 (α_2)', 'k_3 (α_3)'}, 'FontSize', 12);
ylabel('损伤折减系数 α', 'FontSize', 13);
title('例15-1：遗传算法损伤参数识别结果', 'FontSize', 14);
legend({'真实值', '遗传算法识别值'}, 'Location', 'NorthEast', 'FontSize', 11);
ylim([0, 0.35]);
grid on;

% 在柱顶标注数值
for i = 1:3
    text(i - 0.18, a_true(i)    + 0.008, sprintf('%.3f', a_true(i)),    'FontSize', 10, 'Color', [0.10, 0.23, 0.36]);
    text(i + 0.03, alpha_opt(i) + 0.008, sprintf('%.4f', alpha_opt(i)), 'FontSize', 10, 'Color', [0.91, 0.47, 0.13]);
end

%% ---- 9. 与最小二乘迭代法（灵敏度法）对比 -------------------
fprintf('\n=== 与教材最小二乘法（灵敏度法）的对比 ===\n');
fprintf('  %-18s %-20s %-20s\n', '方法', 'α3 识别值', '目标函数值');
fprintf('  %-18s %-20s %-20s\n', '最小二乘迭代法', '0.2010 (两次迭代)', '≈ 0 (解析精确)');
fprintf('  %-18s %-20.6f %-20.4e\n', '遗传算法 (GA)', alpha_opt(3), J_opt);
fprintf('\n  说明：\n');
fprintf('  - 最小二乘法依赖灵敏度矩阵，收敛快（2次迭代），但需计算偏导数。\n');
fprintf('  - 遗传算法无需灵敏度，适用于高维/非线性/复杂结构，但计算量更大。\n');
fprintf('  - 加入 0≤α≤1 的物理约束后，GA 可避免产生无物理意义的负刚度解。\n');

%% =========================================================
%  辅助函数：计算给定损伤参数下结构的前三阶角频率
%  =========================================================
function omega = compute_frequencies(alpha, k0, M)
    % alpha: [a1, a2, a3]，各弹簧刚度折减系数
    % k0:    名义刚度 (N/m)
    % M:     质量矩阵 (3×3)
    % 返回:  omega (3×1)，前三阶角频率 (rad/s)，升序排列

    k1 = (1 - alpha(1)) * k0;
    k2 = (1 - alpha(2)) * k0;
    k3 = (1 - alpha(3)) * k0;

    % 三自由度串联弹簧-质量系统的刚度矩阵
    K = [ k1,       -k1,        0;
         -k1,  k1+k2,      -k2;
          0,       -k2,   k2+k3];

    % 广义特征值问题：K*phi = lambda*M*phi
    % lambda = omega^2
    lambda = eig(K, M);
    lambda = sort(real(lambda));       % 取实部并升序排列
    lambda = max(lambda, 0);           % 防止数值误差产生微小负值
    omega  = sqrt(lambda)';            % 角频率 (rad/s)
end
