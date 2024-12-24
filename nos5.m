% ================================ %
% 連立一次方程式に対する共役勾配法 %
% ================================ %

% SuiteSparse Matrix Collection: (https://sparse.tamu.edu/)
load('nos5.mat');
% load('nos7.mat');
% load('LFAT5.mat');   % condition number = 2.374510e+09

% 係数行列の情報
%（SuiteSparse Matrix Collectionからダウンロードしたもの）
% A = Problem.A;
A = sparse(Problem.A);
[n, ~] = size(A);

% 真の解の設定
%（成分をすべて1のベクトルとする）
x_true = ones(n, 1);

% 右辺項の設定
%（係数行列 * 真の解 とする）
b = A * x_true; % 右辺項

% 初期値の設定
%（成分をすべて0のベクトルとする）
x = zeros(n, 1);    % 初期値 x_0

% 許容誤差と最大反復回数の設定
eps = 10^(-12); % 許容誤差
% max_iter = 2 * n^2; % for nos7
max_iter = 2 * n; % for others

% 相対残差2ノルム，相対誤差2ノルム，相対誤差Aノルム、誤差の履歴を保存先
hist_relres_2 = zeros(max_iter + 1, 1);
hist_relerr_2 = zeros(max_iter + 1, 1);
hist_relerr_A = zeros(max_iter + 1, 1);
hist_solutions = zeros(n, max_iter + 1);   % 各反復での解 x (x_k) を保存
hist_err = zeros(n, max_iter + 1);

%%%%%%%%%%%%%%%%%%%% 問題の設定を出力 %%%%%%%%%%%%%%%%%%%%
fprintf('========================== \n');
fprintf('Coefficient Matrix A.\n');
fprintf('========================== \n');
fprintf('Name: %s\n', Problem.name);
fprintf('Kind: %s\n', Problem.kind);
fprintf('Size: %d×%d\n', n, n);
fprintf('Nnz : %d\n', nnz(A));
fprintf('\n');

fprintf('========================== \n');
fprintf('Exact Solution x_true. \n');
fprintf('========================== \n');
fprintf('x_true := [1,1,...,1]^T\n');
fprintf('\n');

fprintf('========================== \n');
fprintf('Right-Hand Side Vector b. \n');
fprintf('========================== \n');
fprintf('b := A * x_true\n');
fprintf('\n');

fprintf('========================== \n');
fprintf('Initial Guess x_0. \n');
fprintf('========================== \n');
fprintf('x_0 := [0,0,...,0]^T\n');
fprintf('\n');

% fprintf('========================== \n');
% fprintf('Spectrum Info. \n');
% fprintf('========================== \n');
% fprintf('lambda_min = %.4e\n', lambda_min);
% fprintf('lambda_max = %.4e\n', lambda_max);
% fprintf('\n');

% 許容誤差と最大反復回数の設定
fprintf('========================== \n');
fprintf('Computational Conditions. \n');
fprintf('========================== \n');
fprintf('Tolerance: Eps. := %.1e\n', eps);
fprintf('Maximum # Iter. := %d\n', max_iter);
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CG法の実行 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('************************** \n');
fprintf('Performing the CG Method. \n');
fprintf('************************** \n');

% 計測開始（収束履歴の記録も含まれることに注意）
tic;

% 右辺項bと真の解の2ノルム，および真の解のAノルムを保存
norm2_b = norm(b);
norm2_x_true = norm(x_true);
normA_x_true = sqrt(x_true' * A * x_true);

% 初期残差の設定
r = b; % 残差r = b - Ax (x_0=0)

% 初期の推定解、誤差、相対残差2ノルム，相対誤差2ノルム，相対誤差Aノルムの保存
hist_solutions(:, 1) = x;
hist_relres_2(1) = norm(r) / norm2_b;
err = x_true - x;
hist_err(:, 1) = err;
hist_relerr_2(1) = norm(err) / norm2_x_true;
hist_relerr_A(1) = sqrt(err' * A * err) / normA_x_true;

% 初期探索方向の設定
p = r;  % 初期残差方向

% rho_oldの計算（alphaの分子およびbetaの分母）
rho_old = dot(r, r); 

% 反復開始
for iter = 1:max_iter    

    % 使い回す共通項を計算
    w = A * p;
    sigma = dot(p, w);

    alpha = rho_old / sigma;    % アルファ係数
    x = x + alpha * p;  % 解の更新
    r = r - alpha * w;  % 残差を更新

    % 推定解、誤差、相対残差2ノルム，相対誤差2ノルム，相対誤差Aノルムの保存
    hist_solutions(:, iter+1) = x;
    hist_relres_2(iter+1) = norm(r) / norm2_b;
    err = x_true - x;
    hist_err(:, iter+1) = err;
    hist_relerr_2(iter+1) = norm(err) / norm2_x_true;
    hist_relerr_A(iter+1) = sqrt(err' * A * err) / normA_x_true;

    % 収束しているか判定
    if hist_relres_2(iter+1) < eps
        fprintf('... Convergence!\n');
        fprintf('\n');
        last_iter = iter;
        break;
    end
    
    rho_new = dot(r, r);
    beta = rho_new / rho_old;   % ベータ係数を計算
    rho_old = rho_new;
    p = r + beta*p; % 探索方向を更新
    
end % 反復終了

% 計測終了
time = toc;

% 収束しない場合（収束判定条件を満たさなかった場合）
if hist_relres_2(iter+1) >= eps
    fprintf('... No convergence!\n');
    fprintf('\n');
    last_iter = iter;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 固有値計算 (最小・最大固有値を取得), めっちゃ大規模な場合はeigsも視野 -> eigen_valeus.m
% lambda = eig(A);
% lambda_min = min(lambda);
% lambda_max = max(lambda);

% % xのAノルムおよび2ノルムの履歴を保存する配列
% norm_z_A = zeros(max_iter + 1, 1);
% norm_Az_2 = zeros(max_iter + 1, 1);
% % 解の存在する範囲を可視化する用の上限、下限を保存する配列
% upper_bound = zeros(max_iter + 1, 1);
% lower_bound = zeros(max_iter + 1, 1);

% for iter=1:last_iter+1
%     z = hist_err(:, iter);
%     norm_z_A(iter) = sqrt(z' * A * z);
%     norm_Az_2(iter) = norm(A * z, 2);
%     upper_bound(iter) = sqrt(lambda_max) * norm_z_A(iter);
%     lower_bound(iter) = sqrt(lambda_min) * norm_z_A(iter);
% end

% 固有値計算 (最小・最大固有値を取得), めっちゃ大規模な場合はeigsも視野 -> eigen_valeus.m
lambda = eig(A);
lambda_min = min(lambda);
lambda_max = max(lambda);

% xのAノルムおよび2ノルムの履歴を保存する配列
norm_err_A = zeros(max_iter + 1, 1);
% 解の存在する範囲を可視化する用の上限、下限を保存する配列
upper_bound = zeros(max_iter + 1, 1);
rel_err_A = zeros(max_iter + 1, 1);

norm_x_true_A = sqrt(x_true' * A * x_true);

kappa_A_2 = lambda_max / lambda_min;

for iter=1:last_iter+1
    err = hist_err(:, iter);
    norm_err_A(iter) = sqrt(err' * A * err); 
    rel_err_A(iter) = norm_err_A (iter) / norm_x_true_A;
    upper_bound(iter) = sqrt(kappa_A_2)  * rel_err_A(iter);
    % fprintf('%e', upper_bound(iter));
end

% 計算結果の表示
fprintf('========================== \n');
fprintf('Numerical Results. \n');
fprintf('========================== \n');
fprintf('# Iter.: %d\n', last_iter);
fprintf('Time[s]: %.3f\n', time);
fprintf('Relres_2norm = %.2e\n', hist_relres_2(last_iter+1));
fprintf('Relerr_2norm = %.2e\n', hist_relerr_2(last_iter+1));
fprintf('Relerr_Anorm = %.2e\n', hist_relerr_A(last_iter+1));
fprintf('========================== \n');
fprintf('\n');

% 収束履歴の表示
x_axis = 0:last_iter;

figure;
hold on, grid on;
plot(x_axis,hist_relres_2(1:last_iter+1),'-*','DisplayName',strcat('||r_k||_2/||b||_2'));
plot(x_axis,hist_relerr_2(1:last_iter+1),'-*','DisplayName',strcat('||e_k||_2/||x_{true}||_2'));
plot(x_axis,hist_relerr_A(1:last_iter+1),'-*','DisplayName',strcat('||e_k||_A/||x_{true}||_A'));
legend, box on;
title(Problem.name);
xlabel('Number of Iterations');
ylabel('Log_{10} of relative norm');
ylim(gca,[1e-13 1e+1]);
set(gca,...
    'FontSize',16,...
    'YScale','log',...
    'YTick',[1e-12 1e-09 1e-06 1e-3 1],...
    'YTickLabel',{'-12','-9','-6','-3','0'});
hold off;

% 実際は減らしたい誤差のA-normの状況はわからないが、このグラフから残差の2ノルムが減少していることで誤差のA-normも減らすことができている事がわかる。
% x_axis = 0:last_iter;

% figure;
% hold on, grid on;
% plot(x_axis, norm_Az_2(1:last_iter+1),'-*','DisplayName',strcat('||r_k||_2'));
% plot(x_axis, upper_bound(1:last_iter+1),'-*','DisplayName',strcat('lambda_{max} * ||x^* - x_k||_A'));
% plot(x_axis, lower_bound(1:last_iter+1),'-*','DisplayName',strcat('lambda_{min} * ||x^* - x_k||_A'));
% legend, box on;
% title(Problem.name);
% xlabel('Number of Iterations');
% ylabel('Log_{10} of relative norm');
% ylim(gca,[1e-9 1e+10]);
% set(gca, ...
%     'YScale', 'log', ...
%     'FontSize', 16);
% hold off;

x_axis = 0:last_iter;

figure;
hold on, grid on;
plot(x_axis, hist_relres_2(1:last_iter+1),'-*','DisplayName',strcat('||r_k||_2/||b||_2'));
plot(x_axis, hist_relerr_A(1:last_iter+1),'-*','DisplayName',strcat('||e_k||_A/||x_{true}||_A'));
plot(x_axis, upper_bound(1:last_iter+1),'-*','DisplayName',strcat('上限'));
% plot(x_axis, lower_bound(1:last_iter+1),'-*','DisplayName',strcat('lambda_{min} * ||x^* - x_k||_A'));
legend, box on;
title(Problem.name);
xlabel('Number of Iterations');
ylabel('Log_{10} of relative norm');
ylim(gca,[1e-14 1e+5]);
set(gca, ...
    'YScale', 'log', ...
    'FontSize', 16);
hold off;
