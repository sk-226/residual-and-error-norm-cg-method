% ================================ %
% 連立一次方程式に対する共役勾配法 %
% ================================ %

% SuiteSparse Matrix Collection: (https://sparse.tamu.edu/)
load('nos5.mat');
% load('LFAT5.mat');   % condition number = 2.374510e+09

% 係数行列の情報
%（SuiteSparse Matrix Collectionからダウンロードしたもの）
fprintf('========================== \n');
fprintf('Coefficient Matrix A.\n');
fprintf('========================== \n');
% A = Problem.A;
A = sparse(Problem.A);
[n, ~] = size(A);
fprintf('Name: %s\n', Problem.name);
fprintf('Kind: %s\n', Problem.kind);
fprintf('Size: %d×%d\n', n, n);
fprintf('Nnz : %d\n', nnz(A));
fprintf('\n');

% 真の解の設定
%（成分をすべて1のベクトルとする）
fprintf('========================== \n');
fprintf('Exact Solution x_true. \n');
fprintf('========================== \n');

x_true = ones(n, 1);

fprintf('x_true := [1,1,...,1]^T\n');
fprintf('\n');

% 右辺項の設定
%（係数行列 * 真の解 とする）
fprintf('========================== \n');
fprintf('Right-Hand Side Vector b. \n');
fprintf('========================== \n');

b = A * x_true; % 右辺項

fprintf('b := A * x_true\n');
fprintf('\n');

% 初期値の設定
%（成分をすべて0のベクトルとする）
fprintf('========================== \n');
fprintf('Initial Guess x_0. \n');
fprintf('========================== \n');

x = zeros(n, 1);    % 初期値 x_0

fprintf('x_0 := [0,0,...,0]^T\n');
fprintf('\n');

% 許容誤差と最大反復回数の設定
fprintf('========================== \n');
fprintf('Computational Conditions. \n');
fprintf('========================== \n');

eps = 10^(-12); % 許容誤差
% max_iter = 2 * n;
max_iter = 2 * n;

fprintf('Tolerance: Eps. := %.1e\n', eps);
fprintf('Maximum # Iter. := %d\n', max_iter);
fprintf('\n');

% CG法の実行
fprintf('************************** \n');
fprintf('Performing the CG Method. \n');
fprintf('************************** \n');

% 相対残差2ノルム，相対誤差2ノルム，相対誤差Aノルムの保存先
hist_relres_2 = zeros(max_iter + 1, 1);
hist_relerr_2 = zeros(max_iter + 1, 1);
hist_relerr_A = zeros(max_iter + 1, 1);

% 計測開始（収束履歴の記録も含まれることに注意）
tic;

% 右辺項bと真の解の2ノルム，および真の解のAノルムを保存
norm2_b = norm(b);
norm2_x_true = norm(x_true);
normA_x_true = sqrt(x_true' * A * x_true);

% 初期残差の設定
r = b; % 残差r = b - Ax (x_0=0)

% 初期の相対残差2ノルム，相対誤差2ノルム，相対誤差Aノルムの保存
hist_relres_2(1) = norm(r) / norm2_b;
err = x_true - x;
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

    % 相対残差2ノルム，相対誤差2ノルム，相対誤差Aノルムの保存
    hist_relres_2(iter+1) = norm(r) / norm2_b;
    err = x_true - x;
    hist_relerr_2(iter+1) = norm(err) / norm2_x_true;
    hist_relerr_A(iter+1) = sqrt(err' * A * err) / normA_x_true;

    % 収束しているか判定
    if hist_relres_2(iter+1) < eps
        fprintf('... Convergence!\n');
        fprintf('\n');
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
end

% 計算結果の表示
fprintf('========================== \n');
fprintf('Numerical Results. \n');
fprintf('========================== \n');
fprintf('# Iter.: %d\n', iter);
fprintf('Time[s]: %.3f\n', time);
fprintf('Relres_2norm = %.2e\n', hist_relres_2(iter+1));
fprintf('Relerr_2norm = %.2e\n', hist_relerr_2(iter+1));
fprintf('Relerr_Anorm = %.2e\n', hist_relerr_A(iter+1));
fprintf('========================== \n');
fprintf('\n');

% 収束履歴の表示
x_axis = 0:max_iter;
hold on, grid on;
plot(x_axis,hist_relres_2,'-*','DisplayName',strcat('||r_k||_2/||b||_2'));
plot(x_axis,hist_relerr_2,'-*','DisplayName',strcat('||e_k||_2/||x_{true}||_2'));
plot(x_axis,hist_relerr_A,'-*','DisplayName',strcat('||e_k||_A/||x_{true}||_A'));
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
