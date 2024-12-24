% SuiteSparse Matrix Collection: (https://sparse.tamu.edu/)
% load('nos5.mat');
load('LFAT5.mat');   % condition number = 2.374510e+09

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


% 固有値の計算
lambda = eig(A);
lambda_min = min(lambda);
lambda_max = max(lambda);

fprintf('lambda_min: %e\n', lambda_min);
fprintf('lambda_max: %e\n', lambda_max);

% eigs (大規模な場合)
lambda_max = eigs(A, 1, 'largestreal');
lambda_min =eigs(A, 1, 'smallestreal'); 

fprintf('lambda_min: %e\n', lambda_min);
fprintf('lambda_max: %e\n', lambda_max);


% eigs (大規模な場合), 絶対値
lambda_max = eigs(A, 1, 'largestabs');
lambda_min =eigs(A, 1, 'smallestabs'); 

fprintf('lambda_min: %e\n', lambda_min);
fprintf('lambda_max: %e\n', lambda_max);


