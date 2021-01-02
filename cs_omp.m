function x=cs_omp(y,phi)
% y = phi*x
% phi：观测矩阵
% N：原信号长度

N = size(phi,2);
x = zeros(1,N);
A_t=[]; % 基矩阵
r_t=y; % 初始残差
inds = [];

for t=1:N
    % 计算每个列向量贡献
    coeff = abs(phi'*r_t);
    % 找到贡献最大的基
    [~,ind] = max(coeff);
    A_t = [A_t, phi(:,ind)];
    % 对相应的基去除
    phi(:,ind) = 0;
    % 最小二乘
    x_tmp = inv((A_t'*A_t))*A_t'*y;
    r_t = y-A_t*x_tmp;
    inds = [inds, ind];
    % 残差为0时停止迭代
    if norm(r_t) < 1e-9
        break
    end
end
x(inds) = x_tmp;