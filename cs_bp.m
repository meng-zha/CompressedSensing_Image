function x=cs_bp(y,phi)
% minimize ||x||_1 s.t. y = phi*x
% equal to
% minimize 1'*v s.t. y=phi*x; -v<=x<=v

N = size(phi,2);
% 组合 x_aug = (x,v),变化对应的系数矩阵
f = ones(2*N,1);
f(1:N) = 0;
Aeq = [phi, zeros(size(phi))];
beq = y;
A = [-eye(N),-eye(N);eye(N),-eye(N)];
b = zeros(2*N,1);

options = optimoptions('linprog','Display','none');
% 使用linprog求解
x_aug = linprog(f,A,b,Aeq,beq,[],[],[],options);
x = x_aug(1:N);
