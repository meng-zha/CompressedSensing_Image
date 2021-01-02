function x=cs_bayes(y,phi)
% Bayes Compressed Sensing
% refer to https://github.com/navaneethakrishna/Bayesian-Compressive-Sensing

[M,N] = size(phi);
% init the sigma, mu, alpha, alpha0
alpha0 = 100;
alpha = 0.01*ones(1,N);

iter = 100;
non_zeros_inds = 1:N;
likelihood = [-99];
for i = 1:iter
    % estimate sigma, mu
    sigma = inv(alpha0*(phi'*phi)+diag(alpha));
    mu = alpha0*sigma*phi'*y;
    
    n = size(sigma,1);
    % update alpha, alpha0
    gamma = ones(1,n)-diag(sigma)'*diag(alpha);
    alpha = gamma./mu.^2;
    alpha0 = (N-sum(gamma))/norm(y-phi*mu)^2;
    
    % remove the zero alpha
    ind = gamma >= 1e-10;
    alpha = alpha(ind);
    phi = phi(:,ind);
    sigma = sigma(ind,ind);
    mu = mu(ind);
    non_zeros_inds = non_zeros_inds(ind);
    
    % cal the likelihood
    C = 1/alpha0*eye(M)+phi*inv(diag(alpha))*phi';
    likelihood = [likelihood, log(det(C))+y'*inv(C)*y]; 
    
    % terminate
    if isempty(ind) || abs(likelihood(i+1)-likelihood(i))<1e-5
        break;
    end
end
x = zeros(1,N);
x(non_zeros_inds) = mu;