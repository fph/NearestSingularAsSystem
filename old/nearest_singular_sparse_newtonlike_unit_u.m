function [AplusDelta, Delta u v] = nearest_singular_sparse_newtonlike(A, P, uv0)
%
% Computes the nearest singular matrix to A with a specified sparsity
% pattern.
%
% uv0 is an initial value (optional), a vector with length sum(size(A)).
%
% P = nonzero pattern of the permutation (default=nonzero pattern of A)
%
% Works only in the real case!

beta = norm(A,'fro'); % this seems a reasonably-scaled choice

if not(exist('P', 'var')) || isempty(P)
    P = double(A ~= 0);
end

n = size(A, 2);

if not(exist('uv0', 'var')) || isempty(uv0)
%    [V,D,W] = eig(full(A)); [~,ind] = min(abs(diag(D))); uv0 = [W(:,ind); D(ind,ind)*V(:,ind)];
%    fprintf('Computed initial value from eigenvalues of A\n');
    [V,D,W] = svd(full(A)); uv0 = [V(:,end); D(end,end)*W(:,end)];
    fprintf('Computed initial value from singular values of A\n');
    % uv0 = randn(2*n,1);
end

u = uv0(1:n);
v = uv0(n+1:end);

for k = 1:50
    % assembling linear system

    K1 = diag(P * (v .* conj(v)));
    K2 = diag(P' * (u .* conj(u)));

    Delta = u .* (v' .* P);
    normu2m1 = u'*u-1;

    rhs = [K1+beta*(normu2m1)*eye(size(K1)) A; A' K2]*[u;v];
    if norm(rhs) / norm(A,1) < 1e-14
        fprintf('norm(rhs)=%g, frobnorm(Delta)=%g\n', norm(rhs), norm(Delta,'fro'));
        break
    end
    mat = [K1+beta*(u*u'*2+normu2m1*eye(size(K1))) A+2*Delta; (A+2*Delta)' K2];

    fprintf('*** k = %d \n', k);
    % u
    % v
    % eigmat = eig(mat)
    Delta = u .* (v' .* P);
    normrhs = norm(rhs);
    if norm(rhs) == 0
        fprintf('System already solved exactly, cannot improve the solution anymore\n')
        break
    end
    
    duv = -mat \ rhs;

    if any(isnan(duv))
        fprintf('Singular matrix, cannot improve the solution anymore\n')
        break
    end

    du = duv(1:n);
    dv = duv(n+1:end);

    % A rudimentary line search

    alpha = 1;
    while(alpha > 1e-10)
        unew = u + alpha*du;
        vnew = v + alpha*dv;
        nrm = norm(unew);
        unew = unew/nrm;
        vnew = vnew*nrm;
        K1new = diag(P * (vnew .* conj(vnew)));
        K2new = diag(P' * (unew .* conj(unew)));
        rhsnew = [K1new A; A' K2new]*[unew;vnew];
        if norm(rhsnew) < norm(rhs)
            fprintf('line search: accepted alpha=%g; stopping=%g, norm(rhs)=%g\n', alpha,norm(rhsnew),norm(rhs));
            u = unew;
            v = vnew;
            break
        end
%        fprintf('line search: rejected alpha=%g; stopping=%g, norm(rhs)=%g\n', alpha,norm(rhsnew),norm(rhs));
        alpha = alpha / 2;
    end
    if alpha <= 1e-10
        fprintf('Line search failed, cannot improve the solution anymore\n')
    end
end

if k == 50
    fprintf('Maximum number of iterations reached\n')
end
Delta = u .* (v' .* P);
AplusDelta = A + Delta;
