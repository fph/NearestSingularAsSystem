function [AplusDelta, Delta] = nearest_singular_sparse_newtonlike_minres(A, P, uv0)
%
% Computes the nearest singular matrix to A with a specified sparsity
% pattern.
%
% uv0 is an initial value (optional), a vector with length sum(size(A)).
%
% P = nonzero pattern of the permutation (default=nonzero pattern of A)
%
% This uses minres rather than backslash
% Works only in the real case!

beta = norm(A,'fro'); % this seems a reasonably-scaled choice

if not(exist('P', 'var')) || isempty(P)
    P = double(A ~= 0);
end

n = size(A, 2);

if not(exist('uv0', 'var')) || isempty(uv0)
%    [V,D,W] = eig(full(A)); [~,ind] = min(abs(diag(D))); uv0 = [W(:,ind); D(ind,ind)*V(:,ind)];
%    fprintf('Computed initial value from eigenvalues of A\n');
    % [V,D,W] = svd(full(A)); uv0 = [V(:,end); D(end,end)*W(:,end)];
    [v, d, w] = svds(A, 1, 'smallest');
    uv0 = [v; d*w];
    fprintf('Computed initial value from singular values of A\n');
    % uv0 = randn(2*n,1);
end

u = uv0(1:n);
v = uv0(n+1:end);

for k = 1:10
    % assembling linear system

    k1 = P * (v .* conj(v));
    k2 = P' * (u .* conj(u));

    Delta = u .* (v' .* P);
    normu2m1 = u'*u-1;

    % rhs = [k1+beta*normu2m1*eye(size(k1)) A; A' k2]*[u;v];
    rhs = [k1.*u + beta*normu2m1*u + A*v; A'*u + k2.*v];

    if norm(rhs) / norm(A,1) < 1e-14
        fprintf('norm(rhs)=%g, frobnorm(Delta)=%g\n', norm(rhs), norm(Delta,'fro'));
        break
    end
    %mat = [diag(k1)+beta*(u*u'*2+normu2m1*eye(length(k1))) A+2*Delta; (A+2*Delta)' diag(k2)];
    % [L, U] = ilu(mat);
    matop = @(uv) [
        k1.*uv(1:n) + beta*(u*(u'*uv(1:n))*2 + normu2m1*uv(1:n)) + (A+2*Delta)*uv(n+1:2*n);
        (A+2*Delta)'*uv(1:n) + k2.*uv(n+1:2*n);
        ];

    fprintf('*** k = %d \n', k);
    % u
    % v
    % eigmat = eig(mat)
    % Delta = u .* (v' .* P);
    % normrhs = norm(rhs);
    if norm(rhs) == 0
        fprintf('System already solved exactly, cannot improve the solution anymore\n')
        break
    end
    
    % duv = -mat \ rhs;
    % duv = -minres(matop, rhs, 1e-10, n);
    % duv = -gmres(matop, rhs, [], 1e-10, n);
    duv = -minres(matop, rhs, 1e-2, n);
    fprintf('TODO')
    fprintf('System solved');

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
        k1new = P * (vnew .* conj(vnew));
        k2new = P' * (unew .* conj(unew));
        % rhsnew = [k1new A; A' k2new]*[unew;vnew];
        rhsnew = [k1new.*unew + A*vnew; A'*unew+k2new.*vnew];
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
