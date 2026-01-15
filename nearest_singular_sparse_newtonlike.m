function [AplusDelta, Delta u v] = nearest_singular_sparse_newtonlike(A, P, uv0, opts)
arguments
    A {mustBeNumeric}
    P {mustBeNumericOrLogical} = []
    uv0 {mustBeNumeric} = []
    opts.DirectSolve logical = false
    opts.maxit {mustBeInteger} = 10
    opts.tol double = 0
    opts.beta double = norm(A, 'fro')  % this seems a reasonably-scaled choice
    opts.minres_tolerance = 1e-2;
end

% Computes the nearest singular matrix to A with a specified sparsity
% pattern.
%
% uv0 is an initial value (optional), a vector with length sum(size(A)).
%
% P = nonzero pattern of the permutation (default=nonzero pattern of A)
%
% Works only in the real case for now!

[m, n] = size(A);

assert(isreal(A) && isreal(uv0));
beta = opts.beta;

if isempty(P)
    P = double(A ~= 0);
end

if isempty(uv0)
    %    [V,D,W] = eig(full(A)); [~,ind] = min(abs(diag(D))); uv0 = [W(:,ind); D(ind,ind)*V(:,ind)];
    %    fprintf('Computed initial value from eigenvalues of A\n');
    [V,D,W] = svd(full(A)); uv0 = [D(end,end)*V(:,end); W(:,end)];
    fprintf('Computed initial value from singular values of A\n');
    % uv0 = randn(2*n,1);
end

u = uv0(1:n);
v = uv0(n+1:end);

for k = 1:opts.maxit
    % assembling linear system

    k1 = P * (v .* conj(v));
    k2 = P' * (u .* conj(u));

    Delta = u .* (v' .* P);
    normv2m1 = v'*v-1;

    rhs = [k1.*u + A*v; A'*u + k2.*v + beta*normv2m1*v];
    if norm(rhs) / norm(A,1) < 1e-14
        fprintf('norm(rhs)=%g, frobnorm(Delta)=%g\n', norm(rhs), norm(Delta,'fro'));
        break
    end



    fprintf('*** k = %d \n', k);
    if opts.DirectSolve
        fullmat = [diag(k1) A+2*Delta; (A+2*Delta)' diag(k2)+beta*(v*v'*2+normv2m1*eye(size(K2)))];
        duv = -fullmat \ rhs;
    else
        Delta = u .* (v' .* P);
        % TODO: avoid assembling Delta here
        matop = @(uv) [
        k1.*uv(1:m) + (A+2*Delta)*uv(m+1:end);
        (A+2*Delta)'*uv(1:m) + k2.*uv(m+1:end) + beta*(v*(v'*uv(m+1:end))*2 + normv2m1*uv(m+1:end));
        ];
        duv = -minres(matop, rhs, opts.minres_tolerance, m+n);
    end

    if any(isnan(duv))
        fprintf('Singular matrix, cannot improve the solution anymore\n')
        break
    end

    du = duv(1:m);
    dv = duv(m+1:end);

    % A rudimentary line search

    alpha = 1;
    while(alpha > 1e-10)
        unew = u + alpha*du;
        vnew = v + alpha*dv;
        nrm = norm(vnew);
        vnew = vnew/nrm;
        unew = unew*nrm;
        k1new = P * (vnew .* conj(vnew));
        k2new = P' * (unew .* conj(unew));
        rhsnew = [k1new.*unew + A*vnew; A'*unew + k2new.*vnew];
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

if k == opts.maxit
    fprintf('Maximum number of iterations reached\n')
end
Delta = u .* (v' .* P);
AplusDelta = A + Delta;
