function [AplusDelta, Delta u v] = nearest_singular_sparse_newtonlike(A, P, uv0, opts)
arguments
    A {mustBeNumeric}
    P {mustBeNumericOrLogical} = []
    uv0 {mustBeNumeric} = []
    opts.DirectSolve logical = false
    opts.DirectSvd logical = false
    opts.NInitialValues = 1
    opts.maxit {mustBeInteger} = 10
    opts.tol double = 0
    opts.beta double = norm(A, 'fro')  % this seems a reasonably-scaled choice
    opts.minres_tolerance = 1e-2;
end

% Computes the nearest singular matrix to A (m x n) with a specified sparsity
% pattern.
%
% uv0 contains initial values in the format [u;v] ((m+n) x 1)
% If uv0 contains multiple columns, the method is run with several initial
% values and the best Delta is returned.
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
    if opts.DirectSvd
        [V,D,W] = svd(full(A));
        niv = opts.NInitialValues;
        uv0 = [V(:,end-niv+1:end)*D(end-niv+1:end,end-niv+1:end); W(:,end-niv+1:end)];
    else
        [V,D,W] = svds(A, opts.NInitialValues, 'smallest');
        uv0 = [V*D; W];
    end
    fprintf('Computed initial values from singular values of A\n');
end

bestdist = inf;
bestu = [];
bestv = [];

for iv = 1:size(uv0, 2)

    if size(uv0, 2) > 1
        fprintf('# Trying initial value #%d\n', iv);
    end

    u = uv0(1:m, iv);
    v = uv0(m+1:end, iv);

    for k = 1:opts.maxit
        % assembling linear system

        k1 = P * (v .* conj(v));
        k2 = P' * (u .* conj(u));

        normv2m1 = v'*v-1;

        rhs = [k1.*u + A*v; A'*u + k2.*v + beta*normv2m1*v];
        if norm(rhs) / norm(A,1) < 1e-14
            fprintf('norm(rhs)=%g\n', norm(rhs));
            break
        end

        fprintf('## k = %d \n', k);
        if opts.DirectSolve
            Delta = u .* (v' .* P);
            fullmat = [diag(k1) A+2*Delta; (A+2*Delta)' diag(k2)+beta*(v*v'*2+normv2m1*eye(m))];
            duv = -fullmat \ rhs;
        else
            matop = @(uv) [
                k1.*uv(1:m) +                                         A*uv(m+1:end) + (P * (uv(m+1:end) .* conj(v))) .* (2*u);
                A'*uv(1:m) + (P' * (uv(1:m) .* conj(u))) .* (2*v) +   k2.*uv(m+1:end) + beta*(v*(v'*uv(m+1:end))*2 + normv2m1*uv(m+1:end));
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

    % equivalent to: Delta = u .* (v' .* P); normDelta = norm(Delta, 'fro')
    normDelta = sqrt((u.*conj(u))' * (P * (v.*conj(v))))
    
    if normDelta < bestdist
        bestdist = normDelta;
        bestu = u;
        bestv = v;
    end

end

u = bestu;
v = bestv;
Delta = u .* (v' .* P);
AplusDelta = A + Delta;
