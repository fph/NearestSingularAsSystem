function [AplusDelta, Delta, u, v] = nearest_singular_structured_dense_newtonlike_minres(A, P, uv0)
% Computes the nearest singular matrix to A with a specified structure.
%
% uv0 is an initial value (optional), a vector with length sum(size(A)).
%
% P = structure (mn x p)

% uncomment for sanity checks of the formulas
% checkgradient
% checkhessian

[m, n] = size(A);

beta = norm(A,'fro');

if not(exist('uv0', 'var')) || isempty(uv0)
%    [V,D,W] = eig(full(A)); [~,ind] = min(abs(diag(D))); uv0 = [W(:,ind); D(ind,ind)*V(:,ind)];
%    fprintf('Computed initial value from eigenvalues of A\n');
    % [V,D,W] = svd(full(A)); uv0 = [V(:,end); D(end,end)*W(:,end)];
    [v, d, w] = svds(A, 1, 'smallest');
    uv0 = [d*v; w];
    fprintf('Computed initial value from singular values of A\n');
    % uv0 = randn(2*n,1);
end

u = uv0(1:m);
v = uv0(m+1:end);


maxit = 10;

for k = 1:maxit
    matop = @(dudv) H(u, v, dudv(1:m), dudv(m+1:end));
    rhs = G(u,v);

    if norm(rhs) / norm(A,1) < 1e-14
        fprintf('norm(rhs)=%g, frobnorm(Delta)=%g\n', norm(rhs), norm(Delta,'fro'));
        break
    end

    
    fprintf('*** k = %d \n', k);

    if norm(rhs) == 0
        fprintf('System already solved exactly, cannot improve the solution anymore\n')
        break
    end

    if isreal(A) && isreal(u) && isreal(v)
        duv = -minres(matop, rhs, 1e-2, m+n);
    else
        % since our operation is R-linear, we need to separate out
        % real and imaginary part to treat it like a linear system
        duvreal = -minres(@(x) c2r(matop(r2c(x))), c2r(rhs), 1e-2, 2*(m+n));
        duv = r2c(duvreal)
    end


    fprintf('System solved\n');

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
        rhsnew = G(unew, vnew);

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

if k == maxit
    fprintf('Maximum number of iterations reached\n')
end

vecDelta = P*(P' * kron(conj(v),u));
Delta = reshape(vecDelta, size(A));
AplusDelta = A + Delta;

    function c = r2c(r)
        c = r(1:end/2) + 1i * r(end/2+1:end);
    end
    function r = c2r(c)
        r = [re(c); im(c)];
    end

    function f = F(u,v)
        vecDelta = P*(P' * kron(conj(v),u));
        f = 1/2*norm(A(:) + vecDelta, 'fro')^2 + beta/4*(norm(v)^2-1)^2 ;
    end
    function g = G(u,v)
        vecDelta = P*(P' * kron(conj(v),u));
        Delta = reshape(vecDelta, size(A));
        g = [(A+Delta)*v; (A+Delta)'*u + beta*(norm(v)^2-1)*v];
    end
    function hd = H(u,v, du,dv)
        [m,n] = size(A);
        vecDelta = P*(P' * kron(conj(v),u));
        Delta = reshape(vecDelta, size(A));
        M = kron(transpose(v), eye(m))*P;
        N = kron(eye(n), transpose(u))*conj(P);
        hd = [M*(M'*du) + (A+Delta)*dv + M*conj(N'*dv);
              (A+Delta)'*du + N*conj(M'*du) + N*(N'*dv) + beta*(v'*v-1)*dv + beta*v*2*real(v'*dv)];
    end

    function checkgradient
        [m,n] = size(A);
        u = randn(m,1) + 1i*randn(m,1); v = randn(n,1) + 1i*randn(n,1);
        %u = randn(m,1); v = randn(n,1); 'TODO:testing only real gradient'
        du = randn(m,1) + 1i*randn(m,1); dv = randn(n,1) + 1i*randn(n,1);
        %du = randn(m,1); dv = randn(n,1); 'TODO:testing only real gradient'
        err = [];
        epsilon = logspace(-1, -10, 50);
        for i = 1:length(epsilon)
            err(i) = abs(F(u+epsilon(i)*du,v+epsilon(i)*dv) - F(u,v) - epsilon(i)*real(G(u,v)' * [du;dv]));
        end
        loglog(epsilon, err, epsilon, epsilon.^2)
        title('Gradient check: are the lines parallel?')
    end
    function checkhessian
        [m,n] = size(A);
        u = randn(m,1) + 1i*randn(m,1); v = randn(n,1) + 1i*randn(n,1);        
        %u = randn(m,1); v = randn(n,1); 'TODO:testing only real Hessian'
        du = randn(m,1) + 1i*randn(m,1); dv = randn(n,1) + 1i*randn(n,1);
        %du = randn(m,1); dv = randn(n,1); 'TODO:testing only real Hessian'
        err = [];
        epsilon = logspace(-1, -6, 50);
        for i = 1:length(epsilon)
            err(i) = abs(F(u+epsilon(i)*du,v+epsilon(i)*dv) - F(u,v) - epsilon(i)*real(G(u,v)' * [du;dv]) - 0.5*epsilon(i)^2*real([du;dv]'*H(u,v,du,dv)));
        end
        loglog(epsilon, err, epsilon, epsilon.^3)
        title('Hessian check: are the lines parallel?')
    end

end
