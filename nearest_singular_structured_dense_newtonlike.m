function [AplusDelta, Delta, u, v] = nearest_singular_structured_dense_newtonlike(A, P, uv0, opts)
arguments
    A {mustBeNumeric}
    P {mustBeNumericOrLogical} = []
    uv0 {mustBeNumeric} = []
    opts.DirectSolve logical = (length(A)<500)
    opts.DirectSvd logical = (length(A)<500)
    opts.NInitialValues = 1
    opts.maxit {mustBeInteger} = 10
    opts.tol double = 1e-14
    opts.beta double = norm(A, 'fro')  % this seems a reasonably-scaled choice
    opts.minres_tolerance = 1e-2;
    opts.renormalization = true;
    opts.checkgradhess logical = false
    opts.armijo_c = 1e-4;
    opts.alpha0 = 1;
    opts.alpha_decrease = 0.5;
end
% Compute the nearest singular matrix to A with a specified structure.
%
% uv0 contains initial values in the format [u;v] ((m+n) x 1)
% If uv0 contains multiple columns, the method is run with several initial
% values and the best Delta is returned.
%
% P = structure (mn x p) or (m x n x p), like the one produced by
% autobasis() in RiemannOracle
%
% if opts.checkgradhess = true, checks gradient and hessian instead
% of computing things.

[m, n] = size(A);

beta = opts.beta;

if isempty(P)
    % this uses a function from RiemannOracle, which we assume is in the
    % path
    P = autobasis(A);
end

if ndims(P) == 3
    [mp, np, pp] = size(P);
    assert(isequal(size(A), [mp,np]));
    P = reshape(P, [mp*np, pp]);
end

if opts.checkgradhess
    figure(1)
    checkgradient
    figure(2)
    checkhessian
    return
end

if isempty(uv0)
    niv = opts.NInitialValues;
    if opts.DirectSvd
        [V,D,W] = svd(full(A));
    else
        [V,D,W] = svds(A, opts.NInitialValues, 'smallest');
    end
    if opts.renormalization
        % renormalized initial values
        uv0 = nan(m+n, niv);
        for k = 1:niv
            u = V(:, end-niv+k); sigma = D(end-niv+k, end-niv+k); v = W(:,end-niv+k);
            vecGdelta = (P' * kron(conj(v),u)); % the projection is vec(G) = P*P'*vec(u*v'), but we can ignore the outer P if we only want the norm
            frobnormdelta_squared = norm(vecGdelta)^2;
            uv0(:,k) = [-u * sigma / frobnormdelta_squared; v];
        end
    else
        uv0 = [V(:,end-niv+1:end)*D(end-niv+1:end,end-niv+1:end); W(:,end-niv+1:end)];
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

%    vecDelta = P*(P' * kron(conj(v),u));
%    Delta = reshape(vecDelta, size(A));
%    min(svd(full(A+Delta)))


    for k = 1:opts.maxit
        Hbetaop = @(dudv) H(u, v, dudv(1:m), dudv(m+1:end));
        Gbeta = G(u,v);

        if norm(Gbeta) / norm(A,1) < opts.tol
            fprintf('norm(rhs)=%g, frobnorm(Delta)=%g\n', norm(Gbeta), norm(Delta,'fro'));
            break
        end


        fprintf('*** k = %d \n', k);

        if norm(Gbeta) == 0
            fprintf('System already solved exactly, cannot improve the solution anymore\n')
            break
        end

        if isreal(A) && isreal(u) && isreal(v)
            %
            %cond(fullmat)
            if opts.DirectSolve
                % warning: very slow
                Hbeta = fullmat;
                duv = -Hbeta \ Gbeta; % condition_number =cond(fullmat)
                HtG = Hbeta' * Gbeta;
            else
                HtG = Hbetaop(Gbeta); % H is symmetric, so we do not need to transpose
                duv = -minres(Hbetaop, Gbeta, opts.minres_tolerance, m+n);
            end
        else
            % since our operation is R-linear, we need to separate out
            % real and imaginary part to treat it like a linear system
            if opts.DirectSolve
                error('Not implemented')
            else
                duvreal = -minres(@(x) c2r(Hbetaop(r2c(x))), c2r(Gbeta), 1e-2, 2*(m+n));
                duv = r2c(duvreal);
            end
        end

        fprintf('System solved\n');

        if any(isnan(duv))
            fprintf('Singular matrix, cannot improve the solution anymore\n')
            break
        end


        du = duv(1:m);
        dv = duv(m+1:end);

        % Armijo line search

        alpha = opts.alpha0;
        while true
            unew = u + alpha*du;
            vnew = v + alpha*dv;
            nrm = norm(vnew);
            vnew = vnew/nrm;
            unew = unew*nrm;
            Gbetanew = G(unew, vnew);
            armijo_inc = 2*opts.armijo_c * alpha * HtG'*duv;
            if norm(Gbetanew)^2 < norm(Gbeta)^2 + armijo_inc;
                fprintf('line search: accepted alpha=%g; ||Gnew||=%g, ||Gold||=%g, armijo_inc=%g\n', alpha,norm(Gbetanew),norm(Gbeta),armijo_inc);
                u = unew;
                v = vnew;
                break
            end
            fprintf('line search: rejected alpha=%g; ||Gnew||=%g, ||Gold||=%g, armijo_inc=%g\n', alpha,norm(Gbetanew),norm(Gbeta),armijo_inc);
            alpha = alpha * opts.alpha_decrease;            
            if alpha <= 1e-10
                fprintf('Line search failed, cannot improve the solution anymore\n');
                break
            end
        end


    end

    if k == opts.maxit
        fprintf('Maximum number of iterations reached\n')
    end


    vecDelta = P*(P' * kron(conj(v),u));
    Delta = reshape(vecDelta, size(A));
    normDelta = norm(Delta, 'fro')
    if normDelta < bestdist
        bestdist = normDelta;
        bestu = u;
        bestv = v;
    end

end

u = bestu;
v = bestv;
vecDelta = P*(P' * kron(conj(v),u));
Delta = reshape(vecDelta, size(A));
AplusDelta = A + Delta;

    function c = r2c(r)
        % convert a real vector of length 2k [re; im] into re + 1i*im
        c = r(1:end/2) + 1i * r(end/2+1:end);
    end
    function r = c2r(c)
        % inverse of r2c
        r = [real(c); imag(c)];
    end

% function, gradient, Hessian matvec
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
        title('Gradient check: the lines must be parallel')
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
        title('Hessian check: the lines must be parallel')
    end
    function M = fullmat
        % compute the full matrix from matop(); useful for debugging Minres
        if isreal(A) && isreal(u) && isreal(v)
            M = eye(m+n);
            for s = 1:(m+n)
                M(:,s) = Hbetaop(M(:,s));
            end
        else
            error('Not implemented');
        end
    end
end
