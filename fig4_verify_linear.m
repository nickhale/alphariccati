% Code to generate Figure 4 of Approximate solutions to a nonlinear
% functional differential equation
%
% Nick Hale, Stellenbosch University, Jan 2024

% This will take some time (because of the extended precision computation
% in the series solution), so I suggest you go make some coffee and begin
% planning your next D&D campaign.

clc, close all

NN = 5:1:700;
MM = 4:7;

fs = 14;
cols = colororder;
set(0, 'DefaultLineLinewidth', 2)

% Use the symbolic toolbox to compute the unstable series solution:
digits(100)

% Initialise storage:
err = []; sx = [];

for n = [1 3 5]
    alpn = 2^(1/n);

    for M = MM

        for N = NN
            % disp([M, N])

            % Laguerre discretisation:
            bet = 2;
            [t, D, wk, vk] = mylagdif(N, bet);
            % Remove node at t = 0:
            t = t(2:end); D = D(2:end,2:end);
            vk(1) = []; vk = vk.*t; wk(1) = [];
        
            % Remove points that will not contribute:
            % idx = (1:M*sqrt(N))';
            idx = 1:min(ceil(M*sqrt(N)), N-1);
            t = t(idx);
            vk = vk(idx);
            wk = wk(idx);
            D = D(idx, idx);
            I = eye(size(D));

            % Scale the quadrature weights to remove weight function:
            wk = wk.*exp(bet*t.');

            % Barycentric matrices:
            tau = alpn*t;
            P = exp(-bet*tau/2).*barymat(tau, t, vk).*exp(bet*t'/2);

            % Linearised problem:
            J = D + I - 2*P;
            [U, S, V] = svd(J);
            v = V(:,end);

            scl = sqrt(wk*(v.^2));
            v = v/scl;
            [~, idx] = max(abs(v));
            v = -v.*sign(v(idx));

            % % Double precision:
            % E = exp(-t);
            % ck = 2/(1-alpn);
            % dE = inf;
            % k = 0;
            % while ( norm(dE) > 1e-16 )
            %     k = k+1;
            %     dE = ck*exp(-(alpn^k)*t);
            %     E = E + dE;
            %     ck = 2*ck/(1-alpn^(k+1));
            %     if ( k > 100 ), break, end
            % end
            % disp([M, N, k])
            % scl = sqrt(wk*(E.^2));
            % E = (-1)^(n+1)*E;
            % E = E/scl;
            % u = E;

            % Extended precision:
            t_ = sym(t, 'f');
            alpn_ = sym(alpn, 'f');
            E_ = exp(-t_);
            ck_ = 2/(1-alpn);
            dE_ = inf;
            k = 0;
            while ( norm(dE_) > 1e-32 )
                k = k+1;
                dE_ = ck_*exp(-(alpn_^k)*t_);
                E_ = E_ + dE_;
                ck_ = 2*ck_/(1-alpn_^(k+1));
                if ( k > 1000 ), break, end
            end
            disp([M, N, k])
            E_ = double(E_);
            scl = sqrt(wk*(E_.^2));
            E_ = -(-1)^(n)*E_;
            u = E_/scl;

            sz(N) = length(t);
            err(N) = min(norm(u-v, inf), norm(u+v, inf));

        end

        figure(n)
        ss = sz(NN);
        semilogy(ss, err(NN)), shg
        hold on

    end

end

%%

% Add axes labels, legends, etc:
for k = 1:n
    figure(k)
    title(['$\alpha_' int2str(k) ' = 2^{1/' int2str(k) '}$'], 'interp', 'latex');
    legend('$M = 4$', '$M = 5$', '$M = 6$', '$M = 7$', 'interpreter', 'latex')
    grid on
    ylim([1e-14, 1e0])
    xlim([0 200])
    xlabel('$\lceil M\sqrt{N}\rceil$', 'interp', 'latex')
    ylabel('Error', 'interp', 'latex')
    print('-depsc', ['fig4_verifyLinear_' int2str(k)])
    set(gca, 'fontsize', fs)
    set(gca, 'ytick', 10.^(-12:4:0))
end
figure(1)
title('$\alpha_1 = 2$', 'interp', 'latex');
print('-depsc', ['figs/fig4_verifyLinear_' int2str(1)])
alignfigs


