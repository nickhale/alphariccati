% Code to generate Figure 3 of Approximate solutions to a nonlinear
% functional differential equation
%
% Nick Hale, Stellenbosch University, Jan 2024

clc, close all

N = 700;
M = 6;
dom = [0 40];
fs = 14;

cols = colororder;
set(0, 'DefaultLineLinewidth', 3)

% Load ctored values C:
load C

for n = 1:3
    disp(['n = ' int2str(n)])
    alpn = 2^(1/n);

    % Alpha and t grids for output:
    alpp = alpn + [1;-1]'*10.^(-n)  ;
    tt = linspace(0, 10, 1000);

    % Laguerre discretisation:
    bet = 2;
    [t, D, wk, vk] = mylagdif(N, bet);
    % Remove node at t = 0:
    t = t(2:end); D = D(2:end,2:end);
    vk(1) = []; vk = vk.*t; wk(1) = [];

    % Remove points that will not contribute:
    idx = (1:M*sqrt(N))';
    t = t(idx);
    vk = vk(idx);
    wk = wk(idx);
    D = D(idx, idx);
    I = eye(size(D));
    
    % Compute characteristic functions using the infinite series:
    E = exp(-t);
    ck = 2/(1-alpn);
    dE = inf;
    k = 0;
    while ( norm(dE) > 1e-16 )
        k = k+1;
        dE = ck*exp(-(alpn^k)*t);
        E = E + dE;
        ck = 2*ck/(1-alpn^(k+1));
        if ( k > 1000 ), break, end
    end

    % Interpolate to a finer grid for plotting:
    EE = exp(-bet*tt/2).*bary(tt, exp(bet*t/2).*E, t, vk);

    % Compute the appropriate scaling
    Cn = C(n);

    % Barycentric resmapling matrix:
    tau = alpn*t;
    P = exp(-bet*tau/2).*barymat(tau, t, vk).*exp(bet*t'/2);

    % Initialise stored values:
    V = [];
    for alp = alpp

        ep = alp - alpn
        v = Cn*E*ep;

        % Barycentyric matrices:
        tau = alp*t;
        P = exp(-bet*tau/2).*barymat(tau, t, vk).*exp(bet*t'/2);
        E0 = barymat(0, t, vk).*exp(bet*t'/2);

        disp('Newton iterates:')
        for k = 1:30
            Pv = P*v;
            F = D*v + v - 2*Pv - (P*v).^2;
            J = D + I - 2*P - 2*P*diag(v);

            dv = J\F;
            disp([norm(dv, inf), norm(F, inf)])
            v = v - dv;

            if ( norm(dv) < 1e-10 ), break, end
        end

        % Interpolate to Chebyshev grid for smoother plots):
        vv = exp(-bet*tt/2).*bary(tt, exp(bet*t/2).*v, t, vk);

        figure(n)
        epstr = num2str(abs(ep));
        if ( ep > 0), epstr = ['$+' epstr '$']; else, epstr = ['$-' epstr '$']; end
        if ( alp < alpn ), col = cols(1,:); else, col = cols(2,:); end
        h1 = plot(t, v, 'color', col, 'displayname', epstr); hold on
        h2 = plot(t, ep*Cn*E, ':', 'HandleVisibility','off'); shg
        set(h2, 'color', get(h1, 'color'))
        title(['$\alpha_ ' int2str(n) ' = 2^{1/', int2str(n), ...
            '}, \ \varepsilon =  \pm',  num2str(abs(ep)), '$'], 'interpreter', 'latex')
        % cols(1,:) = [];
        xlabel('$t$', 'interp', 'latex', 'FontSize', fs)
        ylabel('$v$', 'interp', 'latex', 'FontSize', fs)
        l = legend;
        set(l, 'interp', 'latex');
        axis tight
        xlim([0 10]), shg
        grid on
        if ( n == 1 )
            title(['$\alpha_1 = 2, \ \varepsilon =  \pm',  ...
                num2str(abs(ep)), '$'], 'interpreter', 'latex')
        end
        yl = max(abs(ylim));
        ylim([-1 1]*yl*(1+eps));
        set(gca, 'FontSize', fs);
    end
    hold off

    switch n
        case 1
            ylim(.16501*[-1 1])
        case 2
            ylim(.1001*[-1 1])
        % case 3
            % ylim(.02501*[-1 1])
    end

    print('-depsc', ['fig3_pert_sol_' int2str(n)])
    
end

alignfigs % https://github.com/nickhale/alignfigs
