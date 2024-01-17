% Code to generate Figure 2 of Approximate solutions to a nonlinear
% functional differential equation
%
% Nick Hale, Stellenbosch University, Jan 2024

clc, close all

% Discretisation parameters
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
    alpp = alpn + [+1;-1]'*10.^(-n)  ;
    tt = linspace(0, 20, 1000);

    % Laguerre discretisation:
    bet = 2;

    I = eye(N);
    [t, ~, vk] = lagpts(N-1, 'glr'); % www.chebfun.org
    t = t/bet;

    % Remove points that will not contribute:
    idx = (1:M*sqrt(N))';
    t = t(idx);
    vk = vk(idx);

    % Compute characteristic functions using the infinite series:
    E = exp(-t);
    Ep = 0;
    ck = 2/(1-alpn);
    dE = inf;
    k = 0;
    while ( norm(dE) > 1e-16 )
        k = k+1;
        dE = ck*exp(-(alpn^k)*t);
        dEp = (1-alpn^(k))*dE;
        E = E + dE;
        Ep = Ep + dEp;
        ck = 2*ck/(1-alpn^(k+1));
        if ( k > 1000 ), break, end
    end
    Ep = Ep - E;

    % Compute the appropriate scaling
    Cn = C(n);

    % Interpolate to finer grid for plotting:
    EE = exp(-bet*tt/2).*bary(tt, exp(bet*t/2).*E, t, vk);
    EEp = exp(-bet*tt/2).*bary(tt, exp(bet*t/2).*Ep, t, vk);

    % Plot:
    figure
    rr = -2*Cn*EEp.*tt/alpn - (Cn*EE).^2;
    h1 = plot(tt/alpn, rr); shg
    title(['$\alpha_ ' int2str(n) ' = 2^{1/', int2str(n), ...
        '}$'], 'interpreter', 'latex')
    if ( n == 1 )
        title('$\alpha_1 = 2$', 'interpreter', 'latex');
    end
    xlabel('$t$', 'interp', 'latex', 'FontSize', fs)
    ylabel('$r$', 'interp', 'latex', 'FontSize', fs)
    axis tight
    xlim([0 10]), shg
    grid on
    set(gca, 'FontSize', fs);

    switch n
        case 1
            ylim([-3 1.5])
        case 2
            ylim([-70 30])
        case 3
            ylim([-500 200])
    end

    print('-depsc', ['fig2_remainder_' int2str(n)])
    hold off
end


