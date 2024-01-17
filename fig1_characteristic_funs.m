% Code to generate Figure 1 of Approximate solutions to a nonlinear
% functional differential equation
%
% Nick Hale, Stellenbosch University, Jan 2024

clc, close all

figure(1)
set(gcf, 'position', [680 1538 1034 420])

fs = 14;
N = 700; M = 6;

cols = colororder;
set(0, 'DefaultLineLinewidth', 3)

% Fine grid for plotting
tt = linspace(0,15,1000);

% Intialise storage:
str = {};

for n = 1:6
    alpn = 2^(1/n);

    % Laguerre discretisation (for normalising):
    bet = 2;
    [t, wk, vk] = lagpts(N-1, 'glr'); % www.chebfun.org
    t = t/bet; wk = wk/bet;
    
    % Remove points that will not contribute:
    idx = (1:M*sqrt(N))';
    t = t(idx);
    vk = vk(idx);
    wk = wk(idx);
    
    % Scale the quadrature weights to remove weight function:
    wk = wk.*exp(bet*t.');
    
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

    % Interpolate to a finer grid for plotting using barycentric interp:
    EE = exp(-bet*tt/2).*bary(tt, exp(bet*t/2).*E, t, vk);

    % Scale for plotting (L_2 norm)
    scl = (-1)^(n+1)*sqrt(wk*(E.^2));
    v0 = E/scl;
    vv0 = exp(-bet*tt/2).*bary(tt, exp(bet*t/2).*v0, t, vk);

    % Plot
    figure(1)
    plot(tt, vv0)
    xlim([0 15]), ylim([-1 1]),     
    hold on, grid on

    % Store legend text:
    str{n} = ['$n$ = ' int2str(n)];

end

figure(1)
% Add axis labels and legend:
xlabel('$t$', 'interpreter', 'latex', 'fontsize', fs), 
ylabel('$u$', 'interpreter', 'latex', 'fontsize', fs);
legend(str{:}, 'interpreter', 'latex', 'fontsize', fs)
set(gca, 'fontsize', fs);

print -depsc fig1_characteristic_funs