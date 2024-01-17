% Code to generate animated version of Figure 5 of Approximate solutions to
% a nonlinear functional differential equation
%
% Nick Hale, Stellenbosch University, Jan 2024

% This code writes the frames to individual png files, which are then
% converted to mp4 using ffmpeg

set(0, 'DefaultLineLinewidth', 3)

ca, clc
col = colororder;
col = [col ; col];
fs = 14; 

figure(1)
subplot(1,2,1)
plot(NaN, NaN); hold on
ABC = {'A', 'B', 'C', 'D', 'E' ,'F', 'G'};
for k = 0:6
    plot((k+.5)*[1 1], [-1, 100], '--', 'color', .8*[1 1 1], 'linewidth', 1.5);
    text((k+.5), 5.15, ['$\mathsf{' ABC{k+1} '}$'], 'horizontalalign', 'center', 'fontsize', fs, 'interp', 'latex');
end
plot([0 20], [0 0], 'k');
xlim([0, 6.75]), ylim([-.1 5])
xlabel('$n$', 'interp', 'latex'); ylabel('$\|v\|_2^2$', 'interp', 'latex')
set(gca, 'fontSize', fs), grid on

% plotDecayingSoln();
nlist = .5:.05:6.5;
for k = 1:numel(nlist)
    figure(k+1);
    plotDecayingSoln(nlist(k));
    ylim([0 1.5])
end

for n = 1:6
    for direction = [1 -1]
        alp0 = 2^(1/n).*(1 + direction*1e-6);
        if ( direction  > 0 ), LS = ':'; else, LS = '-'; end
        [vv, aalp, mvec] = followpath(alp0, ...
            'direction', direction, ... 
            'stepinit', 1e-6, ...
            'plotting', false, ...
            'plotdata', {'color', col(n,:), 'LineStyle', LS}, ... ,
            'nlist', nlist);
        figure(1)
        plot(1./log2(aalp), mvec, 'color', col(n,:), 'linestyle', LS); 
        hold on
        idx = find(nlist == n);
        if (~isempty(idx))
            figure(idx+1);
            plot([0 20], [1 1], 'color', col(n,:), 'linestyle', LS); hold on
            title(['$n = ' num2str(n,'%2.2f') '$'], 'interp', 'latex');
            set(gca, 'fontsize', fs);
        end
    end
end


% This gets missed for some reason
figure(12)
xlabel('$t$', 'interp', 'latex');
ylabel('$u$', 'interp', 'latex');
set(gca, 'fontsize', fs);

%%

% Collate two plots into a single figure using subplots:
clc
for k = 1:numel(nlist)
    figure(999999), clf
    set(gcf, 'position', [214        1538        1026         420])
    ax = zeros(2,1);
    for i = 1:2
        ax(i)=subplot(1,2,i);
    end

    figure(1)
    if ( k > 1 ), delete(h1); end
    h1 = plot(nlist(k)*[1 1], [0 10], '--m', 'markersize', 3);
    
    figs = [1 k+1];
    for i = [1 2]
        figure(figs(i))
        h = get(gcf,'Children');
        newh = copyobj(h,999999);
        for j = 1:length(newh)
            posnewh = get(newh(j),'Position');
            possub  = get(ax(i),'Position');
            set(newh(j),'Position',...
                [possub(1) posnewh(2) possub(3) posnewh(4)]);
            [possub(1) posnewh(2) possub(3) posnewh(4)];
        end
        delete(ax(i));
    end
    figure(999999), drawnow, pause(.1)
    exportgraphics(figure(999999),['movie/frame_' num2str(k, '%3.3d') '.png'],'Resolution',600)
end 
delete(h1)


function [vv, aalp, mvec] = followpath(alp0, varargin)
% DEVELOPER NOTE - Adapted from Chebfun: www.chebfun.org

% Set default values
plotdata = {};
stepmax = 1e-2;               % Maximum steplength
stepmin = 1e-15;            % Mininum steplength
stepinit = [];            % Initial steplength
maxstepno = 1000;             % Maximum number of steps taken
alpmin = 2^(1/6.75);
alpmax = 2^(1/.3);
stopfun = @(v, alp) alp > alpmax || alp < alpmin;
measure = @(v, w) (w*v.^2);
fs = 14;

nlist = .5:.5:6.5;

% Parse VARARGIN
while ( ~isempty(varargin) ) % Go through all elements
    val = varargin{2};
    switch lower(varargin{1})
        case 'direction'
            direction = val;
        case 'measure'
            measure = val;
        case 'plotdata'
            plotdata = val;
        case 'maxstepno'
            maxstepno = val;
        case 'stepmax'
            stepmax = val;
        case 'stepmin'
            stepmin = val;
        case 'stepinit'
            stepinit = val;
        case 'stopfun'
            stopfun = val;
        case 'nlist'
            nlist = val;
    end

    % Throw away option name and argument and move on
    varargin(1:2) = [];
end

if ( isempty(stepinit) )
    stepinit = sqrt(stepmin*stepmax);
end

%%
% Compute the initial solution
n = round(1/log2(alp0));
alpn = 2^(1/n);

% Laguerre discretisation:
N = 700;
M = 6;
bet = 2;
[t, D, wk, vk] = mylagdif(N, bet);
% Remove node at t = 0:
t = t(2:end); D = D(2:end,2:end);
vk(1) = []; vk = vk.*t; wk(1) = [];
% Remove points that will not contribute:
idx = (1:M*sqrt(N))';
% idx = 1:min(ceil(M*sqrt(N)), N-1);
t = t(idx);
vk = vk(idx);
wk = wk(idx);
D = D(idx,idx);
I = eye(size(D));

% Scale the quadrature weights to remove weight function:
wk = (wk.*exp(bet*t).');

% Barycentric resmapling matrix at alp = alpn:
tau = alpn*t;
P = exp(-bet*tau/2).*barymat(tau, t, vk).*exp(bet*t'/2);

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

% Initial guess:
ep = alp0 - alpn;
data = load('C'); Cn = data.C(n); % Import C value from file.
v = Cn*E*ep;

% Barycentyric matrix for intitial alpha:
tau = alp0*t;
P = exp(-bet*tau/2).*barymat(tau, t, vk).*exp(bet*t'/2);

% Newton iteration
warning('off', 'MATLAB:nearlySingularMatrix');
for k = 1:10
    Pv = P*v;
    F = D*v + v - 2*Pv - (Pv).^2;
    J = D + I - 2*P - 2*P*diag(v);

    dv = J\F;
    v = v - dv;

    if ( norm(dv, inf) < 1e-10 ), break, end
end
warning('on', 'MATLAB:nearlySingularMatrix');

% Store initial values
vv = v;
aalp = alp0;
mvec = measure(v, wk);

% Set up initial tangent and tau we want to be orthogonal to.
dv_old = 0*v;
da_old = 1;

% Initialise:
retract = false;
counter = 1;
v_old = v;
alp_old = alp0;
sl = stepinit;

% Start the pseudo-arclength path following algorithm. It proceeds as follows:
%   1. Find a tangent direction.
%   2. Move in the direction of the tangent for the given steplength.
%   3. Compute the Newton correction.
%   4. If Newton was happy, accept the new point on the curve. If not, shrink
%      the steplength by a factor of 4 and go back to step 2.

while ( counter < maxstepno )

    % Find a tangent direction, but only if we were told by Newton not to
    % retract. At the start of the while loop, recall that RETRACT = false.
    if ( ~retract )
        % Compute the tangent
        [dv, da] = tangentBVP(v_old, alp_old, dv_old, da_old);

        % At the start, we need to ensure we're going in right direction:
        if ( counter == 1 )
            dv = direction*dv;
            da = direction*da;
        end

        % Move in the direction of the tangent
        v = v_old + sl*dv;
        alp = alp_old + sl*da;
    end

    % Find a Newton correction to get back on the solution curve:
    [v, alp, newtonIter, retract] = newtonBVP(v, alp, dv, da);

    % If the Newton correction algorithm told us we were trying to take too
    % long tangent steps, we decrease the steplength.
    if ( retract )
        % Move in the direction of the current tangent, but only with
        % quarter of the steplength
        sl = sl/4;
        v = v_old + sl*dv;
        alp = alp_old + sl*da;

        % Have we reached the mininum steplength approved?
        if ( sl < stepmin )
            disp('FAILED: sl < stepmin')
            break
        end

        % Go back to the start of the while loop, haven taken a smaller
        % tangent step:
        continue
    end

    % If we're experiencing good Newton convergence, we try to get the
    % steplength closer to the maximum steplength allowed:
    if ( newtonIter <= 3 )
        sl = min(2*sl, stepmax);
    end

    % We've found a new point, update counter:
    counter = counter + 1;

    % If we've been successful to get here, update old variable values:
    dv_old = dv;
    da_old = da;
    v_old = v;
    alp_old = alp;

    % And store the variables for output:
    aalp(end+1) = alp;
    vv(:,end+1) = v;
    mvec(end+1) = measure(v, wk);

    if ( numel(aalp)  > 2 )
        for kk = 1:numel(nlist)
            tt = linspace(0,20,1001)';
            
            if ( (aalp(end) == nlist(kk) ) || ...
                    (1/log2(aalp(end)) > nlist(kk) && 1/log2(aalp(end-1)) < nlist(kk) ) || ...
                 ( 1/log2(aalp(end-1)) > nlist(kk)  && 1/log2(aalp(end)) < nlist(kk) ) )
                figure(kk+1)
                % Solve at precisely half-integer values of n
                u = newtonSolve(2.^(1./nlist(kk)), v);
                % Interpolate for smooth plots:
                uu = exp(-bet*tt/2).*bary(tt, u.*exp(bet*t/2), t, vk);
                % Plot
                plot(tt, 1+uu, plotdata{:}); hold on
                title(['$n = ' num2str(nlist(kk),'%2.2f') '$'], 'interp', 'latex') 
                xlabel('$t$', 'interp', 'latex');
                ylabel('$u$', 'interp', 'latex');
                set(gca, 'fontsize', fs);
            end
        end
    end
    % alignfigs

    % Is STOPFUN telling us to stop?
    if ( stopfun(v, alp) )
        break
    end

end

    function [dv, dalp] = tangentBVP(v, alp, v_old, a_old)

        % Barycentyric matrices:
        tau = alp*t;
        P = exp(-bet*tau/2).*barymat(tau, t, vk).*exp(bet*t'/2);

        % Linearisation:
        Ju = D + I - 2*P - 2*P*diag(v);
        % Ja = -2*diag(t)*P*D*v - diag(t)*P*D*(v.^2);
        Ja = diag(2/alp.*t.*(-1 - Pv))*D*Pv;       % dF/da
        Jm = [v_old'.*wk, a_old];
        J = [Ju, Ja ; Jm];
        
        % Right hand side
        rhs = [0*v; 1];

        % Solve for the tangent:
        sol = J\rhs;

        % Extract the function and scalar:
        dv = sol(1:end-1);
        dalp = sol(end);

        % Normalize the tangent returned:
        scale = sqrt((dv'.*wk)*dv + dalp^2);
        dv = dv/scale;
        dalp = dalp/scale;

    end


    function [v, alp, newtonIter, retract] = newtonBVP(v, alp, dv, da)
        
        % Barycentyric matrices:
        tau = alp*t;
        P = exp(-bet*tau/2).*barymat(tau, t, vk).*exp(bet*t'/2);

        % Functional condition for tangent step:
        Jm = [dv'.*wk, da];

        % Intialise values:
        normv = norm(v, 2);
        accept = false;
        newtonIter = 0;
        retract = false;

        % Newton iterations:
        while ( ~accept )
            Pv = P*v;
            F = D*v + v - 2*Pv - (Pv).^2;              % Residual 
            Jv = D + I + 2*P*diag(-1-v);               % dF/dv 
            Ja = diag(2/alp.*t.*(-1 - Pv))*D*Pv;       % dF/da
            J = [Jv, Ja ; Jm];

            % RHS is the residual of the differential equation, we then add
            % a 0 at the bottom for the functional condition:
            rhs = [-F; 0];

            % Solve the linear system to obtain a Newton correction:
            sol = J\rhs;

            % Extract the function and scalar parts:
            dv = sol(1:end-1);
            da = sol(end);

            % Take a Newton step. Since we should be close to a solution on
            % the curve, we take a full Newton step.
            v = v + dv;
            alp = alp + da;

            % Increase the iteration counter.
            newtonIter = newtonIter + 1;

            if ( norm(dv, 2)/normv < 5e-4 )
                % Hoorayh. We've converged, so accept the tangent step!
                accept = 1;
            elseif ( newtonIter >= 5 )
                % We wanted to take too many iterations, so tell the
                % followpath algorithm to retract and take a smaller
                % tangent step.
                retract = 1;
                return
            end
        end
    end

    function v = newtonSolve(alp, v)
        tauk = alp*t;
        P_ = exp(-bet*tauk/2).*barymat(tauk, t, vk).*exp(bet*t'/2);
        
        % Newton iteration
        warning('off', 'MATLAB:nearlySingularMatrix');
        for l = 1:10
            Pv_ = P_*v;
            F = D*v + v - 2*Pv_ - (Pv_).^2;
            J = D + I - 2*P_ - 2*P_*diag(v);
        
            dv_ = J\F;
            v = v - dv_;
        
            if ( norm(dv_, inf) < 1e-10 ), break, end
        end
    end


end

function P = barymat(tau,t,w)
  P = w'./(tau-t');
  P = P./sum(P, 2);
  P(isnan(P)) =  1;
end


function plotDecayingSoln(n)

N = 700;
M = 6;

% Laguerre discretisation:
bet = 2;
[t, D, vk] = mylagdif(N, bet);
I = eye(N);
% Remove points that will not contribute:
idx = (1:M*sqrt(N)+1)';
t = t(idx);
vk = vk(idx);
I = I(idx, idx);
D = D(idx, idx);

alp = 2^(1./(n));

% u = exp(-t);
u = 1-exp(t)./(exp(t) + exp(1./((alp-1)*t)));

% Barycentyric matrices:
tau = alp*t;
P = exp(-alp*bet*t/2).*barymat(tau, t, vk).*(exp(bet*t'/2));
E = barymat(0, t, vk).*exp(bet*t'/2);
% Newton
for k = 1:30
    F = D*u + u -  P*u.^2;
    J = D + I - 2*P*diag(u);
    J(1,:) = E; F(1) = 0;
    du = J\F;
    u = u - du;
    if ( norm(du) < 1e-10 ), break, end
end

plot(t, u, ':k'); hold on;  shg
plot(t, 1+0*t, '-k'); 
xlim([0 20]), 
grid on

end


function [x, D, v, w] = mylagdif(N, b)
[x, w, v] = lagpts(N-1, 'glr');
v = v./x; v = [-sum(v) ; v]; x = [0; x]; w = [0, w]; % Include x = 0.
v = v./norm(v,inf);
ii = 1:N+1:N^2;   % Diagonal entries
c = exp(-x/2)./v; C = c./c';
X = 1./(x-x');    X(ii) = -.5;
D = X.*C;         D(ii) = sum(X, 2); 
if ( nargin > 1 )
    x = x/b; D = b*D; w = w/b; % Scale nodes by the factor b.
end
end




