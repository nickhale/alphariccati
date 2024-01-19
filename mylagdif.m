function [x, D, w, v] = mylagdif(N, b)
% MYLAGDIF - Compute Gauss-Laguerre points and differentiation matrix
%   X = MYLAGDIF(N) returns the N-vector [0 ; XL], where GL are the N-1
%   Gauss-Laguerre points in (0, inf).
%
%   [X, D] = MYLAGDIF(N) returns the corresponding differentiation matrix.
%
%   [X, D, W, V] = MYLAGDIF(N) returns also the corresponding quadrature
%   weights, W, and barycentric weights, V. (Note W are not precisely the
%   optimatal quadrature weights for these nodes; they are [0, WL], where
%   WL are the Gauss-Laguerre quadrature weights corresponding to XL.)
%
%   [X, D, W, V] = MYLAGDIF(N, B) incorporates the scaling X->X/B in all
%   the outputs.
%
%   See also lagdif.m from DMSUITE (J.A.C. Weideman & S.C. Reddy 1998)

% Nick Hale, Stellenbosch University, 2024

try 
    % This version uses a not-yet released branch of Chebfun which avoids
    % underflow and hence NaNs in Laguerre differentiation matrix.
    [x, w, ~, ders] = lagpts(N-1, 'glr'); % www.chebfun.org
    v = [1; exp(-x/2)./(x.*ders)]; x = [0; x]; w = [0, w]; % Add pt at 0
    c = [1; ders.*x(2:end)];       % c = exp(-x/2)./v
catch
    [x, w, v] = lagpts(N-1, 'glr');         % www.chebfun.org
    v = v./x; v = [-sum(v) ; v]; x = [0; x]; w = [0, w];   % Add pt at 0
    c = exp(-x/2)./v;
end
C = c./c'; 

ii = 1:N+1:N^2;             % Diagonal entries
X = 1./(x-x');              X(ii) = -.5;
D = X.*C;                   D(ii) = sum(X, 2); 

if ( nargin > 1 ) % Scaling:
    x = x/b; D = b*D; w = w/b; % Scale nodes by the factor b.
end

end