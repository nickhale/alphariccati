function [x, D, w, v] = mylagdif(N, b)
% MYLAGDIF - Compute Gauss-Laguerre points and differentiation matrix
%   X = MYLAGDIF(N) returns the N-vector [0 ; XL], where GL are the N-1
%   Gauss-Laguerre points.
%
%   [X, D] = MYLAGDIF(N) returns the corresponding differentiation matrix.
%
%   [X, D, WL, VL] = MYLAGDIF(N) returns also the corresponding barycentric
%   weights, VL, and quadrature weights, W. (Note W are not precisely the
%   optimatal quadrature weights for these nodes; they are [0, WL], where
%   WL are the Gauss-Laguerre quadrature weights corresponding to XL.)
%
%   [X, D, W, V] = MYLAGDIF(N, B) incorporates the scaling X->X/B in all
%   the outputs.

% Nick Hale, Stellenbosch University, 2024

[x, w, ~, ders] = lagpts(N-1, 'glr'); % www.chebfun.org
v = [1; exp(-x/2)./(x.*ders)]; x = [0; x]; w = [0, w]; % Add pt at 0
c = [1; ders.*x(2:end)];    C = c./c';          % c = exp(-x/2)./v
ii = 1:N+1:N^2;             % Diagonal entries
X = 1./(x-x');              X(ii) = -.5;
D = X.*C;                   D(ii) = sum(X, 2); 

if ( nargin > 1 ) % Scaling:
    x = x/b; D = b*D; w = w/b; % Scale nodes by the factor b.
end

end