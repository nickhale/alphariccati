% Code to generate Table 1 of Approximate solutions to a nonlinear
% functional differential equation
%
% Nick Hale, Stellenbosch University, Jan 2024

C = sym(zeros(6, 1));
L2 = C;
digits(100);

for n = 1:6
    
    alpn = 2^(1/sym(n, 'd'));

    % Compute the terms in the series (50 is enough):
    k = (0:49)';
    terms = 2./(1-alpn.^k); terms(1) = 1;
    ck = cumprod(terms);

    % Keep only those terms larger than eps (for speed)
    idx = find(abs(ck)<eps, 1, 'first');
    ck = ck(1:idx); k = k(1:idx);

    % Compute the (n-1)th moment of E:
    mu = factorial(n-1)*sum(ck./alpn.^(((n-1)+1)*k));

    % Compute the first n moments of E^2
    num = ck*ck';
    den = alpn.^k+alpn.^k';
    A = num(:)./den(:).^((0:n-1)+1);
    nu =  factorial(0:n-1).*sum(A);

    % Cpmpute the wj
    wj = sym(ones(1,n), 'd');
    for j = n-1:-1:1
        wj(j) = alpn^(j+1)*j/(alpn^j-2)*wj(j+1);
    end

    % Compute the Cn:
    C(n,1) = (2*n/alpn) * mu/sum(wj.*nu);
    L2(n,1) = abs(C(n,1))*sqrt(nu(1));

end

% Display table:
n = (1:n)';
C = double(C);
L2 = double(L2);
disp(table(n, C, L2))

save Cn C
