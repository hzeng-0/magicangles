om = exp(2i*pi/3);
K = 4*pi/3;
zS = 1i/sqrt(3);

e1 = om^2;                 % basis of Lambda
e2 = -om;

f1 = (4i*pi/sqrt(3))*om;   % dual basis of Lambda^*
f2 = (4i*pi/sqrt(3))*om^2;

N = 64;
M = 300;

W = 1i*K^2 ...
    *(fourier_shift(N,-1,-1)+fourier_shift(N,1,1) ...
    + om*fourier_shift(N,1,0)+om*fourier_shift(N,-1,0) ...
    + om'*fourier_shift(N,0,1)+om'*fourier_shift(N,0,-1));

% Scalar magic angles
Ak = Inv(2*Dbar(N,0.1,f1,f2))^2 * W^2;
Alphas = 1./sqrt(eigs(Ak, 400));

realAlphas = Alphas(abs(imag(Alphas)) < 0.01);
disp(realAlphas)
disp(diff(realAlphas(1:2:end)))

figure
plot(real(Alphas), imag(Alphas), '.')

% Dbar = 1/(2i)(\partial_x + i\partial_y)
function Dbar = Dbar(N,k,f1,f2)
    D0 = spdiags((-N:1:N)', 0, 2*N+1, 2*N+1);
    E = speye(2*N+1, 2*N+1);
    Dbar = 0.5*( kron(E, f1 * D0) + kron(f2 * D0, E) + k*kron(E,E) );
end

% Invert a diagonal matrix
function A = Inv(B) 
    n = size(B,1);
    A = spdiags(1./diag(B), 0, n, n);
end

% Multiply by e^{i<z, n1 f1 + n2 f2>}
function U=fourier_shift(N,n1,n2)
    N = 2*N + 1;
    U = kron(spdiags(ones(N,1), -n2, N, N),  ...
        spdiags(ones(N,1), -n1, N, N));
end