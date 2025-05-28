function u=u_protected(z)

om = exp(2i*pi/3);
K = 4*pi/3;
zS = 1i/sqrt(3);
e1 = om^2;                     % basis of Lambda
e2 = -om;
f1 = (4i*pi/sqrt(3))*om;     % dual basis of Lambda^* (<e_i, f_j> = 2pi 1_{i==j})
f2 = (4i*pi/sqrt(3))*om^2;

N = 32;
Up = sym_potential(N,0,0);
Um = Up.';

M = 600;
alpha = 20;
[~,~,U] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2*Up*Um, 1, 'smallest');
[~,uu] = F2X(U,0,M,e1,e2);

u = interpolate(e1,e2,uu,M,real(z(:))',imag(z(:))');
u = reshape(u, size(z));


end



% FUNCTIONS

% Fourier representation: a vector is a length (2N+1)^2 column vector V, represent 
% v(z) = \sum_{a,b\in{-N,...,N}} V((2N+1)*(b+N)+a+N+1) e^{i<z, a f1 + b f2+ k>} / sqrt(area_of_(e1,e2))

% Position-space representation: v(z) for
% z  a MxM matrix  with values   z(a,b) = e1*(a-1)/M + e2*(b-1)/M

% Derivatives --------------

function Dbar = Dbar(N,k,f1,f2)
    Dbar = 0.5 * (Dx(N,k,f1,f2) + 1i .* Dy(N,k,f1,f2));
end

function Dx = Dx(N,k,f1,f2)
    D0 = spdiags((-N:1:N)', 0, 2*N+1, 2*N+1); E = speye(2*N+1, 2*N+1);
    Dx = kron(E, real(f1) * D0) + kron(real(f2) * D0, E) + real(k)*kron(E,E);
end

function Dy = Dy(N,k,f1,f2)
    D0 = spdiags((-N:1:N)', 0, 2*N+1, 2*N+1); E = speye(2*N+1, 2*N+1);
    Dy = kron(E, imag(f1) * D0) + kron(imag(f2) * D0, E) + imag(k)*kron(E,E);
end

function A = Inv(B) % invert a diagonal matrix
    n = size(B,1);
    A = spdiags(1./diag(B), 0, n, n);
end

% Multiply by potential ------------

function U=fourier_shift(N,n1,n2)
    N = 2*N + 1;
    U = kron(spdiags(ones(N,1), -n2, N, N),  spdiags(ones(N,1), -n1, N, N));
end

% (2K+) Potential with TBG rotation symmetry containing K + n1 f1 + n2 f2
function U=sym_potential(N,n1,n2)
    w = exp(2i*pi/3);
    U = (-4i*pi/3)* (fourier_shift(N,n1-1,n2+1) + w*fourier_shift(N,-n2-1,n1-n2) + w^2*fourier_shift(N,n2-n1, -n1+1));
end

% Fourier representation of constant 1 function
% WARNING- this is for TBG lattice. for general lattice factor should be sqrt(abs(imag(e1' * e2)))
function V=const_fun(N)
    V = kron((-N:N).'==0, (-N:N).'==0) * sqrt(sqrt(3)/2); 
end


function [z,v] = F2X(V,k,M,e1,e2)
    N = round((sqrt(length(V))-1)/2);
    V0 = reshape(V,2*N+1, 2*N+1);
    V = zeros(M);
    V(1:2*N+1, 1:2*N+1)=V0;
    V = circshift(V, [-N,-N]);

    v = ifft2(V) * M^2  / sqrt(abs(imag(e1' * e2)));

    [y2,y1]=meshgrid(0:M-1, 0:M-1);
    z=(e1*y1+e2*y2)/M;

    v = v .* exp(0.5i*(k'*z+k*conj(z)));
end

function fQ=interpolate(e1,e2,f,M,xQ,yQ)
    % given function on fund domain (output of K2X), interpolate to z=(x,y)
    T = [real(e1), real(e2); imag(e1), imag(e2)];
    v = T^(-1) * [xQ; yQ];
    v = v - floor(v); v = ceil(v*M); v(v==0) = 1;
    ind = v(1,:) + M*(v(2,:)-1);
    fQ = f(ind);
end

