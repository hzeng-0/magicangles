function result=phi(z)

om = exp(2i*pi/3);
K = 4*pi/3;
zS = 1i/sqrt(3);
e1 = om^2;                     % basis of Lambda
e2 = -om;
f1 = (4i*pi/sqrt(3))*om;     % dual basis of Lambda^* (<e_i, f_j> = 2pi 1_{i==j})
f2 = (4i*pi/sqrt(3))*om^2;

N = 32;
W = 1i*K^2 ...
    *(fourier_shift(N,-1,-1)+fourier_shift(N,1,1) ...
    + om*fourier_shift(N,1,0)+om*fourier_shift(N,-1,0) ...
    + om'*fourier_shift(N,0,1)+om'*fourier_shift(N,0,-1));

ww = sparse(2*N+1, 2*N+1); ww(N+1,N+1) = 1;
Phi = Inv(2*Dbar(N,0,f1,f2)+kron(ww,ww)) * W * const_fun(N);
phi = F2Vect(Phi,z,f1,f2,N);
result = imag(phi);

end



% FUNCTIONS

% Fourier representation: a vector is a length (2N+1)^2 column vector V, represent 
% v(z) = \sum_{a,b\in{-N,...,N}} V((2N+1)*(b+N)+a+N+1) e^{i<z, a f1 + b f2+ k>} / sqrt(area_of_(e1,e2))

% Position-space representation: v(z) for
% z  a MxM matrix  with values   z(a,b) = e1*(a-1)/M + e2*(b-1)/M

function Dbar = Dbar(N,k,f1,f2)
    D0 = spdiags((-N:1:N)', 0, 2*N+1, 2*N+1); E = speye(2*N+1, 2*N+1);
    Dbar = 0.5 * (kron(E, f1 * D0) + kron(f2 * D0, E) + k*kron(E,E));
end

function A = Inv(B) % invert a diagonal matrix
    n = size(B,1);
    A = spdiags(1./diag(B), 0, n, n);
end

function U=fourier_shift(N,n1,n2)
    N = 2*N + 1;
    U = kron(spdiags(ones(N,1), -n2, N, N),  spdiags(ones(N,1), -n1, N, N));
end

function V=const_fun(N)
    V = kron((-N:N).'==0, (-N:N).'==0) * sqrt(sqrt(3)/2); 
end

function v = F2Vect(V,z,f1,f2,Nmax)
    N = round((sqrt(length(V))-1)/2);
    V = reshape(V,2*N+1, 2*N+1);
    v = 0 * z;

    for ii = -Nmax:Nmax
        for jj = -Nmax:Nmax
            k = ii * f1 + jj * f2;
            v = v + V(ii+N+1, jj+N+1) * exp(0.5i*(k'*z + k*conj(z)));
        end

    end
    v = v / sqrt(sqrt(3)/2);
end


