om = exp(2i*pi/3);
K = 4*pi/3;
zS = 1i/sqrt(3);

e1 = om^2;                 % basis of Lambda
e2 = -om;

f1 = (4i*pi/sqrt(3))*om;   % dual basis of Lambda^*
f2 = (4i*pi/sqrt(3))*om^2;

N = 32;
M = 300;

W = @(n) 1i*K^2 ...
    *(fourier_shift(N,-n,-n)+fourier_shift(N,n,n) ...
    + om*fourier_shift(N,n,0)+om*fourier_shift(N,-n,0) ...
    + om'*fourier_shift(N,0,n)+om'*fourier_shift(N,0,-n));
W = 0.5*W(1)+0.5*W(2);

% Protected state
alpha = 4;
P = [2*Dbar(N,0,f1,f2), speye(size(W)); alpha^2*W^2, 2*Dbar(N,0,f1,f2)];
[~,s,U] = svds(P, 1, 'smallest');
disp(s)

U2 = P \ U;
U = U(1:length(U)/2);
U2 = U2(1:length(U2)/2);
R = ROT(N,0);
U2= (R^0 + om*R^1 + om^2*R^2)* U2/3;

[z,u] = F2X(U,M,e1,e2);
u = u / exp(1i*angle(u(1,1)));
[~,u2] = F2X(U2,M,e1,e2);
u2 = u2 / exp(1i*angle(u2(1,1)));

disp(norm(imag(u)) / norm(u))
figure
contourf(real(z), imag(z), real(u), 32)
colorbar
axis equal


% Phase function \varphi
ww = sparse(2*N+1, 2*N+1); ww(N+1,N+1) = 1;
Phi = Inv(2*Dbar(N,0,f1,f2) + kron(ww,ww))*W*const_fun(N);
[~,phi] = F2X(Phi,M,e1,e2);

disp(norm(real(phi))/ norm(phi))
phi = imag(phi);
figure
contourf(real(z), imag(z), phi, 32)
colorbar
axis equal

bb = diag(u); bb = bb(round(M*4/9):round(M*5/9));
gg = diag(cos(alpha*phi)); gg = gg(round(M*4/9):round(M*5/9));
C = gg'*bb / (gg'*gg);

figure
title("Comparison of $\rm{Re}\,u$ and $C(\alpha)\cos " + ...
    "(\alpha\varphi)$ on the segment $[-\sqrt{3}i, 0]$",...
    "Interpreter", "Latex");
hold on
plot(imag(diag(z)), real(diag(u)))
plot(imag(diag(z)), C*cos(alpha*diag(phi)))

%% Thingy for u2
bb = diag(u2); bb = bb(round(M*4/9):round(M*5/9));
gg = diag(cos(alpha*phi)); gg = gg(round(M*4/9):round(M*5/9));
ff = diag(sin(alpha*phi)); ff = ff(round(M*4/9):round(M*5/9));

C1 = gg'*bb / (gg'*gg);
C2 = ff'*bb / (ff'*ff);
disp(norm(bb))
disp(norm(bb - C1*gg - C2*ff))
figure
contourf(real(z), imag(z), abs(u2-C1*cos(alpha*phi)-C2*sin(alpha*phi)))
axis equal
colorbar

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

% Transform from fourier basis to position basis
function [z,u] = F2X(U,M,e1,e2)
    N = round((sqrt(length(U))-1)/2);
    V0 = reshape(U,2*N+1, 2*N+1);
    U = zeros(M);
    U(1:2*N+1, 1:2*N+1)=V0;
    U = circshift(U, [-N,-N]);

    u = ifft2(U) * M^2  / sqrt(abs(imag(e1' * e2)));

    [y2,y1]=meshgrid(0:M-1, 0:M-1);
    z=(e1*y1+e2*y2)/M;
end

% Fourier representation of constant 1 function for the lattice \Lambda
function V=const_fun(N)
    V = kron((-N:N).'==0, (-N:N).'==0) * sqrt(sqrt(3)/2); 
end


% For k in {0,K,-K}, R[u](x) = u(wx) acts on L^2_k(C;C^1)
function R = ROT(N,s)         % k = sK
    F1 = @(n1,n2) s+n2-n1;     % multiplication by w^-1 on k+Lambda^*
    F2 = @(n1,n2) -n1;
    R = PERMK(N,F1,F2);
end


% General operator permuting fourier basis
% (FV)(F1(n1,n2),F2(n1,n2)) = V(n1,n2)
function F = PERMK(N,F1,F2)
    indx =  @(m1,m2) (2*N+1)*(m2+N) + m1+N + 1;
    [n2,n1]=meshgrid(-N:N,-N:N); n1=n1(:); n2=n2(:);

    fm1 = F1(n1,n2);
    fm2 = F2(n1,n2);

    mask = fm1 <= N & fm1 >= -N & fm2 <= N & fm2 >= -N;
    F = sparse(indx(fm1(mask),fm2(mask)), indx(n1(mask),n2(mask)), ones(length(find(mask)),1), (2*N+1)^2, (2*N+1)^2);
end