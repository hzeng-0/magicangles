% Hexagonal Lattice %

w = exp(2i*pi/3);
K = 4*pi/3;
zS = 1i/sqrt(3);

e1 = w^2;                     % basis of Lambda
e2 = -w;

f1 = (4i*pi/sqrt(3))*w;     % dual basis of Lambda^*
f2 = (4i*pi/sqrt(3))*w^2;



% Hamiltonian %

N = 32;
k = 0;

% Dbar

Db_1 = Dbar(N, k-K, f1, f2);
Db_2 = Dbar(N, k+K, f1, f2);

Dx_1 = Dx(N,k-K,f1,f2);
Dy_1 = Dy(N,k-K,f1,f2);

Dx_2 = Dx(N,k+K,f1,f2);
Dy_2 = Dy(N,k+K,f1,f2);

% Potential

Up = sym_potential(N,0,0);
Um = Up.';


D = [2*Db_1, 0*Up; 0*Up, 2*Db_2];
U = [0*Up, Up; Um, 0*Up];



%% Eigenvalues of H: fixed alpha

H = [0*D, (D+5*U)'; (D+5*U), 0*D];

E = eigs(H, 400, 'smallestabs');
E = E(E > 0);


figure(1); hold on;
plot(E, (1:200)');
plot(E, sqrt(3)/(4*pi) * E.^2);

figure(2); hold on;
plot(E, (sqrt(3)/(4*pi) * E.^2  - (1:200)' ) ./ E);
yline(0);


%% Eigenvalues of H: varying alpha

alphas = 0:0.1:15;

EE = zeros(length(alphas),40);
for id=1:length(alphas)
    disp(alphas(id))
    EE(id,:) = svds(D+alphas(id)*U, 40, 'smallest');
end

figure(10); hold on;
plot(alphas,EE,'Color','b');
xlabel("\alpha")
ylabel("Positive eigenvalues of $H_0(\alpha)$",'Interpreter','latex')


%% Potential U

inp = @(x,y) (x'*y+x*y')/2;
U_fun = @(x) exp(i*inp(x,K)) + w*exp(i*inp(x,w*K)) +  w^2*exp(i*inp(x,w^2*K));

V_fun = @(x) U_fun(x) * U_fun(-x);

V = @(x) abs(exp(i*x*K) - exp(-i*x*K/2));
xx = 0:0.001:1;
VV = arrayfun(V, xx);
disp(sum(VV) / length(VV));


%% Calculate magic angles with T_k

Ak = Inv(2*Db_1) * Up * Inv(2*Db_2) * Um;

Alphas = 1./sqrt(eigs(Ak, 500));

figure; 
hold on;
scattermult([real(Alphas), imag(Alphas)], 5);



% Gap of real magic angles
RealAlphas = Alphas(abs(imag(Alphas)) < 0.01);

fprintf('Gap of real magic angles is about %d\n', RealAlphas(6)-RealAlphas(5));




%% Eigenfunctions


[~,s,V] = svds(D + 20*U + 0*speye(size(D)), 1, 'smallest');
V1 = V(1:length(V)/2); V2 = V(length(V)/2 + 1:end);
[z,v1] = K2X3(V1,K,200,e1,e2);

figure; 
hold on;
contourf(real(z), imag(z), log(abs((v1))), 30);
hex(zS);
axis equal; 
colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]);
hold off;




%% Bracket

figure; hold on;

[boundary,unsure,areas,count,zz,bra]=bracket_area([0 0 1],600,3,4);
disp([areas, count]);

title('$|\{q, \bar q\}|$','Interpreter','latex');
surf(real(zz),imag(zz),abs(bra),'EdgeColor','none');

scatter3(real(boundary), imag(boundary), 1+ 0*boundary,  3,'red','filled');
scatter3(real(unsure), imag(unsure), 1+ 0*unsure, 4,'red','filled');

axis equal;


%% FUNCTIONS --------------------------------------------------------------

% Derivatives

function Dbar = Dbar(N,k,f1,f2)

    D0 = spdiags((-N:1:N)', 0, 2*N+1, 2*N+1); E = speye(2*N+1, 2*N+1);

    Dbar = 0.5 * (kron(D0 * f1, E) + kron(E, D0 * f2) + k*kron(E,E));
 
end

function A = Inv(B)

    n = size(B,1);

    A = spdiags(1./diag(B), 0, n, n);

end

% Multiply by potential


function U=sym_potential(N,n1,n2) % potential with symmetries for TBG

    U = (-4i*pi/3)* (fourier_shift(N,n1-1,n2+1) + w*fourier_shift(N,-n2-1,n1-n2) + w^2*fourier_shift(N,n2-n1, -n1+1));
end

function U=fourier_shift(N,n1,n2)

    N = 2*N + 1;

    U = kron(spdiags(ones(N,1), -n1, N, N),  spdiags(ones(N,1), -n2, N, N));

end



% Get position-space vector v corresponding to k-space vector V:

% v(z) = sum V(n1,n2) e^(i<n1 f1 + n2 f2 + k, z>)

function [z,v] = K2X(V,k,M,e1,e2)

    N = round((sqrt(length(V))-1)/2);

    V0 = reshape(V,2*N+1, 2*N+1); V0 = V0.';
    V = zeros(M);
    V(1:2*N+1, 1:2*N+1)=V0;
    V = circshift(V, [-N,-N]);

    v = fft2(V);

    [y1,y2]=meshgrid(0:M-1, 0:M-1);
    z=(e1*y1+e2*y2)/M;

    v = v .* exp(0.5i*(k'*z+k*conj(z)));

end

% For graphing (makes 4 copies)

function [z,v] = K2X3(V,k,M,e1,e2)

    [z0,v0]=K2X(V,k,M,e1,e2);

    w = exp(2*pi*1i/3);
    z = [z0-e1-e2, z0-e2; z0-e1, z0];
    
    p2 = exp(-0.5i*(k'*e2+k*e2')); p1 = exp(-0.5i*(k'*e1+k*e1'));
    v = [v0 *p1 *p2, v0 *p2; v0 *p1, v0];

end



% Operator permuting fourier basis:

% (FV)(F1(m1,m2),F2(m1,m2)) = V(m1,m2)  

function F = PermK(N,F1,F2)

    indx =  @(m1,m2) (2*N+1)*(m1+N) + m2+N + 1;

    [m1,m2]=meshgrid(-N:N,-N:N); 
    m1=m1(:); m2=m2(:);

    fm1 = F1(m1,m2);
    fm2 = F2(m1,m2);

    mask = fm1 <= N & fm1 >= -N & fm2 <= N & fm2 >= -N;

    F = sparse(indx(fm1(mask),fm2(mask)), indx(m1(mask),m2(mask)), ones(length(find(mask)),1), (2*N+1)^2, (2*N+1)^2);

end

