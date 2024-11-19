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

n1 = 0; n2 = 0;

Up = (-4i*pi/3)* (fourier_shift(N,n1-1,n2+1) + w*fourier_shift(N,-n2-1,n1-n2) + w^2*fourier_shift(N,n2-n1, -n1+1));
% 
% n1 = 1; n2 = -1;
% Up = Up + 0.1 * (-4i*pi/3)* (fourier_shift(N,n1-1,n2+1) + w*fourier_shift(N,-n2-1,n1-n2) + w^2*fourier_shift(N,n2-n1, -n1+1));
% 
Um = Up.';


D = [2*Db_1, 0*Up; 0*Up, 2*Db_2];
U = [0*Up, Up; Um, 0*Up];



% Scalar model
D_scalar = @(k,alpha) 4*Dbar(N,k,f1,f2)^2 - alpha^2*Up*Um;
D_scalar2 = @(k,alpha) 4*Dbar(N,k+K,f1,f2)*Dbar(N,k-K,f1,f2) - alpha^2*Up*Um;


%%

% Singular values %

% Weyl law as lambda change, fixed alpha

H = [0*D, (D+5*U)'; (D+5*U), 0*D];

E = eigs(H, 400, 'smallestabs');
E = E(E > 0);


figure(1); hold on;
plot(E, (1:200)');
plot(E, sqrt(3)/(4*pi) * E.^2);

figure(2); hold on;
plot(E, (sqrt(3)/(4*pi) * E.^2  - (1:200)' ) ./ E);
yline(0);



%%
% Varying alpha

alphas = 0:0.1:15;

EE = zeros(length(alphas),40);
for id=1:length(alphas)
    disp(alphas(id))
    EEE = eigs([0*D, (D+alphas(id)*U)'; (D+alphas(id)*U), 0*D], 90, 'smallestabs');
    EEE=EEE(EEE>0);
    EE(id,:)=EEE(1:40);
end

figure(10); hold on;
plot(alphas,EE);
xline(5); xline(6);

%%
alphas = 0:0.1:15;

EE = zeros(length(alphas),40);
for id=1:length(alphas)
    disp(alphas(id))
    EE(id,:) = svds(D+alphas(id)*U, 40, 'smallest');
end

%%
figure(10); hold on;
plot(alphas,EE,'Color','b');
xlabel("\alpha")
ylabel("Positive eigenvalues of $H_0(\alpha)$",'Interpreter','latex')
%%xline(5); xline(6);


%%
[V,DD] = eigs([0*D, (D+10*U)'; (D+10*U), 0*D], 60, 'smallestabs');
% V=V(:,59); % eigenfunction for 2.13 eigenvalue (3.6-19th)
V = V(:,5);
V=V(1:length(V)/2);


V1 = V(1:length(V)/2); V2 = V(length(V)/2 + 1:end);

[z,v1] = K2X3(V2,K,200,e1,e2);

figure(8); hold on;
contourf(real(z), imag(z), log(abs(v1)), 10);
hex(zS);
axis equal;

%%



[V,DD] = svds(D_scalar(0,10), 3, 'smallest');
% [V,DD] = svds(D + 10*U + K*speye(size(D)), 3, 'smallest');
V = V(:,3);
% V = V(1:length(V)/2);
V = V / V(round(length(V) / 2));
disp(DD);

[z,v1] = K2X3(V,0,400,e1,e2);
figure; hold on;
contourf(real(z), imag(z), ((imag(v1))), 10);
hex(zS);
axis equal; colorbar;




%%

a = 0:0.01:5;
s = a;

for id = 1:length(a)
    disp(id)
s(id) = svds(D_scalar(0.3,a(id)), 1, 'smallest');


end

figure;
plot(a, s);

%%

s2 = a;

for id = 1:length(a)
    disp(id)
s2(id) = svds(D + a(id)*U, 1, 'smallest');


end

figure;
plot(a, s2);

%%

inp = @(x,y) (x'*y+x*y')/2;
U_fun = @(x) exp(i*inp(x,K)) + w*exp(i*inp(x,w*K)) +  w^2*exp(i*inp(x,w^2*K));

V_fun = @(x) U_fun(x) * U_fun(-x);

V = @(x) abs(exp(i*x*K) - exp(-i*x*K/2));
xx = 0:0.001:1;
VV = arrayfun(V, xx);
disp(sum(VV) / length(VV));


%%


% Magic angles as spectrum of compact operator

Ak = Inv(2*Db_1) * Up * Inv(2*Db_2) * Um;

Alphas = 1./sqrt(eigs(Ak, 500));

figure(3); hold on;
scattermult([real(Alphas), imag(Alphas)], 5);



% Gap of real magic angles
RealAlphas = Alphas(abs(imag(Alphas)) < 0.01);

fprintf('Gap of real magic angles is about %d\n', RealAlphas(6)-RealAlphas(5));




%%

% Eigenfunctions

% Protected state

[~,s,V] = svds(D + 20*U + 0*speye(size(D)), 1, 'smallest');

V1 = V(1:length(V)/2); V2 = V(length(V)/2 + 1:end);

[z,v1] = K2X3(V1,K,200,e1,e2);

figure(7); hold on;
contourf(real(z), imag(z), log(abs((v1))), 30);
hex(zS);
axis equal; colorbar; 
xlim([-0.63,0.63]), ylim([-0.63,0.63]);
hold off;



%% D log v 11/13/2024

alpha = 5;

[~,s,V] = svds(D + alpha*U + 0*speye(size(D)), 5, 'smallest');
V = V(:,5);

V = V ./ V(round(length(V)/4));

V1 = V(1:length(V)/2); V2 = V(length(V)/2 + 1:end);



[z,v1] = K2X3(V1,0,200,e1,e2);
[z,v2] = K2X3(V2,-K,200,e1,e2);

[z,dxv1] = K2X3(Dx_1*V1 ./alpha,0,200,e1,e2);
[z,dyv1] = K2X3(Dy_1*V1 ./alpha,0,200,e1,e2);


[z,dxv2] = K2X3(Dx_2*V2 ./alpha,-K,200,e1,e2);
[z,dyv2] = K2X3(Dy_2*V2 ./alpha,-K,200,e1,e2);

levels=linspace(-200,8,50);

% figure(8); hold on;
% contourf(real(z), imag(z), log(abs(dxv1)), linspace(-5,3,20));
% hex(zS);
% axis equal; colorbar;  
% xlim([-0.63,0.63]), ylim([-0.63,0.63]);
% hold off;

% figure(8); hold on;
% contourf(real(z), imag(z), (imag(dyv1 ./ v1)), levels);
% hex(zS);
% axis equal; colorbar; 
% xlim([-0.63,0.63]), ylim([-0.63,0.63]);
% hold off;

figure(9); hold on;
surf(real(z),imag(z),angle(v1),'EdgeColor','none');
colormap(gca,wheelmap)
hex(zS,10);
colorbar;
axis equal;
xlim([-0.63,0.63]), ylim([-0.63,0.63]);
clim([-pi pi])

figure(10); hold on;
surf(real(z),imag(z),angle(dxv1),'EdgeColor','none');
colormap(gca,wheelmap)
hex(zS,10);
colorbar;
axis equal;
xlim([-0.63,0.63]), ylim([-0.63,0.63]);
clim([-pi pi])

figure(11); hold on;
surf(real(z),imag(z),angle(dyv1),'EdgeColor','none');
colormap(gca,wheelmap)
hex(zS,10);
colorbar;
axis equal;
xlim([-0.63,0.63]), ylim([-0.63,0.63]);
clim([-pi pi])



%% 

% Bracket

figure(5); hold on;

[boundary,unsure,areas,count,zz,bra]=bracket_area([0 0 1],600,3,4);
disp([areas, count]);

title('$|\{q, \bar q\}|$','Interpreter','latex');
surf(real(zz),imag(zz),abs(bra),'EdgeColor','none');

scatter3(real(boundary), imag(boundary), 1+ 0*boundary,  3,'red','filled');
scatter3(real(unsure), imag(unsure), 1+ 0*unsure, 4,'red','filled');

axis equal;


% Integrate bracket??


disp( sum(bra, 'all') / length(zz) )
disp( sum(abs(bra), 'all') / length(zz) )
disp( sum(abs(bra).^2, 'all') / length(zz) )
disp( sum(sqrt(abs(bra)), 'all') / length(zz) )


%%

% Scalar model

inp = @(x,y) (x'*y + x*y')/2;
U = @(z) exp(i*inp(z,K))+ w*exp(i*inp(z,w*K)) + w^2 * exp(i*inp(z,w^2 * K));

VV = @(z) U(z) * U(-z);

z = 0.5 + 1i*(-0.5/sqrt(3):0.0001:0.5/sqrt(3));

val = sum(arrayfun(@(z) sqrt(VV(z)),z)) / (1/sqrt(3))


%%

% FUNCTIONS --------------------------------------------------------------

% Derivatives


function Dx = Dx(N,k,f1,f2)

    D0 = spdiags((-N:1:N)', 0, 2*N+1, 2*N+1); E = speye(2*N+1, 2*N+1);

    Dx = 0.5 * (kron(D0 * real(f1), E) + kron(E, D0 * real(f2)) + k*kron(E,E));
 
end

function Dy = Dy(N,k,f1,f2)

    D0 = spdiags((-N:1:N)', 0, 2*N+1, 2*N+1); E = speye(2*N+1, 2*N+1);

    Dy = 0.5 * (kron(D0 * imag(f1), E) + kron(E, D0 * imag(f2)) + k*kron(E,E));
 
end


function Dbar = Dbar(N,k,f1,f2)

    D0 = spdiags((-N:1:N)', 0, 2*N+1, 2*N+1); E = speye(2*N+1, 2*N+1);

    Dbar = 0.5 * (kron(D0 * f1, E) + kron(E, D0 * f2) + k*kron(E,E));
 
end

function A = Inv(B)

    n = size(B,1);

    A = spdiags(1./diag(B), 0, n, n);

end

% Multiply by potential

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

% For graphing (makes a copy)

function [z,v] = K2X3(V,k,M,e1,e2)

    [z0,v0]=K2X(V,k,M,e1,e2);

    w = exp(2*pi*1i/3);
    z = [z0-e1-e2, z0-e2; z0-e1, z0];
    
    p2 = exp(-0.5i*(k'*e2+k*e2')); p1 = exp(-0.5i*(k'*e1+k*e1'));
    v = [v0 *p1 *p2, v0 *p2; v0 *p1, v0];

end



% Operator permuting fourier basis:

% (FV)(F1(m1,m2),F2(m1,m2)) = V(m1,m2)  

function F = Perm(N,F1,F2)

    indx =  @(m1,m2) (2*N+1)*(m1+N) + m2+N + 1;

    [m1,m2]=meshgrid(-N:N,-N:N); 
    m1=m1(:); m2=m2(:);

    fm1 = F1(m1,m2);
    fm2 = F2(m1,m2);

    mask = fm1 <= N & fm1 >= -N & fm2 <= N & fm2 >= -N;

    F = sparse(indx(fm1(mask),fm2(mask)), indx(m1(mask),m2(mask)), ones(length(find(mask)),1), (2*N+1)^2, (2*N+1)^2);

end

