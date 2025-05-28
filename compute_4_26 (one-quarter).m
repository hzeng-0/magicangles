% Let's look at \bar V^2 potential


om = exp(2i*pi/3);
K = 4*pi/3;
zS = 1i/sqrt(3);
e1 = om^2;                     % basis of Lambda
e2 = -om;
f1 = (4i*pi/sqrt(3))*om;     % dual basis of Lambda^* (<e_i, f_j> = 2pi 1_{i==j})
f2 = (4i*pi/sqrt(3))*om^2;


N = 32;
M = 600;

Up = sym_potential2(N,0,0,1);
Um = Up.';

Up1 = sym_potential2(N,0,0,2);
Um1 = Up1.';

% W = 1i*Up1 * Um1;
W = @(n) 1i*K^2 ...
    *(fourier_shift(N,-n,-n)+fourier_shift(N,n,n) ...
    + om*fourier_shift(N,n,0)+om*fourier_shift(N,-n,0) ...
    + om'*fourier_shift(N,0,n)+om'*fourier_shift(N,0,-n));
W = W(1);


% 
% [z,dzw] = F2X3(1i*Dbar(N,0,f1,f2)'*W * const_fun(N),0,M,e1,e2);
% [z,w] = F2X3(W * const_fun(N),0,M,e1,e2);
% 
% backet = real(conj(w) .* dzw .* w);
% figure
% hold on
% contourf(real(z), imag(z), abs(conj(w) .*dzw .* w), 32);
% % contourf(real(z), imag(z), imag(dzw), 32);
% % 
% axis equal
% xlim([-0.63,0.63]), ylim([-0.63,0.63])
% hex(zS)
% colorbar
% 
% figure
% hold on
% plot(1:length(diag(backet)), diag(backet))


%% Magic angles
Tk = Inv(2 * Dbar(N,0.3, f1, f2))^2 * W^2;
% Tk = Inv(2 * Dbar(N,-K, f1, f2)) * Up1' *  Inv(2 * Dbar(N,0.3- K, f1, f2)) * Um1';
% Tk = -Inv(2*Dbar(N, -K, f1, f2)) * Up'^2 * Inv(2*Dbar(N, +K, f1, f2)) * Um'^2;
Alphas = 1./sqrt(eigs(Tk, 400));
Alphas = [Alphas; -Alphas];

imagAlphas = real(Alphas(abs(imag(Alphas))< 0.1));
imagAlphas = imagAlphas(imagAlphas > 0);
imagAlphas = imagAlphas(1:2:end);

% space = imagAlphas(13) - imagAlphas(12);

points = 0.125 + 0.25 * (0:10);

%%
figure
hold on
scattermult2([real(Alphas), imag(Alphas)], 16, 'blue')
xlim([-0.01 1.6])
ylim([-1.5 1.5])
% plot(real(points), imag(points), 'x');

%%



%% Protected state


%--------


% alpha = 1;
fig=figure;
% vid=VideoWriter('.\results\bababoi.mp4','MPEG-4'); open(vid)
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% for alpha = 0:0.02:2
alphas = 0:0.1:3;
norms = 0*alphas;
for id = 1:length(alphas)
    alpha = alphas(id);
    disp(alpha)
 P = (2*Dbar(N,0,f1,f2))^2 - alpha^2 * W^2;
[~,~,U1] = svds(P, 1, 'smallest');
[z,u1] = F2X3(U1,0,M,e1,e2);
phase = exp(1i*angle(u1(1,1)));
u1 = u1 ./ phase;
norms(id) = norm(imag(u1)) / norm(u1);

end
% P = [2*Dbar(N,0,f1,f2), speye(size(Up)); alpha^2*W^2, 2*Dbar(N,0,f1,f2)];
%--------
 % P = [2*Dbar(N,0,f1,f2), speye(size(Up)); alpha^2*W^2, 2*Dbar(N,0,f1,f2)];

%%
figure
plot(alphas(2:end), log(norms(2:end)))

% U2 = P \ U1;
% U1 = U1(1:length(U1)/2);
% U2 = U2(1:length(U2)/2);

% R = ROT(N,0);
% U2= (R^0 + om*R^1 + om^2*R^2)* U2/3;

%%
U1 = U1 / norm(U1);
% U2 = U2 / norm(U2); % rotation symmetrize
[z,u1] = F2X3(U1,0,M,e1,e2);
phase = exp(1i*angle(u1(1,1)));

u1 = u1 ./ exp(1i*angle(u1(1,1)));

[z,u2] = F2X3(conj(U1),0,M,e1,e2);

u2 = u2 ./ exp(1i*angle(u2(1,1)));
% [~,u2] = F2X3(U2,0,M,e1,e2);
% u2 = u2 ./ exp(1i*angle(u2(10,1)));
%%

P2 = (2*Dbar(N,0,f1,f2)')^2 - alpha^2 * W'^2;
norm(P2 * U1)
%% movie
fig=figure;
vid=VideoWriter('.\results\one_fourth_protected_1.mp4','MPEG-4'); open(vid)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

for alpha=0:0.02:2.99
% alpha=2;

    P = (2*Dbar(N,0,f1,f2))^2 - alpha^2 * W^2;
[~,~,U1] = svds(P, 1, 'smallest');
[z,u1] = F2X3(U1,0,M,e1,e2);
phase = exp(1i*angle(u1(1,1)));
u1 = u1 ./ phase;
disp(norm(imag(u1))/norm(u1))
% %%

clf(fig)
tl=tiledlayout(2,2,"TileSpacing","compact");

nexttile
hold on
title( ['$\rm{Re}\, u\,(\approx u)$ for $u \in \ker D(\alpha = ' sprintf('%.2f', alpha) ')$'], 'Interpreter', 'Latex')
contourf(real(z), imag(z), real(u1), 20)
% clim([-0.5 0.5])
hex(zS); axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63])
% clim([-1 1])
% %

% % 
% figure
% hold on
% title( ['$u \in \ker D(\alpha = 2)$'], 'Interpreter', 'Latex')
% contourf(real(z), imag(z), real(u2), 20)
% % clim([-0.5 0.5])
% hex(zS); axis equal; colorbar
% xlim([-0.63,0.63]), ylim([-0.63,0.63])
% % clim([-1 1])
% % %






% nexttile
% hold on
% contourf(real(z), imag(z), max(log(abs(u1)) / alpha,-1.2), 32)
% hex(zS); axis equal; colorbar
% xlim([-0.63,0.63]), ylim([-0.63,0.63])
% clim([-0.7 0.1])
% 
% nexttile
% hold on
% contourf(real(z), imag(z), max(log(abs(u2)) / alpha,-1.2), 32)
% hex(zS); axis equal; colorbar
% xlim([-0.63,0.63]), ylim([-0.63,0.63])
% clim([-0.7 0.1])

% Phase function
ww = sparse(2*N+1, 2*N+1); ww(N+1,N+1) = 1;
Phi = Inv(2*Dbar(N,0,f1,f2) + kron(ww,ww))*W * const_fun(N);
[z,phi] = F2X3(Phi,0,M,e1,e2);
phi = imag(phi);



a = max(phi, [], 'all') / pi;

% figure

bb = diag(u1);
bb = bb(250:350);
gg = diag(phi);
gg = cos(alpha*gg(250:350));
C = gg' * bb / (gg' * gg);


ff = diag(u1);
ff = ff(1:round(length(ff)/2));
len = length(ff);
ff = ff(round(len / 3): 2*round(len/3));
ff = real(ff);
gg = diag(phi);
gg = gg(1:round(length(gg)/2));
len = length(gg);
gg = gg(round(len / 3): 2*round(len/3));
gg = C*cos(alpha*gg);
nexttile
hold on
plot((1+(0:length(ff)-1)/(length(ff)-1))*abs(zS), ff, 'LineWidth', 2)
plot((1+(0:length(ff)-1)/(length(ff)-1))*abs(zS), gg, 'LineWidth', 2)
xlim([1/sqrt(3), 2/sqrt(3)]);
title('$u$ and $C(\alpha)\cos(\alpha\varphi)$ on the hexagon edge',...
    'Interpreter', 'latex')
legend(["$u$", "$C(\alpha)\cos(\alpha\varphi)$"], 'Interpreter', 'latex')
ylim([-0.8 0.8])




nexttile
hold on
title( ['$\cos(\alpha\varphi), \alpha = ' sprintf('%.2f', alpha) '$'], 'Interpreter', 'Latex')
contourf(real(z), imag(z), cos(alpha*phi), 10)
hex(zS); axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63])

nexttile
hold on
title( ['$|u(\alpha) - C(\alpha) \cos(\alpha \varphi)|$'], 'Interpreter', 'Latex')
contourf(real(z), imag(z), min(abs(real(u1) - C*cos(alpha*phi)),0.3), 32)
hex(zS); axis equal; colorbar   
xlim([-0.63,0.63]), ylim([-0.63,0.63])


drawnow
frame=getframe(fig);
writeVideo(vid,frame)
end

close(vid)
%% imaginary part for current potential

alphas = 0.1:0.1:3;

imagnorm = zeros(3, length(alphas));
Ns = [32 64];

for id =1:length(alphas)
    alpha = alphas(id)
    for id2=1:2
        tic
        N = Ns(id2);
        Up1 = sym_potential2(N,0,0,2);
        Um1 = Up1.';
        W = @(n) 1i*K^2 ...
    *(fourier_shift(N,-n,-n)+fourier_shift(N,n,n) ...
    + om*fourier_shift(N,n,0)+om*fourier_shift(N,-n,0) ...
    + om'*fourier_shift(N,0,n)+om'*fourier_shift(N,0,-n));
W = 0.5*(W(1)+W(2));
        P = (2*Dbar(N,0,f1,f2))^2 - alpha^2 * W^2;
        [~,s,U1] = svds(P, 1, 'smallest');
        disp(s)
        [~,u1] = F2X3(U1,0,M,e1,e2);
        phase = exp(1i*angle(u1(1,1)));
        u1 = u1 ./ phase;
        imagnorm(id2,id) = norm(imag(u1)) / norm(u1);
        toc 
    end
end

%% plot

figure
tl=tiledlayout(2,1,'TileSpacing', 'compact');

title(tl, "$P(\alpha) =(2D_{\bar z})^2 - \alpha^2 W^2$, for $W(z) =0.5(iU(z)U(-z)+iU(2z)U(-2z))$", ...
    'Interpreter', 'latex')

nexttile
hold on
title("$\|\rm{Im}\,u(\alpha)\|_2$ for $u \in \ker P(\alpha)$, $\|u\|_2 = 1$, $\rm{Im}\,u(z=0) = 0$, calculated with $N = 32$", 'Interpreter', 'latex')
plot(alphas, imagnorm(1,:), 'LineWidth', 2)

nexttile
hold on
title("difference in $\|\rm{Im}\,u(\alpha)\|_2$ when calculated with $N = 32, 64$", 'Interpreter', 'latex')
plot(alphas, imagnorm(1,:) - imagnorm(2,:), 'LineWidth', 2)

% nexttile
% hold on
% title("difference in $\|\rm{Im}\,u(\alpha)\|_2$ when calculated with $N = 64, 128$", 'Interpreter', 'latex')
% plot(alphas, imagnorm(2,:) - imagnorm(3,:), 'LineWidth', 2)

%% ratio

fig=figure;
vid=VideoWriter('.\results\ratio.mp4','MPEG-4'); open(vid)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

for alpha=0:0.02:2.99

    P = (2*Dbar(N,0,f1,f2))^2 - alpha^2 * W^2;
[~,~,U1] = svds(P, 1, 'smallest');
[z,u1] = F2X3(U1,0,M,e1,e2);
phase = exp(1i*angle(u1(1,1)));
u1 = u1 ./ phase;


clf(fig)

nexttile
hold on
title( ['$\rm{Re}\, u$ for $u \in \ker D(\alpha = ' sprintf('%.2f', alpha) ')$'], 'Interpreter', 'Latex')
contourf(real(z), imag(z), real(u1), 20)
% clim([-0.5 0.5])
hex(zS); axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63])




% Phase function
ww = sparse(2*N+1, 2*N+1); ww(N+1,N+1) = 1;
Phi = Inv(2*Dbar(N,0,f1,f2) + kron(ww,ww))*W * const_fun(N);
[z,phi] = F2X3(Phi,0,M,e1,e2);
phi = imag(phi);



a = max(phi, [], 'all') / pi;

% figure

bb = diag(u1);
bb = bb(250:350);
gg = diag(phi);
gg = cos(alpha*gg(250:350));
C = gg' * bb / (gg' * gg);


nexttile
hold on
title( ['$\cos(\alpha\varphi), \alpha = ' sprintf('%.2f', alpha) '$'], 'Interpreter', 'Latex')
contourf(real(z), imag(z), cos(alpha*phi), 10)
hex(zS); axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63])

nexttile
hold on
title( ['$|u(\alpha) - C(\alpha) \cos(\alpha \varphi)|$'], 'Interpreter', 'Latex')
contourf(real(z), imag(z), min(abs(real(u1) - C*cos(alpha*phi)),0.3), 32)
hex(zS); axis equal; colorbar   
xlim([-0.63,0.63]), ylim([-0.63,0.63])


drawnow
frame=getframe(fig);
writeVideo(vid,frame)
end

close(vid)
%% 2nd wkb

[~,aa] = F2X(W*const_fun(N),0,M,e1,e2);
[~,bb] = F2X(2*Dbar(N,0,f1,f2)*W*const_fun(N), 0, M,e1,e2);
rat = X2F(bb ./ aa, 0, N, e1,e2);
ww = sparse(2*N+1, 2*N+1); ww(N+1,N+1) = 1;
Phi1 = -0.5*Inv(2*Dbar(N,0,f1,f2) + kron(ww,ww))*rat;
[z,phi1] = F2X(Phi1, 0,M,e1,e2);
%%
figure
contourf(real(z), imag(z), real(phi1)/3, 100)
axis equal
colorbar
%%
figure
contourf(real(z), imag(z), imag(phi1), 32)
axis equal
colorbar

%%
[z,tr,logu]=triangle(10i,0.05);
logu(~tr) = NaN;

figure
contourf(real(z), imag(z), real(logu), 32)
axis equal
colorbar
%%

[xx,yy] = meshgrid(-0.5:0.01:0.5, -0.6:0.01:0.6);
z = xx + 1i*yy;

sv = SV(z);
ep = 0.02;
sv1 = SV(z + ep);
sv2 = SV(z + ep*1i);

dx = (sv1 - sv)/ep;
dy = (sv2 - sv)/ep;

dx(abs(dx) > 50) = 0;
dy(abs(dy) > 50) = 0;

figure
tiledlayout(2,2,"TileSpacing",'compact')

nexttile
hold on
contourf(real(z), imag(z), real(dx), 32)
hex(zS)
axis equal
colorbar

nexttile
hold on
contourf(real(z), imag(z), imag(dx), 32)
hex(zS)
axis equal
colorbar

nexttile
hold on
contourf(real(z), imag(z), real(dy), 32)
hex(zS)
axis equal
colorbar

nexttile
hold on
contourf(real(z), imag(z), imag(dy) + real(dx), 32)
hex(zS)
axis equal
colorbar

%%
% Fourier representation: a vector is a length (2N+1)^2 column vector V, represent 
% v(z) = \sum_{a,b\in{-N,...,N}} V((2N+1)*(b+N)+a+N+1) e^{i<z, a f1 + b f2+ k>} / sqrt(area_of_(e1,e2))

% Position-space representation: v(z) for
% z  a MxM matrix  with values   z(a,b) = e1*(a-1)/M + e2*(b-1)/M


% Dbar and multiply by potential in fourier basis

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

function U=sym_potential2(N,n1,n2,p)
    w = exp(2i*pi/3);
    U = (-4i*pi/3)* (fourier_shift(N,n1-1,n2+1) + w^p*fourier_shift(N,-n2-1,n1-n2) + w'^p*fourier_shift(N,n2-n1, -n1+1));
end

% Fourier representation of constant 1 function
% WARNING- this is for TBG lattice. for general lattice factor should be sqrt(abs(imag(e1' * e2)))
function V=const_fun(N)
    V = kron((-N:N).'==0, (-N:N).'==0) * sqrt(sqrt(3)/2); 
end


% Fourier transform

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

% Inverse
function V=X2F(v,k,N,e1,e2)
    M = size(v, 1);
    [y2,y1]=meshgrid(0:M-1, 0:M-1);
    z=(e1*y1+e2*y2)/M;

    v = v ./ exp(0.5i*(k'*z+k*conj(z)));
    
    V = fft2(v) * sqrt(abs(imag(e1' * e2))) / M^2;
    
    V = circshift(V, [N,N]);
    V = V(1:2*N+1, 1:2*N+1);
    V = V(:);
end

% For specific set of points
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

% Same as F2X but do for two layers (D(alpha) is 2x2 matrix), and correct
% for phase
function [z,v1,v2] = F2X2(V,k,M,e1,e2)
    K = 4*pi/3;

    [z,v1] = F2X(V(1:length(V)/2),k-K,M,e1,e2);
    [~,v2] = F2X(V(1+length(V)/2:end),k+K,M,e1,e2);

    phase=exp(1i * angle(v1(1)) );
    v1=v1 ./ phase; v2=v2./phase;
end

% Same as F2X, F2X2 but repeat over some copies of fundamental domain,
% and choose a phase

function [z,v] = F2X3(V,k,M,e1,e2, nophase)
    [z0,v0]=F2X(V,k,M,e1,e2);

    z = [z0-e1-e2, z0-e1; z0-e2, z0];
    
    p2 = exp(-0.5i*(k'*e2+k*e2')); p1 = exp(-0.5i*(k'*e1+k*e1'));
    v = [v0 *p1 *p2, v0 *p1; v0 *p2, v0];
    if ~exist('nophase', 'var')
        phase = exp(1i * angle(v0(1)));
        v = v ./  phase;
    end
end

function [z,v1,v2] = F2X4(V,k,M,e1,e2)
    K = 4*pi/3;

    [z,v1] = F2X3(V(1:length(V)/2),k-K,M,e1,e2);
    [~,v2] = F2X3(V(1+length(V)/2:end),k+K,M,e1,e2);
    phase=exp(1i*angle(v1(length(v1)/2)));
    v1=v1 ./ phase; v2=v2./phase;
end



% Symmetries of TBG

% For k in {0,K,-K}, C[u](x) = u(wx) acts on L^2_k(C;C^2)
function C = ROT2(N,s)         % k = sK
    C = [ROT(N,s-1),                 sparse((2*N+1)^2,(2*N+1)^2); 
        sparse((2*N+1)^2,(2*N+1)^2), ROT(N,s+1)];
end

% Project onto eigenspace with Cu = w^r u
function [V0,V1,V2] = PROJC(N,V,s)
    om = exp(2i*pi/3); R = ROT2(N,s);

    V0= (R^0 + om'^0 .* R + om^0 .* R^2)/3;
    V1= (R^0 + om'^1 .* R + om^1 .* R^2)/3;
    V2= (R^0 + om'^2 .* R + om^2 .* R^2)/3;
    V0 = V0 * V; V1 = V1 * V; V2 = V2 * V;
end

% Inclusion of V, L^2_{sK}(C;C^1) = V + RV + R^2V + C.
function Inc = ROT_INC(N,s,f1,f2)
    [yy,xx] = meshgrid(-N:N, -N:N);
    k = xx * f1 + yy * f2 + s*4*pi/3;
    k = k(:);
    vv = find( (abs(k) > 0.0001) & (angle(k) >= -pi/3-0.00001) & (angle(k) < pi/3 - 0.000001) );
    m = length(vv);
    Inc = sparse(vv, 1:m, ones(m,1), (2*N+1)^2, m);
end

function [Inc, Inc0] = ROT_INC2(N,s,f1,f2)
    [yy,xx] = meshgrid(-N:N, -N:N);
    k = xx * f1 + yy * f2 + s*4*pi/3;
    k = k(:);
    vv = find( (abs(k) > 0.0001) & (angle(k) >= -pi/3-0.00001) & (angle(k) < pi/3 - 0.000001) );
    m = length(vv);
    Inc = sparse(vv, 2:m+1, ones(m,1), (2*N+1)^2, m+1);
    Inc0 = sparse(find(abs(k) < 0.0001), 1, 1, (2*N+1)^2, m+1);
end


% For k in {0,K,-K}, R[u](x) = u(wx) acts on L^2_k(C;C^1)
function R = ROT(N,s)         % k = sK
    F1 = @(n1,n2) s+n2-n1;     % multiplication by w^-1 on k+Lambda^*
    F2 = @(n1,n2) -n1;
    R = PERMK(N,F1,F2);
end


% % r1: reflect on x axis
% % r2: reflect on y axis
% 
% % To act on L^2_k(C;C), must have R(k)-k = s1 e1 + s2 e2 for some integers s1,s2
% 
% function R = REF(N,s1,s2,r1,r2)
%     if mod(r1+r2,2)==0
%         fun1 = @(m1,m2) s1+m1*(-1)^r1;
%         fun2 = @(m1,m2) s1+m2*(-1)^r1;
%     else
%         fun1 = @(m1,m2) s1+m2*(-1)^r1;
%         fun2 = @(m1,m2) s2+m1*(-1)^r1;
%     end
%     R = PERMK(N,fun1,fun2);
% end

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


% Dbar and multiply by potential in position basis

% Next two functions are position space (finite element) representations
% of Dbar and V.
% Not exactly equivalent to fourier representations.

% Finite element for Dbar on periodic functions
function Db=FE_Dbar(M, f1, f2)
    D0 = spdiags(ones(M,1), 1, M, M); D0(M, 1) = 1;
    D1 = kron(speye(M, M), D0 - D0') * M/(2i);
    D2 = kron(D0 - D0', speye(M, M)) * M/(2i);

    Db = (1/2)*(f1 * D1 + f2 * D2) ./ (2*pi);
end

% Finite element for Dbar on functions with boundary conditions given by functions H1, H2
% H1 and H2 are holomorphic multipliers defining a line bundle:
%  sections are { u satisfying  u(z + e1) = H1(z) u(z),  u(z + e2) = H2(z) u(z)  }
% H1, H2 must satisfy compatibility condition
function Db=FE_Dbar2(M, f1, f2, e1, e2, H1, H2)
    A = spdiags(ones(M,1), 1, M, M); B = sparse(M,M); B(M,1) = 1;
    D1_front = kron(speye(M,M), A) + kron(spdiags(H1(e2*(0:M-1)' ./ M) ,0,M,M),   B);
    D1_back = kron(speye(M,M), A') + kron(spdiags(1./H1(e1*(-1)/M + e2*(0:M-1)' ./ M) ,0,M,M),   B');

    D2_front = kron(A, speye(M,M)) + kron(B, spdiags(H2(e1*(0:M-1)' ./ M) ,0,M,M));
    D2_back = kron(A', speye(M,M)) + kron(B', spdiags(1./H2(e2*(-1)/M + e1*(0:M-1)' ./ M) ,0,M,M));

    D1 = (D1_front - D1_back) * M/(2i);
    D2 = (D2_front - D2_back) * M/(2i);

    Db = (f1 * D1 + f2 * D2) ./ (4*pi);
end

% Multiply pointwise by V
function Mat=Mult(V, M)
    Mat = spdiags(V, 0, M^2, M^2);
end


% Miscellaneous


function result=theta0(z, tau)
    result = 0 * z;
    N = 200;
    for n=-N:N
        result = result + exp(pi*1i*n^2*tau + 2i*pi*n*z);
    end
end

function result=theta1(z, tau)
    result = 0 * z;
    N = 200;
    for n=-N:N
        result = result + exp(pi*1i*(n+0.5)^2*tau + 2i*pi*(n+0.5)*(z+0.5));
    end
end

function hex(v,height)
    % draw hexagon in current plot
    
    vv=exp((0:6) .* 2i*pi/6) * v;
    if ~exist('height', 'var')
        for k=1:6
                plot([real(vv(k)),real(vv(k+1))],[imag(vv(k)),imag(vv(k+1))], 'blue');
        end

    else
        for k=1:6
                line([real(vv(k)),real(vv(k+1))],[imag(vv(k)),imag(vv(k+1))], [height height], 'Color', 'blue');
        end
    
    end
end


function scattermult(A,dotsize) 
    % scatter with multiplicity
    
    [C,~,ic]=uniquetol(A, 0.001, 'ByRows', true);
    count=accumarray(ic,1); 
    set1=C(count==1,:); set2=C(count==2,:); set3=C(count==3,:);
    set4=C(count==4,:); set5=C(count>4,:);
    
    scatter(set1(:,1),set1(:,2),dotsize,'blue','filled');
    scatter(set2(:,1),set2(:,2),dotsize,'red','filled');
    scatter(set3(:,1),set3(:,2),dotsize,'black','filled');
    scatter(set4(:,1),set4(:,2),dotsize,'green','filled');
    scatter(set5(:,1),set5(:,2),dotsize,'magenta','filled');

end

function fQ=interpolate(e1,e2,f,M,xQ,yQ)
    % given function on fund domain (output of K2X), interpolate to z=(x,y)
    T = [real(e1), real(e2); imag(e1), imag(e2)];
    v = T^(-1) * [xQ; yQ];
    v = v - floor(v); v = ceil(v*M); v(v==0) = 1;
    ind = v(1,:) + M*(v(2,:)-1);
    fQ = f(ind);
end


function scattermult2(A,dotsize, color) 
    % scatter with multiplicity
    
    [C,~,ic]=uniquetol(A, 0.001, 'ByRows', true);
    count=accumarray(ic,1); 
    set1=C(count==1,:); set2=C(count==2,:); set3=C(count==3,:);
    set4=C(count==4,:); set5=C(count>4,:);
    
    scatter(set1(:,1),set1(:,2),dotsize,color,'filled');
    scatter(set2(:,1),set2(:,2),dotsize,color,'filled');
    scatter(set2(:,1),set2(:,2),3*dotsize,color);

    scatter(set3(:,1),set3(:,2),dotsize,color,'filled');
    scatter(set3(:,1),set3(:,2),3*dotsize,color);
    scatter(set3(:,1),set3(:,2),6*dotsize,color);

end

% from online
function map = wheelmap(n,sat,rev)
% 100x3 colormap that goes in color wheel order.  Uses hsv2rgb.
% colors can be circularly shifted, and order can be reversed, see below.
% bottom (row 1) to top (row 100) of colorbar is
% red -> yellow -> green -> cyan -> blue -> magenta->red
%
% optional inputs:
% integer n produces a circular shift in the colors.
% Positive vals circularly shift colors up on the colorbar, negative vals
% shift them down.  Default is 0.
% Values of n for particular colors at bottom and top: 
%   n = 0  red->red
%     -17  yellow->yellow
%     -33  green->green
%      50  cyan->cyan
%      33  blue->blue
%      17  magenta->magenta
%
% sat is the hsv saturation.  Default is 2/3.
%
% if rev is 'r' (in quotes) then colors go in reverse order.
% To bypass unneeded options use [], for example wheelmap([],[],'r')
%
% for angle(z) in the complex plane, the 2pi discontinuity does not appear
% and the colors go around z = 0 in counterclockwise order. 
% for a general complex function, use
% caxis([-pi pi])
% to bypass autoscaling.  Then branch cut discontinuities will be correct.
if nargin == 0
    n = 0; sat = 2/3; rev = 'x'; 
elseif nargin == 1
    sat = 2/3; rev = 'x';  
elseif nargin == 2
    rev = 'x'  ;
end
if isempty(n)
    n = 0;
end
if isempty(sat)
    sat  = 2/3;
end
h = linspace(0,1,100)';
s = sat*ones(size(h));
v = ones(size(h));
map = hsv2rgb([h s v]);
map = circshift(map,n,1);
if rev == 'r'
    map = flipud(map);
end
end