% Hexagonal Lattice %

w = exp(2i*pi/3);
K = 4*pi/3;
zS = 1i/sqrt(3);

e1 = w^2;                     % basis of Lambda
e2 = -w;

f1 = (4i*pi/sqrt(3))*w;     % dual basis of Lambda^* (<e_i, f_j> = 1_{i==j})
f2 = (4i*pi/sqrt(3))*w^2;



% Hamiltonian %

N = 32;

% Dbar

Db_1 = Dbar(N, -K, f1, f2);
Db_2 = Dbar(N, +K, f1, f2);

% Potential

Up = sym_potential(N,0,0);
Um = Up.';


D = [2*Db_1, 0*Up; 0*Up, 2*Db_2];
U = [0*Up, Up; Um, 0*Up];




%% Calculate magic angles with T_k

Ak = Inv(2*Db_1) * Up * Inv(2*Db_2) * Um;

Alphas = 1./sqrt(eigs(Ak, 500));

figure; 
hold on;
scattermult([real(Alphas), imag(Alphas)], 5);



% Gap of real magic angles
RealAlphas = Alphas(abs(imag(Alphas)) < 0.01);

fprintf('Gap of real magic angles ~ %d\n', RealAlphas(6)-RealAlphas(5));


%% Eigenfunctions


[~,s,V] = svds(D + 20*U + K*speye(size(D)), 1, 'smallest');
[z,v1,v2] = K2X4(V,K,200,e1,e2);

fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
til=tiledlayout(1,2,'TileSpacing','compact');
minlevel=-16;
levels=linspace(minlevel,1,28);


nexttile; hold on;
contourf(real(z), imag(z), max(log(abs((v1))),minlevel), levels);
hex(zS); axis equal; colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off;

nexttile; hold on;
contourf(real(z), imag(z), max(log(abs((v2))),minlevel), levels);
hex(zS); axis equal; colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off;



%% Eigenvalues of H: fixed alpha

H = [0*D, (D+5*U)'; (D+5*U), 0*D];

E = svds(D+5*U, 200, 'smallest');

figure; hold on;
plot(E, (1:200)');
plot(E, sqrt(3)/(4*pi) * E.^2);

figure; hold on;
plot(E, (sqrt(3)/(4*pi) * E.^2  - (1:200)' ) ./ E);
yline(0);


%% Eigenvalues of H: varying alpha

alphas = 0:0.1:15;

EE = zeros(length(alphas),40);
for id=1:length(alphas)
    disp(alphas(id))
    EE(id,:) = svds(D+alphas(id)*U, 40, 'smallest');
end

figure; hold on;
plot(alphas,EE,'Color','b');
xlabel("\alpha")
ylabel("Positive eigenvalues of $H_0(\alpha)$",'Interpreter','latex')


%% Potential U

inp = @(x,y) (x'*y+x*y')/2;
U_fun = @(x) exp(1i*inp(x,K)) + w*exp(1i*inp(x,w*K)) +  w^2*exp(1i*inp(x,w^2*K));

V_fun = @(x) U_fun(x) * U_fun(-x);

V = @(x) abs(exp(1i*x*K) - exp(-1i*x*K/2));
xx = 0:0.001:1;
VV = arrayfun(V, xx);
disp(sum(VV) / length(VV));


%% Bracket

figure; hold on;

[boundary,unsure,areas,count,zz,bra]=bracket_area([0 0 1],600,3,4);
disp([areas, count]);

title('$|\{q, \bar q\}|$','Interpreter','latex');
surf(real(zz),imag(zz),abs(bra),'EdgeColor','none');

scatter3(real(boundary), imag(boundary), 1+ 0*boundary,  3,'red','filled');
scatter3(real(unsure), imag(unsure), 1+ 0*unsure, 4,'red','filled');

axis equal;


%% Animation comparing lowest eigenstate at k=K,0, and 3, as \alpha \rightarrow \infty

vid=VideoWriter('.\results\different_k_state_2.mp4','MPEG-4'); open(vid);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);


[~,~,~,~,zz,bra]=bracket_area([0 0 1],600,3,4);


alphas = 0:0.01:12;
for id=1:length(alphas)

    clf(fig);
    tl=tiledlayout(2,3,'TileSpacing','compact');
    title(tl, ['$\alpha = ' sprintf('%.2f', alphas(id))  '$'], 'Interpreter', 'latex');

    [~,s1,V1] = svds(D + alphas(id)*U + 0*speye(size(D)), 5, 'smallest');
    [~,s2,V2] = svds(D + alphas(id)*U + 2*speye(size(D)), 5, 'smallest');
    [~,s3,V3] = svds(D + alphas(id)*U + K*speye(size(D)), 1, 'smallest');
    V1 = V1(:,5); V2 = V2(:,5); V3 = V3(:,1);
    s1 = diag(s1); s2=diag(s2); s3 = diag(s3);
    V10=PROJC(N,V1,0,0);V11=PROJC(N,V1,0,1);V12=PROJC(N,V1,0,2);
    sizes = '0:' +  string(norm(V10) >0.001) + ', 1:' +  string(norm(V11) >0.001) + ', 2:' +  string(norm(V12) >0.001);
    if norm(V10)> 0.001
        V1 = V10;
    elseif norm(V11) > 0.001
        V1 = V11;
    else
        V1 = V12;
    end
    [z,v11,v12] = K2X4(V1,K,200,e1,e2);
    [~,v21,v22] = K2X4(V2,0,200,e1,e2);
    [~,v31,v32] = K2X4(V3,3,200,e1,e2);


    minlevel=-8.5;
    levels=linspace(minlevel,0.5,28);


    nexttile; hold on;
    contourf(real(z), imag(z), max(log(abs(v11)),minlevel), levels);
    title(['$k=0,\, E=' sprintf('%.3f', s1(5)) '$'],'Interpreter','latex');
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 1]); hold off;

    nexttile; hold on;
    contourf(real(z), imag(z), max(log(abs(v21)),minlevel), levels);
    title(['$k=2,\, E=' sprintf('%.3f', s2(5)) '$'],'Interpreter','latex');
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 1]); hold off;

    nexttile; hold on;
    contourf(real(z), imag(z), max(log(abs(v31)),minlevel), levels);
    title(['$k=K,\, E=' sprintf('%.3f', s3(1)) '$'],'Interpreter','latex');
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 1]); hold off;


    nexttile; hold on;
    title(sizes);
    surf(real(zz),imag(zz),abs(bra),'EdgeColor','none');
    axis equal; xlim([-0.63,0.63]), ylim([-0.63,0.63]); 
    hold off;

    nexttile; hold on;
    S = [[0*s1; 0*s2+1; 0*s3+2], [s1; s2; s3]];
    [C,~,ic]=uniquetol(S, 0.001, 'ByRows', true);
    count=accumarray(ic,1);
    S1 = C(count==1,:); S2 = C(count==2,:); S3 = C(count==3,:); S4 = C(count>=4,:);
    scatter(S1(:,1),S1(:,2),8,'black','filled');
    scatter(S2(:,1),S2(:,2),10,'blue','filled');
    scatter(S3(:,1),S3(:,2),12,'red','filled');
    scatter(S4(:,1),S4(:,2),14,'red','filled');




    drawnow;
    frame=getframe(fig);
    writeVideo(vid,frame);

end

close(vid);



%%



%% FUNCTIONS --------------------------------------------------------------

% Below we provide functions related to the space L^2_k(C/(Z e1 + Z e2); C) 
% (equivalently L^2(Z f1 + Z f2 + k; C) ).

% Some functions will also be specifically for TBG.



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


function U=fourier_shift(N,n1,n2)

    N = 2*N + 1;

    U = kron(spdiags(ones(N,1), -n1, N, N),  spdiags(ones(N,1), -n2, N, N));

end

% Potential with rotation symmetry for TBG

function U=sym_potential(N,n1,n2)

    w = exp(2i*pi/3);

    U = (-4i*pi/3)* (fourier_shift(N,n1-1,n2+1) + w*fourier_shift(N,-n2-1,n1-n2) + w^2*fourier_shift(N,n2-n1, -n1+1));
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

function [z,v1,v2] = K2X2(V,k,M,e1,e2)

    K = 4*pi/3;

    [z,v1] = K2X(V(1:length(V)/2),k-K,M,e1,e2);
    [~,v2] = K2X(V(1+length(V)/2:end),k+K,e1,e2);

    phase=v1(1) / abs(v1(1));
    v1=v1 ./ phase; v2=v2./phase;

end

% Next two functions are for graphing

function [z,v] = K2X3(V,k,M,e1,e2)

    [z0,v0]=K2X(V,k,M,e1,e2);

    w = exp(2*pi*1i/3);
    z = [z0-e1-e2, z0-e2; z0-e1, z0];
    
    p2 = exp(-0.5i*(k'*e2+k*e2')); p1 = exp(-0.5i*(k'*e1+k*e1'));
    v = [v0 *p1 *p2, v0 *p2; v0 *p1, v0];

end

function [z,v1,v2] = K2X4(V,k,M,e1,e2)

    K = 4*pi/3;

    [z,v1] = K2X3(V(1:length(V)/2),k-K,M,e1,e2);
    [~,v2] = K2X3(V(1+length(V)/2:end),k+K,M,e1,e2);

    phase=v1(length(v1)/2) / abs(v1(length(v1)/2));
    v1=v1 ./ phase; v2=v2./phase;

end






% Symmetries of TBG

% For k in {0,K,-K}, C[u](x) = u(wx) acts on L^2_k(C;C^2)

function C = ROT2(N,s)         % k = sK

    C = [ROT(N,s-1),                 sparse((2*N+1)^2,(2*N+1)^2); 
        sparse((2*N+1)^2,(2*N+1)^2), ROT(N,s+1)];

end

% Project onto eigenspace with Cu = w^r u

function V2 = PROJC(N,V,s,r)

    w = exp(2i*pi/3); R = ROT2(N,s);

    M= (R^0 + w'^r .* R + w^r .* R^2)/3;
    V2 = M* V;

end


% For k in {0,K,-K}, R[u](x) = u(wx) acts on L^2_k(C;C^1)

function R = ROT(N,s)         % k = sK

    f1 = @(m1,m2) s+m2-m1;     % multiplication by w^-1 on k+Lambda^*
    f2 = @(m1,m2) -m1;
    R = PERMK(N,f1,f2);

end



% r1: reflect on x axis
% r2: reflect on y axis

% To act on L^2_k(C;C), must have R(k)-k = s1 e1 + s2 e2 for some integers s1,s2

function R = REF(N,s1,s2,r1,r2)

    if mod(r1+r2,2)==0
        fun1 = @(m1,m2) s1+m1*(-1)^r1;
        fun2 = @(m1,m2) s1+m2*(-1)^r1;
    else
        fun1 = @(m1,m2) s1+m2*(-1)^r1;
        fun2 = @(m1,m2) s2+m1*(-1)^r1;
    end
    R = PERMK(N,fun1,fun2);

end


% Operator permuting fourier basis:

% (FV)(F1(m1,m2),F2(m1,m2)) = V(m1,m2)  

function F = PERMK(N,F1,F2)

    indx =  @(m1,m2) (2*N+1)*(m1+N) + m2+N + 1;

    [m1,m2]=meshgrid(-N:N,-N:N); 
    m1=m1(:); m2=m2(:);

    fm1 = F1(m1,m2);
    fm2 = F2(m1,m2);

    mask = fm1 <= N & fm1 >= -N & fm2 <= N & fm2 >= -N;

    F = sparse(indx(fm1(mask),fm2(mask)), indx(m1(mask),m2(mask)), ones(length(find(mask)),1), (2*N+1)^2, (2*N+1)^2);

end

