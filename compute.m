% Hi

% This is basically Professor Zworski's code that I got a long time ago.

% First let's set up some constants:

om = exp(2i*pi/3);
K = 4*pi/3;
zS = 1i/sqrt(3);

e1 = om^2;                     % basis of Lambda
e2 = -om;

f1 = (4i*pi/sqrt(3))*om;     % dual basis of Lambda^* (<e_i, f_j> = 2pi 1_{i==j})
f2 = (4i*pi/sqrt(3))*om^2;


N = 32;

% Potential

Up = sym_potential(N,0,0);
Um = Up.';
%for f_{-2}, use   Up2 = sym_potential(N,1,-1); Um2 = Up2.';

%% Potential function V

V = Um*Up*const_fun(N);
[z,v] = F2X3(V,0,300,e1,e2);

figure
tiledlayout(1,2,'TileSpacing','compact')

% |V|
nexttile
hold on 
contourf(real(z), imag(z), abs(v), 32)
hex(zS); axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63])

% angle(V)
nexttile
hold on
surf(real(z), imag(z), angle(v), 'EdgeColor', 'none');
axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63])
clim([-pi,pi]); colormap(gca,wheelmap);


figure
hold on
sv= sqrt(v);
zz=z(1:30:end,1:30:end); 
sv=sv(1:30:end,1:30:end);
quiver(real(zz(:)), imag(zz(:)), real(sv(:)), imag(sv(:)))
axis equal
xlim([-0.63,0.63]), ylim([-0.63,0.63])
hex(zS)

%% Magic angles using compact operator

% for chiral model
% Ak = Inv(2*Dbar(N, -K, f1, f2)) * Up * Inv(2*Dbar(N, +K, f1, f2)) * Um;

% for "scalar model"
Ak = Inv(4*Dbar(N,-K,f1,f2)^2) *  Up * Um;

Alphas = 1./sqrt(eigs(Ak, 400));
Alphas = [-Alphas; Alphas];



figure
hold on
scattermult2([real(Alphas), imag(Alphas)], 16, 'red')

% % Retrieve magic angles
% save_angles=load('.\angles\1_-2_0.0002_circle.mat', 'save_angles').('save_angles');


%% Protected state for scalar model

[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - 15^2*Up*Um, 1, 'smallest');
[z,u] = F2X3(V,0,300,e1,e2);
u = u ./ exp(1i*angle(u(1,1)));

figure
tiledlayout(1,2,'TileSpacing', 'compact')

% Plot log(abs(u_0))

nexttile
hold on
contourf(real(z), imag(z), abs(u), 32)
hex(zS); axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off

% Plot angle(u_0)

nexttile; hold on
surf(real(z), imag(z), angle(u), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
hex(zS,10);
view(2); axis equal;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off;
colorbar;

%% When is there solution for other rotation eigenspace at k=0 (scalar model)

% s: k = sK.
% t: rotation eigenspace u(\om z) = \om^{-t} u(z}.
% when s =/=0, you are looking at subset of magic angles.
% when s=0, only t=1 valid. then it's when there's an element of H^0 with
% the rotation symmetry t=1.
s=0;
t=1;

% s=1, t = 0, 2 coincide


R = ROT(N,s);
Inc = ROT_INC(N,s,f1,f2);

Proj = @(t) (R^0 + om^t*R^1 + om'^t * R^2) / 3;

D = Inc' * Dbar(N,s*K,f1,f2) * Inc; % rotation dbar
V = Inc'*(Up*Um)*(3*Proj(t))*Inc;   % rotation V

T = Inv(4*D^2) * V;

Alphas = 1./sqrt(eigs(T, 100));
Alphas = [Alphas; -Alphas];

figure
hold on
scattermult([real(Alphas), imag(Alphas)], 10)
yline(0)
axis equal

%% Position representation of P(alpha)

V = sym_potential(10,0,0) * sym_potential(10,0,0).' * const_fun(10);
M=300;
[z,v_pot] = F2X(V, 0, M, e1, e2);
alpha = 10;

P = 4*FE_Dbar(M, f1, f2)^2 - alpha^2 * Mult(v_pot(:), M);


[~,s,v] = svds(P, 3, 'smallest');
disp(diag(s))


w = v(:,3);
w = reshape(w, [M M]);
z = z((1:2:M), (1:2:M)); w = w((1:2:M), (1:2:M));

w = w / norm(w(:));
w= w * (M/2) / (sqrt(sqrt(3)/2));
w = w ./ exp(1i*angle(w(1,1)));


figure
tiledlayout(1,2,'TileSpacing','compact')

minlevel=-16;
levels=linspace(minlevel,1,28);

% Plot log(abs(u_0))

nexttile
hold on
contourf(real(z), imag(z), max(log(abs(w)), minlevel), levels)
hex(zS); axis equal; colorbar


% Plot angle(u_0)

nexttile; hold on
surf(real(z), imag(z), angle(w), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
hex(zS,10);
view(2); axis equal;
colorbar

%% Unsuccessfully trying to see -1 degree line bundle at magic angle

% Still unclear how to see it.

M = 400;
[z,up0] = F2X(sym_potential(10,0,0)*const_fun(10), -2*K, M, e1, e2);
[~,um0] = F2X(sym_potential(10,0,0).'*const_fun(10), 2*K, M, e1, e2);

k=0;
n = 1; % degree of line bundle we are tensoring with

Db_maker = @(k) FE_Dbar2(M, f1, f2, e1, e2, ...  % k = 0 --> zeros/poles at z=0
    @(z) 0*z+exp(1i*(k*e1'+k'*e1)/2 + n*1i*pi), @(z) 0*z+exp( 1i*(k*e2'+k'*e2)/2 +n*(-2*pi*1i*z/e1 -1i*pi*om) ));


Db_1 = Db_maker(k-K);
Db_2 = Db_maker(k+K);

alpha = 9.82907;
DD = [2*Db_1, alpha*Mult(up0(:),M); alpha*Mult(um0(:),M), 2*Db_2];
[~,s,vvv] = svds(DD, 8, 'smallest');
disp(s)


vv = vvv(:,7);
v1 = vv(1:length(vv)/2); 
v2 = vv(1+length(vv)/2:end);
v1 = reshape(v1, [M M]); 
v2 = reshape(v2, [M M]);

z = z((1:2:M), (1:2:M));
v1 = v1((1:2:M), (1:2:M));
v2 = v2((1:2:M), (1:2:M));


v1 = v1 ./ theta1(z./e1, om);
v1(1,1) = v1(1,2);
v1 = v1 / norm(v1);


v2 = v2 ./ theta1(z./e1, om);
v2(1,1) = v2(1,2);
v2 = v2 / norm(v2);


figure

tiledlayout(1,3,'TileSpacing','compact')
minlevel=-16;

nexttile
hold on
contourf(real(z), imag(z), max(log(abs(v1)), minlevel), 30)
axis equal; colorbar

nexttile
hold on
contourf(real(z), imag(z), max(log(abs(v2)), minlevel), 30)
axis equal; colorbar

nexttile; hold on
surf(real(z), imag(z), angle(v1), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
view(2)
axis equal
colorbar

%% Theta function

M = 300;
[z,up0] = F2X(const_fun(10), 0, M, e1, e2);
Th = theta1(z./e1, om);

figure
tiledlayout(1,2,'TileSpacing','compact')

nexttile
hold on
contourf(real(z), imag(z), max(log(abs(Th)),-5), 64)
axis equal
colorbar

nexttile
hold on
surf(real(z), imag(z), angle(Th), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
hex(zS,10);
view(2); axis equal;
colorbar

%% log |u_K| / alpha as alpha grows

vid=VideoWriter('.\results\log_norm_asymptotics_4.mp4','MPEG-4'); open(vid)
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

alphas = 20:0.1:20;

for id=1:length(alphas)
    clf(fig);
    tl=tiledlayout(2,2,'TileSpacing','compact');
    title(tl, ['$\alpha = ' sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex')

    % Get protected state
    [~,~,V] = svds((4*Dbar(N,0,f1,f2)^2 -alphas(id)^2 * Up*Um), 1, 'smallest');
    [z,v1] = F2X3(V,K,300,e1,e2);

    % Plot in hexagon
    nexttile(1, [2 1]);
    hold on
    title('$\alpha^{-1} \log |u|$', 'Interpreter', 'latex')
    minlevel=-1;
    levels=linspace(minlevel,0.25,40);
    contourf(real(z), imag(z), max(log(abs(v1)) ./ alphas(id), minlevel), levels)
    hex(zS); axis equal; colorbar
    xlim([-0.63,0.63]), ylim([-0.63,0.63]), clim([minlevel 0.25])

    % Plot on two lines
    MM = 600;
    [~,v1] = F2X(V,K,MM,e1,e2);
    
    nexttile;   
    hold on
    title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex')
    plot((0:MM-1)'./MM, log(abs(diag(v1))) ./ alphas(id))
    xlim([0 1]), ylim([-1 0.2])
    xline(1/3), xline(2/3), yline(0, 'Color', 'red')
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]

    nexttile;
    hold on
    title('restricted to line from $0$ to $\omega$', 'Interpreter', 'latex')
    plot((MM-1:-1:0)'./MM, log(abs(v1(1,:))) ./ alphas(id))
    xlim([0 1]), ylim([-1 1])
    xline(1/2), yline(0, 'Color', 'red')
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]

    drawnow
    frame=getframe(fig);
    writeVideo(vid,frame)
end

close(vid)

%% fourier transform of u_K restricted to edge of hexagon

% vid=VideoWriter('.\results\edge_fourier_asymptotics_3.mp4','MPEG-4'); open(vid);
fig=figure;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

alphas = 16:0.1:16;

for id=1:length(alphas)

    clf(fig);
    tl=tiledlayout(1,2,'TileSpacing','compact');
    title(tl, ['$\alpha = ' sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex');

    [~,~,V] = svds((4*Dbar(N,0,f1,f2)^2 -alphas(id)^2 * Up*Um), 1, 'smallest');
    M = 600; MM = 3*M;
    [~,v1] = F2X(V,K,MM,e1,e2);
    v0 = diag(v1); v0 = v0(M: 2*M-1);
    v0 = v0 / norm(v0); v0 = v0 / sign(v0(1));
    v0 = [v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0]; MM = 13*M;

    % fourier transform
    V0 = fft(v0);
    V0 = V0 * sqrt(alphas(id)/13);
    V0 = circshift(V0, round(MM/2));

    nexttile;   
    hold on
    title('protected state on line from $\frac{\sqrt{3}i}{3}$ to $\frac{2\sqrt{3}i}{3}$', 'Interpreter', 'latex');
    plot((0:(MM/13)-1)' ./(MM/13) .*(2/sqrt(3)), v0(1:MM/13));
    xlim([0 2/sqrt(3)]);
    ylim([-0.2 0.2]);
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 1 1]; % [width, height, depth]


    nexttile;
    hold on
    title('fourier transform', 'Interpreter', 'latex');
    xlabel('$\xi / \alpha$', 'Interpreter', 'latex');
    plot(((0:MM-1)' - round(MM/2)) .* (2*pi*sqrt(3)/(13*2*alphas(id))), abs(V0));
    xlim([-10 10]);
    ylim([-10 20]);
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 1 1]; % [width, height, depth]

%     drawnow;
%     frame=getframe(fig);
%     writeVideo(vid,frame);
end
% 
% close(vid);

%% log |u_K| / alpha as potential varies

vid=VideoWriter('.\results\log_norm_varying_potential_3.mp4','MPEG-4'); open(vid);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

Up2 = sym_potential(N,1,-1); Um2 = Up2.';
U2 = [0*Up, Up2; Um2, 0*Up];

alphas = @(id) 10;
thetas = 0:0.002:1;
theta = thetas * pi;

for id=1:length(theta)
    
    clf(fig)
    tl=tiledlayout(3,2,'TileSpacing','compact');
    title(tl, ['$U = \cos\theta f_1 + \sin\theta f_{-2},\,\,\theta = ' ...
        sprintf('%.3f', theta0(id)) '\pi$'], 'Interpreter', 'latex')
    
    % get protected state
    [~,~,V] = svds(D + alphas(id)*(cos(theta(id))*U + sin(theta(id))*U2) + K*speye(size(D)), 1, 'smallest');
    [z,v1,~] = F2X4(V,K,300,e1,e2);

    % points where max is achieved
    v1_max = max(abs(v1), [], 'all');
    max_inds = find(abs(abs(v1) - v1_max) < 0.0001);
    
    % graph in hexagon
    nexttile(1, [3 1])
    hold on
    title(['$\alpha^{-1} \log |u_K|,\quad\quad \alpha = '  sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex')
    minlevel=-1;
    levels=linspace(minlevel,0.25,40);
    contourf(real(z), imag(z), max(log(abs(v1)) ./ alphas(id), minlevel), levels)
    hex(zS); axis equal; colorbar
    plot(real(z(max_inds)), imag(z(max_inds)), '.', 'color', 'red', 'MarkerSize', 10)
    xlim([-0.63,0.63]), ylim([-0.63,0.63]), clim([minlevel 0.25])

    % graph on two lines
    MM = 600;
    [~,v1,~] = F2X2(V,K,MM,e1,e2);
    
    nexttile;   
    hold on
    title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex')
    plot((0:MM-1)'./MM, log(abs(diag(v1))) ./ alphas(id))
    xlim([0 1]), ylim([-1 0.2])
    xline(1/3), xline(2/3), yline(0, 'Color', 'red')
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]
    
    nexttile;
    hold on
    title('restricted to line from $0$ to $\omega$', 'Interpreter', 'latex')
    plot((MM-1:-1:0)'./MM, log(abs(v1(1,:))) ./ alphas(id))
    xlim([0 1]), ylim([-1 0.5]);
    xline(1/2), yline(0, 'Color', 'red')
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]
    
    % graph magic angles
    nexttile;
    Alphas = save_angles(:, 1+10*(id-1));
    
    hold on
    scattermult([real(Alphas), imag(Alphas)], 13)
    xlim([-0.1 11]), ylim([-8 8]), yline(0)

    drawnow;
    frame=getframe(fig);
    writeVideo(vid,frame)

end

close(vid)

%% does restriction of u_0 to edge of hexagon, approximately solve semiclassical ode?

alpha = 20;

[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2 * Up * Um, 1, 'smallest');
[~,v] = F2X(V,0,300,e1,e2);
phase = exp(1i*angle(v(1,1)));
V = V ./ phase;

Pot = Up*Um*const_fun(N);
DxV = Dx(N,0,f1,f2) * V;
DyV = Dy(N,0,f1,f2) * V;

DxxV = Dx(N,0,f1,f2) * DxV;
DxyV = Dy(N,0,f1,f2) * DxV;
DyyV = Dy(N,0,f1,f2) * DyV;


MM = 500;
[~,v] = F2X(V,0,3*MM,e1,e2);
[~,dxv] = F2X(DxV,0,3*MM,e1,e2);
[~,dyv] = F2X(DyV,0,3*MM,e1,e2);
[~,dxxv] = F2X(DxxV,0,3*MM,e1,e2);
[~,dxyv] = F2X(DxyV,0,3*MM,e1,e2);
[~,dyyv] = F2X(DyyV,0,3*MM,e1,e2);
[~,pot] = F2X(Pot,0,3*MM,e1,e2);


v = diag(v); dxv = diag(dxv).* (2 / alpha); dyv = diag(dyv).* (2 / alpha);
dxxv = diag(dxxv).* (2 / alpha)^2; dyyv = diag(dyyv).* (2 / alpha)^2; 
dxyv = diag(dxyv).* (2 / alpha)^2; pot = diag(pot);

%1 v      real
%2 Dxv    real   (Dx = -i dx)
%3 Dyv    imag
%4 DxDxv   real
%5 DyDyv   real
%6 DxDyv   imag
A = [real(v), real(dxv), imag(dyv), real(dxxv), real(dyyv), imag(dxyv), pot].';
names = ["u", "dxu", "dyu", "dxxu", "dyyu", "dxyu", "potential V"];

mm = 20;
A = A(:,MM+1+mm:2*MM-mm);

% B = normalize(A,2, "norm");
% ind = [1 2 4];
% [c,s,~] = svd(B(ind,:), "econ");
% 
% disp(names(ind))
% disp("singular values")
% disp(s)
% disp("coefficients")
% disp(c)
% 
% 
% figure
% hold on
% plot((mm:M-1-mm)'./M, abs(B(ind,:)))
% yline(0)
% xlim([0 1])
% legend(names(ind))

figure
hold on
title("$V$ vs second derivatives of $u$  on vertical edge of hexagon", 'Interpreter', 'Latex')
plot((mm:MM-1-mm)'./MM, A(4,:) ./ A(1,:))
plot((mm:MM-1-mm)'./MM, A(5,:) ./ A(1,:))
plot((mm:MM-1-mm)'./MM, A(6,:) ./ A(1,:))
plot((mm:MM-1-mm)'./MM, real(A(7,:)))
ylim([-100, 100])
legend(["$(4h^2 \partial_x^2 u) / u$", "$(4h^2 \partial_y^2 u) / u$", "$\Im (4h^2 \partial_x\partial_y u) / u$", "V"], 'Interpreter', 'Latex')


% disp( (imag(dyv)' * dxv) / (norm(dxv) * norm(dyv)) )
% disp( -norm(dxv) / norm(dyv) )



%% ... calculate c (coefficient in ode) for some values of alpha

alphas=5:0.1:20;
xx_yy = alphas * 0;
xx_xy = alphas * 0;
x_y = alphas * 0;
cor11 = alphas * 0;
cor22 = alphas * 0;

cor33 = alphas * 0;


for id=1:length(alphas)
   
alpha = alphas(id);
% alpha = 20;
disp(alpha)



[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2 * Up * Um, 1, 'smallest');
[~,v] = F2X(V,0,300,e1,e2);
phase = exp(1i*angle(v(1,1)));
V = V ./ phase;

Pot = Up*Um*const_fun(N);
DxV = Dx(N,0,f1,f2) * V;
DyV = Dy(N,0,f1,f2) * V;

DxxV = Dx(N,0,f1,f2) * DxV;
DxyV = Dy(N,0,f1,f2) * DxV;
DyyV = Dy(N,0,f1,f2) * DyV;


MM = 500;
[~,v] = F2X(V,0,3*MM,e1,e2);
[~,dxv] = F2X(DxV,0,3*MM,e1,e2);
[~,dyv] = F2X(DyV,0,3*MM,e1,e2);
[~,dxxv] = F2X(DxxV,0,3*MM,e1,e2);
[~,dxyv] = F2X(DxyV,0,3*MM,e1,e2);
[~,dyyv] = F2X(DyyV,0,3*MM,e1,e2);
[~,pot] = F2X(Pot,0,3*MM,e1,e2);


v = diag(v); dxv = diag(dxv).* (2 / alpha); dyv = diag(dyv).* (2 / alpha);
dxxv = diag(dxxv).* (2 / alpha)^2; dyyv = diag(dyyv).* (2 / alpha)^2; 
dxyv = diag(dxyv).* (2 / alpha)^2; pot = diag(pot);


A = [real(v), real(dxv), imag(dyv), real(dxxv), real(dyyv), imag(dxyv), pot].';
names = ["u", "dxu", "dyu", "dxxu", "dyyu", "dxyu", "potential V"];

mm = 20;
A = A(:,MM+1+mm:2*MM-mm);


ind = [2 3];
[c,d,~] = svd(A(ind,:), "econ");
x_y(id) = c(3)/c(4);
cor1 = corrcoef(A(ind, :)');
cor11(id) = cor1(2);
% disp(c)
% disp(d)

ind = [4 5];
[c,d,~] = svd(A(ind,:), "econ");
xx_yy(id) = c(3)/c(4);
cor2 = corrcoef(A(ind, :)');
cor22(id) = cor2(2);
% disp(c)
% disp(d)

ind = [4 6];
[c,d,~] = svd(A(ind,:), "econ");
xx_xy(id) = c(3)/c(4);
cor3 = corrcoef(A(ind, :)');
cor33(id) = cor3(2);
% disp(c)
% disp(d)

end

figure
tiledlayout(2,1,'TileSpacing', 'compact')

nexttile
hold on
title("correlation between derivatives restricted to the vertical edge", 'interpreter', 'latex')

plot(alphas, -cor11, 'LineWidth', 2)
plot(alphas, -cor22, 'LineWidth', 2)
plot(alphas, -cor33, 'LineWidth', 2)
legend(["$\partial_{yy} u$ and $\partial_{xx} u$",...
        "$\textrm{Im}\,\partial_{xy}  u $ and $\partial_{xx} u$", ...
        "$\partial_y u$ and $\textrm{Im}\,\partial_x u$"],...
       'Interpreter', 'Latex')


nexttile
hold on
title("approximate ratios between derivatives", 'interpreter', 'latex')
plot(alphas, sqrt(xx_yy), 'LineWidth', 2)
plot(alphas, xx_xy, 'LineWidth', 2)
plot(alphas, x_y, 'LineWidth', 2)

legend(["$\sqrt{\partial_{yy} u / \partial_{xx} u} \approx$",...
        "$\textrm{Im}\,\partial_{xy}  u / \partial_{xx} u \approx$", ...
        "$\partial_y u / \textrm{Im}\,\partial_x u \approx$"],...
       'Interpreter', 'Latex')

%% finite element for a hexagon dirichlet problem, Laplacian phi = {q,bar q}/|xi|^2
model = createpde;

rho = exp(1i*pi/3);
hex_z = zS * rho.^(0:5);
g = decsg([2 6 real(hex_z) imag(hex_z)]');
geometryFromEdges(model,g);

applyBoundaryCondition(model,"dirichlet", ...
                       "Edge",1:model.Geometry.NumEdges, ...
                       "u",0);

% f= {q, bar q} / |xi|^2
V = const_fun(N);
V = Up*Um*V;

dzV = 1i* Dbar(N,0,f1,f2)' * V;

MM=1000;
[~,v] = F2X(V,0,MM,e1,e2);
[~,dzv] = F2X(dzV,0,MM,e1,e2);
bracket = imag( sqrt(conj(v)) .* dzv );
f = abs(bracket) ./ abs(v);
f(1,1) = f(2,1); % lol

% % sanity checks
% figure
% hold on
% surf(real(z),imag(z),abs(bracket) ./ 1500, 'EdgeColor', 'none')
% view(2); axis equal; colorbar;
% figure 
% hold on
% contourf(real(z),imag(z),f, 18)
% view(2); axis equal; colorbar;


ffun = @(z,etc) interpolate(e1,e2,f,MM,z.x,z.y);
specifyCoefficients(model, "m",0,"d",0,"c",1,"a",0,"f",ffun);

generateMesh(model,"Hmax",0.01);
result=solvepde(model);

figure
hold on
pdeplot(result.Mesh, XYData=result.NodalSolution, Contour="on", ColorMap="default", Levels=18); 
axis equal
title("$\Phi$", 'Interpreter', 'latex')
% saveas(gcf,'./Phi_solution.png')

%% log eigenfunction / alpha   smoothened, then take laplacian

alpha = 20;
[~,~,V] = svds(D + alpha*U + K*speye(size(D)), 1, 'smallest');
V = V(1:length(V)/2); V = V / norm(V);
[z,v] = F2X(V,K,300,e1,e2);


rad = 8;
chi = ones(rad,rad); chi(300,300) = 0; %chi=circshift(chi, [-rad/2, -rad/2]);
v_conv =  ifft2(fft2(abs(v)) .* fft2(chi)) ./ rad^2;
v_conv = log(abs(v_conv)) ./ alpha;
% v_conv = max(log(abs(v)) ./ alpha, -0.5);

V_conv =X2F(v_conv,0,N,e1,e2);
[~,Lap_v_conv]=F2X(Dbar(N,0,f1,f2) * Dbar(N,0,f1,f2)' * V_conv, 0, 300, e1,e2);



figure
tiledlayout(1,4,'TileSpacing','compact');

minlevel=-1;
levels=linspace(minlevel,0.25,40);

nexttile;
hold on
title('$\chi$ where $\int \chi = 1$', 'Interpreter', 'latex');
contourf(real(z), imag(z), -chi ./ 2, levels);
axis equal; clim([minlevel 0.25]); 


nexttile;
hold on
title('$\alpha^{-1} \log |u_K|$', 'Interpreter', 'latex');
contourf(real(z), imag(z), max(log(abs(v)) ./ alpha, minlevel), levels);
clim([minlevel 0.25]); axis equal;

nexttile;
hold on
title('$\alpha^{-1} \log (|u_K| \ast \chi)$', 'Interpreter', 'latex');
contourf(real(z), imag(z), max(v_conv, minlevel), levels);
axis equal; clim([minlevel 0.25]); colorbar;

% nexttile;
% hold on
% title('$\Delta(\alpha^{-1} \log (|u_K| \ast \chi))$', 'Interpreter', 'latex');
% surf(real(z), imag(z), real(Lap_v_conv),'EdgeColor', 'None');
% axis equal; 
% view(2);
% % hex(zS); axis equal; colorbar;
% ylim([-1.721,0]), xlim([-0.5,0.5]); 
%  colorbar;

nexttile;
hold on
title('max$(-3,\Delta(\alpha^{-1} \log (|u_K| \ast \chi)) )$', 'Interpreter', 'latex');
contourf(real(z), imag(z), min(max(real(Lap_v_conv),-3), 18),15);
axis equal; 
view(2);
% hex(zS); axis equal; colorbar;
ylim([-1.721,0]), xlim([-0.5,0.5]); 
 colorbar;

% saveas(gcf,'./results/Laplacian_of_phi_1.png')

%% Phase function branch cut WKB (from AM's Phys Rev Letter paper on WKB)

MM = 300;

% Get scalar potential as a position space vector:
V = Up*Um*const_fun(N);
[z,v] = F2X(V,0,MM,e1,e2);


% Take square root
% Matlab sqrt returns number with  -pi <= arg <= pi
v = sqrt(-v);
v = 1i*v;

figure
hold on
contourf(real(z), imag(z), real(v), 32)
hex(zS)
axis equal
colorbar

figure
hold on
contourf(real(z), imag(z), imag(v), 32)
hex(zS)
axis equal
colorbar


% Take fourier transform
A = X2F(v,0,N,e1,e2); 

% Invert Dbar
Db = Dbar(N,0,f1,f2);
ww = sparse(2*N+1, 2*N+1); ww(N+1,N+1) = 1;
DD = 2i*Db + kron(ww, ww);

% 
A = Inv(DD) * A;
a0 = A(N*(2*N+1) + N+1) * sqrt(sqrt(3)/2);
A(N*(2*N+1) + N+1) = 0;


% transform back to position basis
[z,psi] = F2X(A,0,MM,e1,e2);
psi = psi + 0.5*(a0'*z+a0*conj(z));

figure
hold on
contourf(real(z), imag(z), real(psi), 32)
hex(zS)
axis equal
colorbar

figure
hold on
contourf(real(z), imag(z), imag(psi), 32)
hex(zS)
axis equal
colorbar

% Get 2 wkb solutions
[~,vv] = F2X(Up*Um*const_fun(N),0,MM,e1,e2); % potential v
vv = vv .^ (-0.25);
% vv(imag(z) > 0) = -1i*vv(imag(z)>0);


alpha = 15;
v1 = exp(1i*psi*alpha) .* vv;
v2 = exp(-1i*psi*alpha) .* vv;

%% ... wronskian with actual solution near edges

% V1 = X2K(v1, 0, N, e1, e2);
% V2 = X2K(v2, 0, N, e1, e2);


% protected state
[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2*Up*Um, 1, 'smallest');
[~,v_temp] = F2X(V,0,300,e1,e2);
V = V / exp(1i*angle(v_temp(1,1)));



% set of z to calculate for
step = 0.003;
[xx, yy] = meshgrid(-0.05:step:0.05, sqrt(3)/3-0.03:step:2*sqrt(3)/3+0.03);
z = xx + 1i*yy;

aa = F2Vect(A,z,f1,f2,32);
vv = F2Vect(Up*Um*const_fun(N),z,f1,f2,32);
vv = vv .^ -0.25;

u = F2Vect(V,z,f1,f2,32);
du = F2Vect(Dbar(N,0,f1,f2) * V,z,f1,f2,32);
u1 = exp(1i*aa*alpha) .* vv;
u2 = exp(-1i*aa*alpha) .* vv;
du1 = 0.5*(u1 - circshift(u1, [0 1]) + 1i* (u1 - circshift(u1, [1 0])) ) ./step;
du2 = 0.5*(u2 - circshift(u2, [0 1]) + 1i* (u2 - circshift(u2, [1 0])) ) ./step;

z = z(2:end, 2:end);
u = u(2:end, 2:end);
u1 = u1(2:end, 2:end);
u2 = u2(2:end, 2:end);
du = du(2:end, 2:end);
du1 = du1(2:end, 2:end);
du2 = du2(2:end, 2:end);
vv = vv(2:end, 2:end);


% u1 = K2Vect(V1,z,f1,f2,32);
% du1 = K2Vect(Dbar(N,0,f1,f2) * V1,z,f1,f2,32);
% u2 = K2Vect(V2,z,f1,f2,32);
% du2 = K2Vect(Dbar(N,0,f1,f2) * V2,z,f1,f2,32);


% wronskians
w1 = du1.*u - u1.*du ;
w2 = du2.*u - u2.*du ;
w3 = du2.*u1 - du1.*u2;

w1 = w1 ./ w3; w2 = w2 ./ w3;


figure
hold on
contourf(real(z), imag(z), min(max(imag(w1), -num), num), 100)
xline(0)
plot(0,2*sqrt(3)/3,'o')
plot(0,sqrt(3)/3,'o')
colorbar
axis equal

%% ... compare actual solution and wkb on edge

alpha =15;
xx = (sqrt(3)/3:0.002:2*sqrt(3)/3) * 1i;
aa = F2Vect(A,xx,f1,f2,32);
aa = aa + 0.5*(a0'*xx+a0*conj(xx));
vv = F2Vect(Up*Um*const_fun(N),xx,f1,f2,32);

u2 = cos(alpha * aa) .* (vv .^ -0.25);
u2 = min(u2, 1);
u2 = u2 / max(u2);

uu = F2Vect(V, xx, f1, f2, 32);
uu = uu / max(uu);


% plot
figure
hold on
plot(imag(xx), u2);
plot(imag(xx), uu);


xline(sqrt(3)/3)
xline(2*sqrt(3)/3)
yline(0)


% ... comparison of exponential growth rate near edge


alpha = 15;


% actual solution
[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2*Up*Um, 1, 'smallest');
fig=figure;
vid=VideoWriter('.\results\physics_growth_comparison_2.mp4','MPEG-4'); open(vid)
% set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

for s=-0.5*sqrt(3)*1/2:0.005:0.5*sqrt(3)*1/2
x = (-0.5:0.01:0.5) + sqrt(3)*1i/2 + 1i*s;

%
v = F2Vect(V,x,f1,f2,N);
psi = F2Vect(A,x,f1,f2,N);
psi = imag(psi);

clf(fig)
hold on
title("$\psi$ and $\log|u|$ for $\alpha = 15$", 'Interpreter', 'latex')
plot(x, log(abs(v)) / alpha)
plot(x, psi);
plot(x, -psi);

drawnow
frame=getframe(fig);
writeVideo(vid, frame)

end

close(vid)


%% comparison near zS with Airy function in w coordinates

% prologue -------------


% vid=VideoWriter('.\results\airy_holo_6.mp4','MPEG-4'); open(vid);
fig=figure;
% set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);


% alphas = 0.1:0.1:19.5;
% 
% for id=1:length(alphas)
% alpha = alphas(id);
alpha = 18;

% ----------------------
disp(id)
clf(fig);
tl=tiledlayout(2,3,'TileSpacing','compact');
title(tl, ['$\alpha = ' sprintf('%.1f', alpha) '$'], 'Interpreter', 'latex')

A = alpha^(-2/3) * -1i * 3^(1/3) / (2*pi);

% Range of w
[x,y] = meshgrid(-1:0.005:1, -1:0.005:1);
w = x + 1i*y; w= 5*w;
z = zS + A * w;

% Get protected state
M=4000;
[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - (alpha)^2*Up*Um, 1, 'smallest');
[~,v] = F2X(V,0,M,e1,e2);


apply = @(f) arrayfun(@(z) interpolate(e1,e2,f,M,real(z),imag(z)), z);
phase = exp(1i*angle(interpolate(e1,e2,v,M,real(zS), imag(zS) ) ) );




[~,dv] = F2X(1i*Dbar(N,0,f1,f2)*V, 0,M,e1,e2);

v = apply(v);
dv = apply(A'*dv); %%% ?????

dv = dv ./ phase;v = v ./ phase; 



% Airy function and derivatives
% c = (2/3)*(16/3)^(1/3);
c = 3^(-2/3) / 2;
% c=0;
ww = conj(w) - c* alpha^(-2/3) * w.^2;
%ww = conj(w);
ai0 = airy(ww);
ai1 = airy(om * ww);
ai2 = airy(om^2 * ww);
dai0 = airy(1,ww);
dai1 = airy(1,om * ww);
dai2 = airy(1,om^2 * ww);

F0 = ai0 + ai1 + ai2;
dF0 = dai0 + om*dai1 + om^2 * dai2;
F2 = ai0 + om^2 * ai1 + om*ai2;
dF2 = dai0 + dai1 + dai2;

% Compute wronskians
wronsF = F0 .* dF2 - F2 .* dF0; % 0.8270
% wronsF = F0 * 0 -3*sqrt(3)/(2*pi);
wrons0 = v .* dF0 - F0 .* dv;
wrons2 = v .* dF2 - F2 .* dv;

g0 = wrons2 ./ wronsF;
g2 = - wrons0 ./ (w .* wronsF);

% % compare
v2 = (wrons2 .* F0 - wrons0 .* F2) ./ wronsF;
disp( max(abs(v), [], 'all') )
disp( max(abs(v-v2), [], 'all') )

%
nexttile(1)
hold on
title("$\alpha^{-1}\log |u(w)|$", 'Interpreter', 'latex')
contourf(real(w), imag(w), max(log(abs(v)) ./ abs(alpha), -1.2), linspace(-1.2, 0.2, 50))
% contourf(real(w), imag(w), abs(v), 32)
axis equal
clim([-1.2, 0])
colorbar

nexttile(4)
hold on
title("$angle(u(w))$", 'Interpreter', 'latex')
surf(real(w), imag(w), angle(v), 'EdgeColor', 'none')
clim([-pi,pi]), colormap(gca,wheelmap)
view(2)
axis equal
xlim([-5, 5]), ylim([-5 5])
colorbar

nexttile(2)
hold on
title("$\alpha^{-1}\log |g_0(w)|$", 'Interpreter', 'latex')
contourf(real(w), imag(w), min(max(log(abs(g0)) ./ alpha, -1.2), 1),  linspace(-1.2, 1, 50))
axis equal
clim([-1.2 1])
colorbar

nexttile(5)
hold on
title("$angle(g_0(w^3))$", 'Interpreter', 'latex')
surf(real(w), imag(w), angle(g0), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
view(2)
axis equal
xlim([-5, 5]), ylim([-5 5])
colorbar

nexttile(3)
hold on
title("$\alpha^{-1}\log |g_2(w)|$", 'Interpreter', 'latex')
contourf(real(w), imag(w), min(max(log(abs(g2)) ./ alpha, -1.2), 1),  linspace(-1.2, 1, 50))
axis equal
clim([-1.2 1])
colorbar

nexttile(6)
hold on
title("$angle(g_2(w))$", 'Interpreter', 'latex')
surf(real(w), imag(w), angle(g2), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
view(2)
axis equal
xlim([-5, 5]), ylim([-5 5])
colorbar

% epilogue -------------
% drawnow
% frame=getframe(fig);
% writeVideo(vid,frame)
% end
% 
% close(vid)

% ----------------------


%% for some set of alpha, compute g_0(0) and g_0(2) from wronskian with Airy function


% prologue -------------

% [xx,yy] = meshgrid(0.1:0.15:15, 0.1:0.15:13);
% [xx,yy] = meshgrid(0.1:0.2:5, 0.1:0.5:5);
% alphas = xx + 1i * yy; alphas = alphas(:);
alphas = 0.1:0.1:18;
g00=0*alphas; g20 = 0*alphas;

for id=1:length(alphas)
alpha = alphas(id);
disp(alpha)
disp(id)

% ----------------------

A = -1i*(8 * pi^3 * alpha^2 / 3)^(-1/3);

% Range of w
[x,y] = meshgrid(-1:0.005:1, -1:0.005:1);
w = x + 1i*y; w= 5*w;
z = zS + A * w;

% Get protected state
M=400;
[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2*Up*Um, 1, 'smallest');
[~,v] = F2X(V,0,M,e1,e2);
[~,dv] = F2X(1i*Dbar(N,0,f1,f2)*V, 0,M,e1,e2);
phase = exp(1i*angle(v(1,1)));
dv = dv ./ phase;v = v ./ phase;


apply = @(f) arrayfun(@(z) interpolate(e1,e2,f,M,real(z),imag(z)), z);
v = apply(v);
dv = apply(A*dv);


% Airy function and derivatives
% c = 2^(2/3)*(pi/3)/(32*pi^3/3)^(1/3);
c = (2/3)*(16/3)^(1/3);
ww = conj(w) - c * alpha^(-2/3) * w.^2;

ai0 = airy(ww);
ai1 = airy(om * ww);
ai2 = airy(om^2 * ww);
dai0 = airy(1,ww);
dai1 = airy(1,om * ww);
dai2 = airy(1,om^2 * ww);

F0 = ai0 + ai1 + ai2;
dF0 = dai0 + om*dai1 + om^2 * dai2;
F2 = ai0 + om^2 * ai1 + om*ai2;
dF2 = dai0 + dai1 + dai2;

% Compute wronskians
wronsF = F0 * 0 -3*sqrt(3)/(2*pi);
wrons0 = v .* dF0 - F0 .* dv;
wrons2 = v .* dF2 - F2 .* dv;

g0 = wrons2 ./ wronsF;
g2 = - wrons0 ./ (w .* wronsF);

% % compare
% v2 = (wrons2 .* F0 - wrons0 .* F2) ./ wronsF;
% disp( max(abs(v), [], 'all') )
% disp( max(abs(v-v2), [], 'all') )


g00(id) = g0(201,201); g20(id) = g2(200,200);

% epilogue -------------

% toc
end
% save("./results/2d_g00.mat", "g00","g20", "alphas") % 4 alpha in 2d grid
% ----------------------

% for alpha on real line, graph (g_0(0), g_2(0)) as a point in R^2

% vid=VideoWriter('.\results\g00_g20_circle_2.mp4','MPEG-4'); open(vid);
% fig=figure;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% 
% 
% for id=2:length(alphas)
%     clf(fig)
%     hold on
%     title(['$(g_0(0), g_2(0))$, for $\alpha = ' sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex')
% 
% plot(real(g00), real(g20))
% xline(0); yline(0);
% plot(real(g00(id)), real(g20(id)), 'o', 'Color', 'red')
% axis equal
% 
%  drawnow;
%     frame=getframe(fig);
%     writeVideo(vid,frame);
% end
% 
% close(vid);

% --

% for alpha on real line, plot log(abs(g_0(0))), log(abs(g_2(0)))
% 
% fig=figure 
% hold on
% plot(alphas, log(g00), 'LineWidth', 2)
% % plot(alphas, log(g20))
% drawnow
% 
% % Ak = Inv(4*Dbar(N,K,f1,f2)^2) *  Up * Um;
% % Alphas = 1./sqrt(eigs(Ak, 500));
% % Alphas = Alphas(abs(imag(Alphas)) < 0.001 & 0 < real(Alphas) & real(Alphas) < 22);
% for id=1:length(Alphas)
%     xline(real(Alphas(id)));
% end
% 
% % legend(["$\log(g_0(0))$", "$\log(g_2(0))$", "magic $\alpha$"], 'Interpreter', 'Latex')
% legend(["$\log(g_0(0))$", "magic $\alpha$"], 'Interpreter', 'Latex', 'FontSize', 15)
% 
% % saveas(fig,'./results/g_00_2.png')

% --


% for alpha in 2d grid

% load("./results/2d_g00.mat", "g00","g20", "alphas")
% alphas=reshape(alphas, [87 100]);
% g00=reshape(g00, [87 100]);
% g20=reshape(g20, [87 100]);
% 
% figure
% tl=tiledlayout(1,2,'TileSpacing','compact');
% title(tl, '$|g_0(0)|$ and $|g_2(0)|$ as functions of $\alpha$', 'Interpreter', 'latex')
% 
% nexttile
% hold on
% title('$|g_0(0)|$', 'Interpreter', 'Latex')
% contourf(real(alphas), imag(alphas), log(abs(g00)), 30)
% colorbar
% axis equal
% 
% nexttile
% hold on
% title('$|g_2(0)|$', 'Interpreter', 'Latex')
% contourf(real(alphas), imag(alphas), log(abs(g20)), 30)
% colorbar
% axis equal 

% %saveas(gcf,'./results/g_00_plane_1.png')


%

%% for ODE, compare WKB and exact


% for start=0:0.02:1

% 
% V = @(x) exp(2i*pi*x);
% DV = @(x) 2i*pi*exp(2i*pi*x);
% SV = @(x) exp(1i*pi*x);

% ep = 2i;
% V = @(x) 1+ep*exp(2i*pi*x) + ep* exp(-2i*pi*x);
% DV = @(x) ep*2i*pi*exp(2i*pi*x) - ep* 2i*pi*exp(-2i*pi*x);
% SV = @(x) sqrt(V(x));


ep = 0.5i;
SV = @(x) 1+ep*exp(2i*pi*x) + ep* exp(-2i*pi*x);
DSV = @(x) ep*2i*pi*exp(2i*pi*x) - ep* 2i*pi*exp(-2i*pi*x);
V = @(x) SV(x).^2;
DV = @(x) 2*SV(x) .* DSV(x); 


start = 0;
dist = 1;

MM = 2500;
x = start+dist*(0:MM)/MM;

V = V(x); DV = DV(x); SV = SV(x);
for id=1:length(SV)-1
    if abs(SV(id+1) - SV(id)) > abs(SV(id))
        SV(id+1) = -SV(id+1);
    end
end



U = @(id, alpha) [1 1; SV(id) -SV(id)] * (eye(2) +0.125/alpha * DV(id) * SV(id)^(-3) * [0 -1; 1 0]);

l0 = sum(SV(1:end-1)) * dist/MM;
l1 = -0.25*sum(DV ./ V) * dist/MM;


l0_cum = cumsum(SV(2:end)) * dist/MM;
l1_cum = cumsum(-0.25 * DV(2:end) ./ V(2:end) ) * dist/MM;

alpha = 30i;

 
figure
hold on
plot(x, real(SV), 'Color', 'red')
plot(x, imag(SV), 'Color', 'blue')


exact_sol = zeros(4,MM);
wkb_sol = zeros(4,MM);
exact_sol(:,1) = [1 0 0 1]';
wkb_sol(:,1) = [1 0 0 1]';

for id=2:MM+1
        exact = reshape(exact_sol(:,id-1), [2 2]);
        exact = ( eye(2) + alpha*[0 1; V(id) 0]*(dist/MM) ) * exact;
        exact_sol(:,id) = reshape(exact, [4 1]);

        wkb = U(id,alpha) * [exp(alpha*l0_cum(id-1) + l1_cum(id-1)), 0; 0, exp(-alpha*l0_cum(id-1) + l1_cum(id-1))] * U(1, alpha)^(-1);
        wkb_sol(:,id) = reshape(wkb, [4 1]);
end



A = reshape( exact_sol(:,end) , [2 2] );
B = reshape( wkb_sol(:,end) , [2 2] );

R = U(length(x),alpha)^(-1) * A * U(1, alpha);
RR = [exp(alpha * l0 + l1), 0; 0, exp(-alpha*l0 + l1)];


% % Eigenfunction
% 
% % [V1,~] = eig(A);
% [V1,~] = eig(R); V1 = U(1,alpha) * V1;
% [V2,~] = eig(RR); V2 = U(1,alpha) * V2;
% 
% sol1 = arrayfun(@(a,b,c,d) [a c; b d] * V1, exact_sol(1,:), exact_sol(2,:), exact_sol(3,:), exact_sol(4,:), 'UniformOutput',false);
% sol2 = arrayfun(@(a,b,c,d) [a c; b d] * V2, wkb_sol(1,:), wkb_sol(2,:), wkb_sol(3,:), wkb_sol(4,:), 'UniformOutput',false);
% 
% get1 = @(sol) sol(1,1);
% get2 = @(sol) sol(1,2);

% end

%% compare for range of alpha

alphas = (1:0.1:20)*1i;

% prologue ---
% vid=VideoWriter('.\results\ODE_asymptotics_3.mp4','MPEG-4'); open(vid)
% fig=figure;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% 
% for len=0.8:0.01:1.9
% 
% clf(fig)
% tl=tiledlayout(1,2,'TileSpacing', 'compact');
% title(tl, ['interval [-0.2, ' sprintf('%.2f', -0.2 + len) ']'])
figure
tiledlayout(1,2,'TileSpacing', 'compact')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
% ---
     
M_exact = zeros(length(alphas), 2);
for ii=1:length(alphas)
    alpha = alphas(ii);
    M = eye(2);
    for jj=2:length(x)
        M = ( eye(2) + alpha*[0 1; V(jj) 0]*(dist/MM) ) * M;
    end
    M_exact(ii,:) = eig(M);
end

M_WKB = zeros(length(alphas), 2);
for ii=1:length(alphas)
    alpha = alphas(ii);
    M = U(length(x),alpha) * [exp(alpha*l0 + l1), 0; 0, exp(-alpha*l0 + l1)] * U(1, alpha)^(-1);
    M_WKB(ii,:) = eig(M);
end

nexttile
hold on
title('WKB')
plot(abs(alphas), real(log(M_WKB)),'LineWidth', 1)
plot(abs(alphas), imag(log(M_WKB)),'LineWidth', 1)
ylim([-15 15])
legend

nexttile
hold on
title('exact')
plot(abs(alphas), real(log(M_exact)),'LineWidth', 1)
plot(abs(alphas), imag(log(M_exact)),'LineWidth', 1)
ylim([-15 15])
legend

% epilogue ---
% drawnow
% frame=getframe(fig);
% writeVideo(vid,frame)
% end
% 
% close(vid)
% ---


%% Compare using spectral characterization of magic alpha


e1 = 1;                     % basis of Lambda
e2 = 1i;

f1 = 2*pi;     % dual basis of Lambda^* (<e_i, f_j> = 2pi 1_{i==j})
f2 = 2i*pi;



% 1
% Pot = 1*(0.5*fourier_shift(N,1,0)-0.5*fourier_shift(N,-1,0))+fourier_shift(N,0,0);
% Pot = Pot^2;

% 2
% Pot = 0.3*om*(0.5*fourier_shift(N,1,0)-0.5*fourier_shift(N,-1,0))+fourier_shift(N,0,0)...
%     + 0.3*fourier_shift(N,3,0) + 0.2*fourier_shift(N,4,0) + 0.1*fourier_shift(N,-3,0);
% Pot = Pot^2;


% 3
N=32;

Pot = 1i*(2*fourier_shift(N,1,0)-0.5*fourier_shift(N,-1,0))+fourier_shift(N,0,0);


% % % 4
%Pot = 0*fourier_shift(N,0,0) + 0.45*(fourier_shift(N,1,0) + fourier_shift(N,-1,0)) + 0close al*(fourier_shift(N,2,0) + fourier_shift(N,-2,0));
% Pot = Pot^2;

% 5
% Pot = 0.5*(0.5*fourier_shift(N,1,0)-0.5*fourier_shift(N,-1,0))+fourier_shift(N,0,0);
% Pot = Pot^4;

% Pot = fourier_shift(N,0,0) + 0.2i*(fourier_shift(N,1,0) + fourier_shift(N,-1,0) + fourier_shift(N,0,1) + fourier_shift(N,0,-1));


% figure
% hold on
% [z,pot] = F2X(Pot*const_fun(N), 0, 300, e1,e2);
% surf(real(z), imag(z), angle(pot), 'EdgeColor', 'none');
% view(2)
% colormap(gca,wheelmap)
% colorbar

% Pot = (0.5*fourier_shift(N,-1,0) + fourier_shift(N,1,0));
% 
k = 0.001;
Ak = Inv(4*Dbar(N,k,f1,f2)^2) * Pot;
Alphas = 1./sqrt(eigs(Ak, 300));

figure
hold on
scattermult([real(Alphas), imag(Alphas)], 15)
axis equal

%% Movie of the above

e1 = 1;                     % basis of Lambda
e2 = 1i;

f1 = 2*pi;     % dual basis of Lambda^* (<e_i, f_j> = 2pi 1_{i==j})
f2 = 2i*pi;

% Pot = (-0.3i*fourier_shift(N,-1,0)+0.3*fourier_shift(N,-2,0)+0.1*fourier_shift(N,3,0)+fourier_shift(N,0,0));



vid=VideoWriter('.\results\ode_angles_2.mp4','MPEG-4'); open(vid)
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);


for t=0:0.04:4
Pot = (0.5*t*fourier_shift(N,1,0)-0.5*t*fourier_shift(N,-1,0)+0.1*fourier_shift(N,2,0)+fourier_shift(N,0,0));

       clf(fig);
    title(['Spec$(2D_{\bar z}+0.01)^{-2} W^2$, where $W = 1+' sprintf('%.2f', t) 'i\sin(2\pi x)+0.1e^{4\pi x}$'], 'Interpreter', 'latex')

% figure
% hold on
% [z,pot] = F2X(Pot*const_fun(N), 0, 300, e1,e2);
% surf(real(z), imag(z), angle(pot), 'EdgeColor', 'none');
% view(2)
% colormap(gca,wheelmap)
% colorbar

% Pot = (0.5*fourier_shift(N,-1,0) + fourier_shift(N,1,0));
% 
k = 0.01;
Ak = Inv(4*Dbar(N,k,f1,f2)^2) * Pot^2;
Alphas = 1./sqrt(eigs(Ak, 600));
Alphas = [Alphas; -Alphas];

hold on
% scattermult([real(Alphas), imag(Alphas); real(-Alphas), imag(-Alphas)], 15)
scatter(real(Alphas), imag(Alphas), 15, 'filled')

axis equal
xlim([-45 45]) 
ylim([-45 45])



    drawnow
    frame=getframe(fig);
    writeVideo(vid,frame)
end

close(vid)



%% For global nonzero potential function, graph magic angles
% Pot = (fourier_shift(N,-3,-3)+fourier_shift(N,1,1) + fourier_shift(N,0,4)+5*fourier_shift(N,0,0));
% Pot = (-0.01*fourier_shift(N,-1,0)+fourier_shift(N,1,0));
Pot = (fourier_shift(N,-3,-3)+fourier_shift(N,1,1) + 2.1*fourier_shift(N,0,0));
% W = (fourier_shift(N,-1,0)+fourier_shift(N,0,-1)+fourier_shift(N,1,1) + 3.1*fourier_shift(N,0,0));
% Pot = W^2;

k=5;
V = Pot*const_fun(N);
[~,v] = F2X(V,0,300,e1,e2);
v= sqrt(v);
A0 = sum(v(:)) / 300^2;
Ak = Inv(4*Dbar(N,k,f1,f2)^2) * Pot;

Alphas = 1./sqrt(eigs(Ak, 300));
Alphas = [Alphas; -Alphas];

[xx,yy] = meshgrid(-10:10,-10:10);
AA = -xx*f1 + yy*f2; AA = AA(:); AA = (AA + k)/A0;
AA = [AA; -AA];

p = k/A0;
disp(Alphas(abs(Alphas - p) < 0.1) - AA(abs(AA - p) < 0.1));


figure 
hold on
% title('$V(z) = 2.1 + e^{2\pi i(a_1+a_2)} + e^{2\pi i (-2a_1 -3a_2)}$ for $z=a_1 \omega^2 - a_2 \omega$', 'Interpreter', 'latex')
% title('$V(z) = -e^{2\pi i(-a_1)} + 1.1e^{2\pi i a_1}$ for $z=a_1 \omega^2 - a_2 \omega$', 'Interpreter', 'latex')

scatter(real(Alphas), imag(Alphas), 5, 'filled', 'Color', 'blue')
plot(real(AA), imag(AA), 'o', 'Color', 'black')
% legend(["$-\frac{1}{\alpha^2} \in Spec_{L^2(\bf{C}/\Lambda; \bf{C})}(2D_{\bar z} + K)^{-2}V$", "$\frac{1}{A_0}(\pm K + \Lambda^*)$"], 'Interpreter', 'latex')
% legend(["$-\frac{1}{\alpha^2} \in Spec_{L^2(\bf{C}/\Lambda; \bf{C})}(2D_{\bar z} + K)^{-2}V$"], 'Interpreter', 'latex')

axis equal

%% ... graph eigenfunction for increasing alpha 

W = (fourier_shift(N,-1,0)+fourier_shift(N,0,-1)+fourier_shift(N,1,1) + 3.1*fourier_shift(N,0,0));
% W = (fourier_shift(N,-3,-3)+fourier_shift(N,1,1) + 2.1*fourier_shift(N,0,0));
% W = fourier_shift(N,0,0);

[~,w] = F2X(W*const_fun(N),0,300,e1,e2);
A0 = sum(w(:)) / 300^2;

%
% vid=VideoWriter('.\results\global_squareroot_eigenfunction_3.mp4','MPEG-4'); open(vid)
% fig=figure;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

fig=figure;
alpha = 30;

% alphas = 10:0.1:20;
% for alpha=alphas
% %

k = alpha * A0;
disp(k)
% %%
% D = 2*[Dbar(N,0,f1,f2), 0*W; 0*W, Dbar(N,0,f1,f2)] - alpha*[0*W, speye(size(W)); W^2, 0*W];
% k = eigs(D, 10, 'smallestabs');
% disp(k)
% % k = k(2);
% 
% %%
[~,s,U] = svds(4*Dbar(N,k,f1,f2)^2 - alpha^2*W^2, 1, 'smallest');
disp(s)
[~,u] = F2X3(U,k,300,e1,e2);

[~,s,U2] = svds(2*Dbar(N,k,f1,f2) - alpha*W, 1, 'smallest');
disp(s)
[z,u2] = F2X3(U2,k,300,e1,e2);

% Db = Dbar(N,0,f1,f2);
% ww = sparse(2*N+1, 2*N+1); ww(N+1,N+1) = 1;
% DD = 2*Db + kron(ww, ww);
% 
% Phi = Inv(DD) * W * const_fun(N);
% Phi(N*(2*N+1) + N+1) = 0;
% [~,phi] = F2X3(Phi, 0, 300, e1, e2, 'nophase');
% phi = alpha*phi + 1i*(z .* k' + conj(z) .* k)/2;

clf(fig)
tl=tiledlayout(1,2,'TileSpacing', 'compact');

title(tl, [['$\alpha =' sprintf('%.1f', alpha) '$']...
    "$k = \alpha (\int W) / (\int 1)$,   $W(a_1\omega^2 - a_2\omega) = 3.1 + e^{2\pi ia_1} + e^{2\pi i a_2} + e^{2\pi i(-a_1-a_2)}$" ...
    ], 'Interpreter', 'latex')


nexttile
hold on
title('$|\alpha|^{-1}\log|u|$, $u$ is lowest eigenfunction of $P^*P$, where $P(\alpha) = (2D_{\bar z}+k)^2 - \alpha^2 W^2$', 'Interpreter', 'latex');
minlevel=-0.8;
levels=linspace(minlevel,0.1,40);
contourf(real(z), imag(z), max(log(abs(u)) / alpha, minlevel), levels)
clim([-0.8, 0.1])
axis equal; colorbar
hex(zS)
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off

nexttile
hold on
title('$\log|u_1|$, $u_1 \in ker_{L^2(\bf{C}/ \Lambda; \bf{C})}2D_{\bar z}+k - \alpha W$', 'Interpreter', 'latex');

contourf(real(z), imag(z), log(abs(u2)), 32)
axis equal; colorbar
hex(zS)
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off


% %
% drawnow
% frame=getframe(fig);
% writeVideo(vid,frame)
% end
% 
% close(vid)
% %

%% ... graph eigenfunction/poisson bracket for changing potential

vid=VideoWriter('.\results\global_squareroot_potential_3.mp4','MPEG-4'); open(vid)
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

for ii=0:0.1:5

    disp("-----")

W = (fourier_shift(N,-1,0)+fourier_shift(N,0,-1)+ii*fourier_shift(N,1,1) + (2.1+ii)*fourier_shift(N,0,0));

[~,w] = F2X(W*const_fun(N),0,300,e1,e2);
A0 = sum(w(:)) / 300^2;


alpha = 13;


k2 = alpha * A0;
k = 0;
disp(k)

[~,s,U] = svds(4*Dbar(N,k,f1,f2)^2 - alpha^2*W^2, 1, 'smallest');
disp(s)
[~,u] = F2X3(U,k,300,e1,e2);

[~,s,U2] = svds(2*Dbar(N,k2,f1,f2) - alpha*W, 1, 'smallest');
disp(s)
[z,u2] = F2X3(U2,k2,300,e1,e2);


clf(fig)
tl=tiledlayout(1,3,'TileSpacing', 'compact');

title(tl, [[' $W(a_1\omega^2 - a_2\omega) = ' sprintf('%.1f', ii+2.1)  ' + e^{2\pi ia_1} + e^{2\pi i a_2} + ' sprintf('%.1f', ii)  'e^{2\pi i(-a_1-a_2)}$']...
    "$\alpha = 13, k = 0$"...%\alpha (\int W) / (\int 1)$" ...
    ], 'Interpreter', 'latex')


nexttile
hold on
title('$|\alpha|^{-1}\log|u|$,  $u$ lowest eigenfunction of $P^*P$', 'Interpreter', 'latex');
minlevel=-0.7 - ii*0.17;
levels=linspace(minlevel,0.2,40);
contourf(real(z), imag(z), max(log(abs(u)) / alpha, minlevel), levels)
clim([-0.7 - ii*0.17, 0.1])
axis equal; colorbar
hex(zS)
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off

nexttile
hold on
title('$\log|u_1|$, $u_1 \in ker_{L^2(\bf{C}/ \Lambda; \bf{C})}2D_{\bar z}+k_2 - \alpha W$', 'Interpreter', 'latex');

contourf(real(z), imag(z), log(abs(u2)), 32)
axis equal; colorbar
hex(zS)
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off



[z,w] = F2X3(W*const_fun(N),0,500,e1,e2);
[~,dzw] = F2X3(1i* Dbar(N,0,f1,f2)' * W*const_fun(N),0,500,e1,e2);

bracket = (abs(w).^2) .* imag( dzw );

nexttile
hold on

title("$|\{q, \bar q\}|_{q^{-1}(0)}$, for $q = 4\bar\zeta^2 - W^2$", 'Interpreter', 'latex')
surf(real(z), imag(z), abs(bracket), 'Edgecolor', 'None')
colorbar
axis equal
view(2)
surf(real(z), imag(z), 0*bracket,'Edgecolor', 'none', 'facecolor', 'red')
hex(zS, 10000)
xlim([-0.63,0.63]), ylim([-0.63,0.63]);


%
drawnow
frame=getframe(fig);
writeVideo(vid,frame)
end

close(vid)
%% Sample use of code for finding bracket area

figure
hold on

[boundary,unsure,areas,count,zz,bra]=bracket_area([0 0 1],600,3,4);

    
disp([areas, count]);

title('$|\{q, \bar q\}|$','Interpreter','latex');
surf(real(zz),imag(zz),abs(bra),'EdgeColor','none');

scatter3(real(boundary), imag(boundary), 1+ 0*boundary,  3,'red','filled');
scatter3(real(unsure), imag(unsure), 1+ 0*unsure, 4,'red','filled');

axis equal;

%% Eigenvalues of H: single alpha

H = [0*D, (D+5*U)'; (D+5*U), 0*D];
E = svds(D+5*U, 200, 'smallest');

figure hold on
plot(E, (1:200)');
plot(E, sqrt(3)/(4*pi) * E.^2);

% Compare with Weyl's law
figure hold on
plot(E, (sqrt(3)/(4*pi) * E.^2  - (1:200)' ) ./ E);
yline(0);

%% Eigenvalues of H for range of alpha

alphas = 0:0.1:15;

EE = zeros(length(alphas),40);
for id=1:length(alphas)
    disp(alphas(id))
    EE(id,:) = svds(D+alphas(id)*U, 40, 'smallest');
end

figure hold on
plot(alphas,EE,'Color','b');
xlabel("\alpha")
ylabel("Positive eigenvalues of $H_0(\alpha)$",'Interpreter','latex')

%% Graph log bands / alpha    as alpha grows

% Dbar
Db_1 = Dbar(N, -K, f1, f2);
Db_2 = Dbar(N, +K, f1, f2);


D = [2*Db_1, 0*Up; 0*Up, 2*Db_2];
U = [0*Up, Up; Um, 0*Up];

alphas=3:0.05:20;

E = zeros(20, length(alphas));

for id=1:length(alphas)

    disp(alphas(id))

    E(:,id) = log( svds(D + alphas(id)*U, 20, 'smallest') ) ./ alphas(id);

end

fig=figure; hold on
tl=tiledlayout(2,1,'TileSpacing', 'compact');
title(tl, '$k=0$', 'Interpreter', 'latex');

nexttile; hold on
title('$\alpha^{-1}\log E_j$ for $j=1,2,3$', 'Interpreter', 'latex');
plot(alphas, E(1:3, :));
xlabel('$\alpha$', 'Interpreter', 'latex');
ylim([-1, 0.5]); xlim([3, 20]);

nexttile; hold on
title('$\alpha^{-1}\log E_j$ for $j=1,\dots, 20$', 'Interpreter', 'latex');
plot(alphas, E(1:20, :));
xlabel('$\alpha$', 'Interpreter', 'latex');
ylim([-1, 0]); xlim([3, 20]);

saveas(fig,'./results/scaled_log_eigenvalues_k0_2.png')

%% Graph the first 10 eigenfunctions of H

p=10;
[~,s,V] = svds(D + 17*U + K*speye(size(D)), p, 'smallest');
disp([(1:p)', diag(s)])

fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
tiledlayout(2,5,'TileSpacing','tight');


minlevel=-1;
levels=linspace(minlevel,0.25,40);


for id=1:10
disp(id)
V2 = V(:,id);
[U0,U1,U2] = PROJC(N,V2,1);   % at least look symmetric bruh
disp([id, norm(U0), norm(U1), norm(U2)])
if norm(U0) > 10^(-5)
    V2 = U0;
elseif norm(U1) > 10^(-5)
    V2 = U1;
else
    V2 = U2;
end
V2 = V2 ./ norm(V2);

[z,v1,v2] = F2X4(V2,K,200,e1,e2);


nexttile; hold on
contourf(real(z), imag(z), max(log(abs(v1)) ./ 17,minlevel), levels);
hex(zS); axis equal; colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]);


end



%% Bands of scalar model

alphas = 2:0.01:10;
% alpha = alphas(abs(alphas - 3.00574) < 0.01);
% alphas = alpha + (-0.1:0.005:0.1);
% alpha = 2.9;
k_range = [-3*K/2:0.2:-0.5, -0.5:0.021:0.5  ,0.5:0.2:3*K/2];
bands = zeros(2, length(k_range));

vid=VideoWriter('.\results\bands_scalar_3.mp4','MPEG-4'); open(vid);
fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);


for alpha=alphas

    for id=1:length(k_range)

        P = 4*Dbar(N,k_range(id),f1,f2)^2 - alpha^2 * Up * Um;
        s = svds(P, 2, 'smallest');
        bands(:,id) = s;

    end

    clf(fig);
    hold on
    title(['$\alpha = $' sprintf('%.4f', alpha)], 'Interpreter', 'latex')
    plot(k_range, bands, 'LineWidth', 3);
    plot(k_range, -bands, 'LineWidth',3);
    xline(0)
    xline(K)
    xline(-K)
    xlim([-3*K/2 3*K/2])
    ym = max(bands(:));
    ylim([-ym*1.1,ym*1.1]);

    drawnow
    frame=getframe(fig);
    writeVideo(vid,frame)
end

close(vid)


%% Mengxuan multilayer n=2 bands

D = @(alpha, k) [2*Dbar(N,k-K, f1,f2), alpha*Up;
                 alpha*Um,  2*Dbar(N,k+K,f1,f2)];
ss = size(Up);
Z = sparse(ss(1), ss(2));
Tp = [speye(ss), Z;
    Z, Z];
Tm = [Z, Z;
    Z, speye(ss)];


% alpha = alphas(abs(alphas - 4.60233 - 3.399i) < 0.01);
% alphas = 0.3:0.01:1.7;
alpha = alphas(abs(alphas - 1.45282) < 0.01);


k_range = -3*K/2:0.2:3*K/2;
k_range2=K + (-0.1:0.01:0.1);
num=2;
num2=2;
bands = zeros(num, length(k_range));
bands2 = zeros(num, length(k_range2));

% vid=VideoWriter('.\results\bands_multilayer_2.mp4','MPEG-4'); open(vid);
fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);


% for alpha=alphas

    for id=1:length(k_range)

        % P = 4*Dbar(N,k_range(id),f1,f2)^2 - alpha^2 * Up * Um;
        D2 = [D(alpha, k_range(id)), Tp;
            Tm,  D(0, k_range(id))];
        s = svds(D2, num, 'smallest');
        bands(:,id) = s;

    end

    for id=1:length(k_range2)

        % P = 4*Dbar(N,k_range(id),f1,f2)^2 - alpha^2 * Up * Um;
        D2 = [D(alpha, k_range2(id)), Tp;
            Tm,  D(0, k_range2(id))];
        s2 = svds(D2, num2, 'smallest');
        bands2(:,id) = s2;

    end

    clf(fig);
    hold on
    tl=tiledlayout(1,2,'TileSpacing','compact');
    title(tl,['$\alpha = $' sprintf('%.2f', alpha)], 'Interpreter', 'latex')

    nexttile
    hold on
    plot(k_range, bands);
    plot(k_range, -bands);
    xline(0)
    xline(K)
    xline(-K)
    ym = max(bands,[],'all');
    ylim([-ym*1.1,ym*1.1]);

    nexttile
    hold on
    plot(k_range2, bands2);
    plot(k_range2, -bands2);
    xline(K)
    ym = max(bands2,[],'all');
    ylim([-ym*1.1,ym*1.1]);



    % drawnow
    % frame=getframe(fig);
    % writeVideo(vid,frame)
% end

% close(vid)

%% FUNCTIONS ----------------------------------

% Fourier representation: a vector is a length (2N+1)^2 column vector V, represent 
% v(z) = \sum_{a,b\in{-N,...,N}} V((2N+1)*(b+N)+a+N+1) e^{i<z, a f1 + b f2+ k>} / sqrt(area_of_(e1,e2))

% Position-space representation: v(z) for
% z  a MxM matrix  with values   z(a,b) = e1*(a-1)/M + e2*(b-1)/M


% Fourier basis operators

function Dbar = Dbar(N,k,f1,f2)
    D0 = spdiags((-N:1:N)', 0, 2*N+1, 2*N+1); E = speye(2*N+1, 2*N+1);
    Dbar = 0.5*(kron(E, f1 * D0) + kron(f2 * D0, E) + k*kron(E,E));
end

function A = Inv(B) % invert a diagonal matrix
    n = size(B,1);
    A = spdiags(1./diag(B), 0, n, n);
end

function U=fourier_shift(N,n1,n2)
    N = 2*N + 1;
    U = kron(spdiags(ones(N,1), -n2, N, N),  spdiags(ones(N,1), -n1, N, N));
end

% (2K+) Potential with TBG rotation symmetry containing K + n1 f1 + n2 f2
function U=sym_potential(N,n1,n2)
    w = exp(2i*pi/3);
    U = (-4i*pi/3)* (fourier_shift(N,n1-1,n2+1) + w*fourier_shift(N,-n2-1,n1-n2) + w'*fourier_shift(N,n2-n1, -n1+1));
end

% Fourier representation of constant 1 function
% This is for TBG lattice. for general lattice factor should be sqrt(abs(imag(e1' * e2)))
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

% For a specific set of points
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


% Same as F2X but repeat over some copies of fundamental domain,
% and choose a phase

function [z,v] = F2X3(V,k,M,e1,e2)
    [z0,v0]=F2X(V,k,M,e1,e2);

    z = [z0-e1-e2, z0-e1; z0-e2, z0];
    
    p2 = exp(-0.5i*(k'*e2+k*e2')); p1 = exp(-0.5i*(k'*e1+k*e1'));
    v = [v0 *p1 *p2, v0 *p1; v0 *p2, v0];
end


% Symmetry operators

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


% Position basis operators

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


% Misc -----------------------


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


function scattermult2(A,dotsize, color)     
    [C,~,ic]=uniquetol(A, 0.001, 'ByRows', true);
    count=accumarray(ic,1); 
    for num=1:10
        alphas = C(count==num,:);
        scatter(alphas(:,1), alphas(:,2), dotsize, color, 'filled')
        for jj=1:num-1
            scatter(alphas(:,1), alphas(:,2), 3*jj*dotsize, color)
        end
    end
end

function fQ=interpolate(e1,e2,f,M,xQ,yQ)
    % given function on fund domain (output of K2X), interpolate to z=(x,y)
    T = [real(e1), real(e2); imag(e1), imag(e2)];
    v = T^(-1) * [xQ; yQ];
    v = v - floor(v); v = ceil(v*M); v(v==0) = 1;
    ind = v(1,:) + M*(v(2,:)-1);
    fQ = f(ind);
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