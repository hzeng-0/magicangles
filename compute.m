% Hi,

% This is basically Professor Zworski's code (that I got a long time ago.)
% We represent our operators in the fourier basis.

% The code is separated into sections- in the Matlab editor, type Ctrl-Enter to run a section.


% First let's set up some constants:

om = exp(2i*pi/3);
K = 4*pi/3;
zS = 1i/sqrt(3);

e1 = om^2;                     % basis of Lambda
e2 = -om;

f1 = (4i*pi/sqrt(3))*om;     % dual basis of Lambda^* (<e_i, f_j> = 2pi 1_{i==j})
f2 = (4i*pi/sqrt(3))*om^2;


% Now, for example, we can make the Hamiltonian for non-scalar model.
% (Using the functions at the end of this file.)

N = 32;

% Dbar

Db_1 = Dbar(N, -K, f1, f2);
Db_2 = Dbar(N, +K, f1, f2);

% Potential

Up = sym_potential(N,0,0);
Um = Up.';
%for f_{-2}, use   Up2 = sym_potential(N,1,-1); Um2 = Up2.';

D = [2*Db_1, 0*Up; 0*Up, 2*Db_2];
U = [0*Up, Up; Um, 0*Up];


%% Retrieve magic angles

save_angles=load('.\angles\1_-2_0.0002_circle.mat', 'save_angles').('save_angles');


%% Calculate magic angles with T_k

% for non-scalar (original) model
Ak = Inv(2*Db_1) * Up * Inv(2*Db_2) * Um;

% for scalar model
% Ak = Inv(4*Dbar(N,K,f1,f2)^2) *  Up * Um;

Alphas = 1./sqrt(eigs(Ak, 500));

figure 
hold on
scattermult([real(Alphas), imag(Alphas)], 5)

%% Protected state for scalar model

[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - 10^2*Up*Um, 1, 'smallest');
[z,v] = K2X3(V,0,600,e1,e2);
v = v ./ exp(1i*angle(v(1,1)));

figure
tiledlayout(1,2,'TileSpacing','compact')
minlevel=-16;
levels=linspace(minlevel,1,28);

% Plot log(abs(u_0))

nexttile
hold on
contourf(real(z), imag(z), max(log(abs(v)), minlevel), levels)
hex(zS); axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off


% Plot angle(u_0)

nexttile; hold on
surf(real(z), imag(z), angle(v), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
hex(zS,10);
view(2); axis equal;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off;
colorbar;

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

%% Graph potential function V

V = Up*Um*const_fun(N);
[z,v] = K2X3(V,0,300,e1,e2);
% v = sqrt(v);

fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
tiledlayout(1,2,'TileSpacing','compact')

% |V|
nexttile; hold on
contourf(real(z), imag(z), abs(v), 32)
hex(zS); axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63])

% angle(V)
nexttile; hold on
surf(real(z), imag(z), angle(v), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
hex(zS,10); view(2); axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63])

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


%% Animation for log |u_K| / alpha as alpha grows

vid=VideoWriter('.\results\log_norm_asymptotics_4.mp4','MPEG-4'); open(vid)
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

alphas = 3:0.1:20;

for id=1:length(alphas)
    clf(fig);
    tl=tiledlayout(2,2,'TileSpacing','compact');
    title(tl, ['$\alpha = ' sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex')

    % Get protected state
    [~,~,V] = svds(D + alphas(id)*U + K*speye(size(D)), 1, 'smallest');
    [z,v1,~] = K2X4(V,K,300,e1,e2);

    % Plot in hexagon
    nexttile(1, [2 1]);
    hold on
    title('$\alpha^{-1} \log |u_K|$', 'Interpreter', 'latex')
    minlevel=-1;
    levels=linspace(minlevel,0.25,40);
    contourf(real(z), imag(z), max(log(abs(v1)) ./ alphas(id), minlevel), levels)
    hex(zS); axis equal; colorbar
    xlim([-0.63,0.63]), ylim([-0.63,0.63]), clim([minlevel 0.25])

    % Plot on two lines
    M = 600;
    [~,v1,~] = K2X2(V,K,M,e1,e2);
    
    nexttile;   
    hold on
    title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex')
    plot((0:M-1)'./M, log(abs(diag(v1))) ./ alphas(id))
    xlim([0 1]), ylim([-1 0.2])
    xline(1/3), xline(2/3), yline(0, 'Color', 'red')
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]

    nexttile;
    hold on
    title('restricted to line from $0$ to $\omega$', 'Interpreter', 'latex')
    plot((M-1:-1:0)'./M, log(abs(v1(1,:))) ./ alphas(id))
    xlim([0 1]), ylim([-1 1])
    xline(1/2), yline(0, 'Color', 'red')
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]

    drawnow
    frame=getframe(fig);
    writeVideo(vid,frame)
end

close(vid)


%% Animation of u_K restricted to edge of hexagon, and fourier transform

vid=VideoWriter('.\results\edge_fourier_asymptotics_3.mp4','MPEG-4'); open(vid);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

alphas = 5:0.1:20;

for id=1:length(alphas)

    clf(fig);
    tl=tiledlayout(1,2,'TileSpacing','compact');
    title(tl, ['$\alpha = ' sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex');

    [~,~,V] = svds(D + alphas(id)*U + K*speye(size(D)), 1, 'smallest');
    MM = 600; M = 3*MM;
    [~,v1,~] = K2X2(V,K,M,e1,e2);
    v0 = diag(v1); v0 = v0(MM: 2*MM-1);
    v0 = v0 / norm(v0); v0 = v0 / sign(v0(1));
    v0 = [v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0; 0*v0]; MM = 13*MM;

    % fourier transform
    V0 = fft(v0);
    V0 = V0 * sqrt(alphas(id)/13);
    V0 = circshift(V0, round(MM/2));

    nexttile;   
    hold on
    title('protected state on line from $-\frac{\sqrt{3}i}{3}$ to $-\frac{2\sqrt{3}i}{3}$', 'Interpreter', 'latex');
    plot((0:(MM/13)-1)' ./(MM/13) .*(2/sqrt(3)), v0(1:MM/13));
    xlim([0 2/sqrt(3)]);
    ylim([-3 3]);
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 1 1]; % [width, height, depth]


    nexttile;
    hold on
    title('scaled fourier transform', 'Interpreter', 'latex');
    xlabel('$\xi / \alpha$', 'Interpreter', 'latex');
    plot(((0:MM-1)' - round(MM/2)) .* (2*pi*sqrt(3)/(13*2*alphas(id))), abs(V0));
    xlim([-10 10]);
    ylim([-10 20]);
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 1 1]; % [width, height, depth]

    drawnow;
    frame=getframe(fig);
    writeVideo(vid,frame);
end

close(vid);

%% Animation for log |u_K| / alpha as potential varies

vid=VideoWriter('.\results\log_norm_varying_potential_3.mp4','MPEG-4'); open(vid);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

Up2 = sym_potential(N,1,-1); Um2 = Up2.';
U2 = [0*Up, Up2; Um2, 0*Up];

alphas = @(id) 10;
theta0 = 0:0.002:1;
theta = theta0 * pi;

for id=1:length(theta)
    
    clf(fig)
    tl=tiledlayout(3,2,'TileSpacing','compact');
    title(tl, ['$U = \cos\theta f_1 + \sin\theta f_{-2},\,\,\theta = ' ...
        sprintf('%.3f', theta0(id)) '\pi$'], 'Interpreter', 'latex')
    
    % get protected state
    [~,~,V] = svds(D + alphas(id)*(cos(theta(id))*U + sin(theta(id))*U2) + K*speye(size(D)), 1, 'smallest');
    [z,v1,~] = K2X4(V,K,300,e1,e2);

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
    M = 600;
    [~,v1,~] = K2X2(V,K,M,e1,e2);
    
    nexttile;   
    hold on
    title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex')
    plot((0:M-1)'./M, log(abs(diag(v1))) ./ alphas(id))
    xlim([0 1]), ylim([-1 0.2])
    xline(1/3), xline(2/3), yline(0, 'Color', 'red')
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]
    
    nexttile;
    hold on
    title('restricted to line from $0$ to $\omega$', 'Interpreter', 'latex')
    plot((M-1:-1:0)'./M, log(abs(v1(1,:))) ./ alphas(id))
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

%% Graph log bands / alpha    as alpha grows

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

%% Animation overlaying scalar and non-scalar log |u| / alpha  as alpha grows

vid=VideoWriter('.\results\log_norm_asymptotics_overlay2_1.mp4','MPEG-4'); open(vid);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

alphas = 10:0.1:20;

for id=1:length(alphas)
    clf(fig);
    tl=tiledlayout(2,2,'TileSpacing','compact');
    title(tl, ['$\alpha = ' sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex');

    [~,~,V] = svds(D + alphas(id)*U + K*speye(size(D)), 1, 'smallest');
    V  = V(1:length(V)/2); V = V ./norm(V);
    [z,v1] = K2X3(V,0,300,e1,e2);

    [~,~,V_scalar] = svds((4*Dbar(N,0,f1,f2)^2 -alphas(id)^2 * Up*Um - K^2), 1, 'smallest');
    [~,u0] = K2X3(V_scalar,0,300,e1,e2);


    nexttile(1);
    hold on
    title('$\alpha^{-1} \log |u_{K,1}|$, where $D(\alpha)u_K = 0$', 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.25,40);
    contourf(real(z), imag(z), max(log(abs(v1)) ./ alphas(id), minlevel), levels);
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]); 

    nexttile(3);
    hold on
    title('$\alpha^{-1} \log |v_0|$, where $v_0$ is state for smallest singular value of $Q(\alpha,0)-K^2$', 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.25,40);
    contourf(real(z), imag(z), max(log(abs(u0)) ./ alphas(id), minlevel), levels);
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]); 


    M = 600;
    [~,v1] = K2X(V,K,M,e1,e2);
    [~,u0] = K2X(V_scalar,0,M,e1,e2);
    
    nexttile
    hold on
    title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex');
    plot((0:M-1)'./M, log(abs(diag(v1))) ./ alphas(id), 'Color', 'blue');
    plot((0:M-1)'./M, log(abs(diag(u0))) ./ alphas(id), 'Color', 'red');
    xlim([0 1])
    ylim([-1 0.2])
    xline(1/3); xline(2/3); yline(0, 'Color', 'black')
    legend('nonscalar', 'scalar')
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]
    

    nexttile
    hold on
    title('restricted to line from $0$ to $\omega$', 'Interpreter', 'latex');
    plot((M-1:-1:0)'./M, log(abs(v1(1,:))) ./ alphas(id), 'Color', 'blue');
    plot((M-1:-1:0)'./M, log(abs(u0(1,:))) ./ alphas(id), 'Color', 'red');
    xlim([0 1])
    ylim([-1 1])
    xline(1/2); yline(0, 'Color', 'black');
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]

    drawnow;
    frame=getframe(fig);
    writeVideo(vid,frame);
end

close(vid);

%% Finite element for a hexagon dirichlet problem
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

M=1000;
[~,v] = K2X(V,0,M,e1,e2);
[~,dzv] = K2X(dzV,0,M,e1,e2);
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


ffun = @(z,etc) interpolate(e1,e2,f,M,z.x,z.y);
specifyCoefficients(model, "m",0,"d",0,"c",1,"a",0,"f",ffun);

generateMesh(model,"Hmax",0.01);
result=solvepde(model);

figure
hold on
pdeplot(result.Mesh, XYData=result.NodalSolution, Contour="on", ColorMap="default", Levels=18); 
axis equal
title("$\Phi$", 'Interpreter', 'latex')
% saveas(gcf,'./Phi_solution.png')


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

[z,v1,v2] = K2X4(V2,K,200,e1,e2);


nexttile; hold on
contourf(real(z), imag(z), max(log(abs(v1)) ./ 17,minlevel), levels);
hex(zS); axis equal; colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]);


end

%% Log eigenfunction / alpha   smoothened, then take laplacian

alpha = 20;
[~,~,V] = svds(D + alpha*U + K*speye(size(D)), 1, 'smallest');
V = V(1:length(V)/2); V = V / norm(V);
[z,v] = K2X(V,K,300,e1,e2);


rad = 8;
chi = ones(rad,rad); chi(300,300) = 0; %chi=circshift(chi, [-rad/2, -rad/2]);
v_conv =  ifft2(fft2(abs(v)) .* fft2(chi)) ./ rad^2;
v_conv = log(abs(v_conv)) ./ alpha;
% v_conv = max(log(abs(v)) ./ alpha, -0.5);

V_conv =X2K(v_conv,0,N,e1,e2);
[~,Lap_v_conv]=K2X(Dbar(N,0,f1,f2) * Dbar(N,0,f1,f2)' * V_conv, 0, 300, e1,e2);



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

%% (FIX) Conjugate by multiplication with e^(phi/h)

M = 300;
figure
tiledlayout(1,2)

% Sine wave in one direction (real function)
Phi = (fourier_shift(N,1,0) + fourier_shift(N,-1,0))*const_fun(N) .* 0.1;
[z, phi] = K2X(Phi, 0, M, e1, e2);
phi = real(phi);

alpha = 10;
ephi = exp(phi .* alpha);


% protected state times e^(phi/h)
alpha = 19;
[~,~,U0] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2 * Up*Um, 1, "smallest");
[~,u0] = K2X(U0, 0, M, e1, e2);

nexttile
hold on
minlevel=-1.5;
levels=linspace(minlevel,1,80);
contourf(real(z), imag(z), max(phi+ (log(abs(u0)) ./ alpha), minlevel), levels)
axis equal
colorbar


% {q, bar q}
Lap_Phi = -Dbar(N,0,f1,f2) * Dbar(N,0,f1,f2)' * Phi;
[~, lap_phi] = K2X(Lap_Phi, 0, M, e1, e2);

V = const_fun(N);
V = Up*Um*V;

dzV = 1i* Dbar(N,0,f1,f2)' * V;

[~,v] = K2X(V,0,M,e1,e2);
[~,dzv] = K2X(dzV,0,M,e1,e2);
bracket = imag( sqrt(conj(v)) .* dzv );
bracket = abs(bracket)-8 .* abs(v) .* real(lap_phi);

nexttile
hold on
contourf(real(z), imag(z), bracket, 30)
axis equal
colorbar


%% Animation for log |u|/alpha and log |2hDbar u / sqrt(V)|/alpha   as alpha grows

vid=VideoWriter('.\results\log_norm_scalar_overlay_4.mp4','MPEG-4'); open(vid);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

V0 = Up * Um * const_fun(N);

alphas = 10:0.1:20;

for id=1:length(alphas)
    clf(fig);
    tl=tiledlayout(2,2,'TileSpacing','compact');
    title(tl, ['$\alpha = ' sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex');

    [~,~,V] = svds((4*Dbar(N,0,f1,f2)^2 -alphas(id)^2 * Up*Um), 1, 'smallest');
    hdbarV = 2*Dbar(N,0,f1,f2) * V ./ alphas(id);
    [z,v] = K2X3(V,0,300,e1,e2);
    [~,dv] = K2X3(hdbarV,0,300,e1,e2);
    [~,v0] = K2X3(V0,0,300,e1,e2);
    dv = dv ./ sqrt(v0);

    nexttile(1);
    hold on
    title('$h\log |u_{0}|$, where $P(\alpha)u_0 = 0$', 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.5,40);
    contourf(real(z), imag(z), max(log(abs(v)) ./ alphas(id), minlevel), levels);
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]); 

    nexttile(3);
    hold on
    title('$h\log |2hD_{\bar z} u_{0} / \sqrt{V}|$', 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.5,40);
    contourf(real(z), imag(z), max(log(abs(dv)) ./ alphas(id), minlevel), levels);
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]); 


    M = 600;
    [~,v] = K2X(V,K,M,e1,e2);
    [~,dv] = K2X(hdbarV,0,M,e1,e2);
    [~,v0] = K2X(V0,0,M,e1,e2);
    dv = dv ./ sqrt(v0);
    dv(201,201) = dv(202,202); dv(401,401)=dv(402,402);
    
    nexttile
    hold on
    title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex');
    plot((0:M-1)'./M, log(abs(diag(v))) ./ alphas(id), 'Color', 'blue');
    plot((0:M-1)'./M, log(abs(diag(dv))) ./ alphas(id), 'Color', 'red');
    xlim([0 1])
    ylim([-1 0.2])
    xline(1/3); xline(2/3); yline(0, 'Color', 'black')
    legend('h log u', 'h log (2hDbar u / sqrt V)')
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]
    

    nexttile
    hold on
    title('restricted to line from $0$ to $\omega$', 'Interpreter', 'latex');
    plot((M-1:-1:0)'./M, log(abs(v(1,:))) ./ alphas(id), 'Color', 'blue');
    plot((M-1:-1:0)'./M, log(abs(dv(1,:))) ./ alphas(id), 'Color', 'red');
    xlim([0 1])
    ylim([-1 1])
    xline(1/2); yline(0, 'Color', 'black');
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]

    drawnow;
    frame=getframe(fig);
    writeVideo(vid,frame);
end

close(vid);


%% Does restriction of u_0 to edge of hexagon, approximately solve semiclassical ode?

alpha = 20;

[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2 * Up * Um, 1, 'smallest');
[~,v] = K2X(V,0,300,e1,e2);
phase = exp(1i*angle(v(1,1)));
V = V ./ phase;

Pot = Up*Um*const_fun(N);
DxV = Dx(N,0,f1,f2) * V;
DyV = Dy(N,0,f1,f2) * V;

DxxV = Dx(N,0,f1,f2) * DxV;
DxyV = Dy(N,0,f1,f2) * DxV;
DyyV = Dy(N,0,f1,f2) * DyV;


M = 500;
[~,v] = K2X(V,0,3*M,e1,e2);
[~,dxv] = K2X(DxV,0,3*M,e1,e2);
[~,dyv] = K2X(DyV,0,3*M,e1,e2);
[~,dxxv] = K2X(DxxV,0,3*M,e1,e2);
[~,dxyv] = K2X(DxyV,0,3*M,e1,e2);
[~,dyyv] = K2X(DyyV,0,3*M,e1,e2);
[~,pot] = K2X(Pot,0,3*M,e1,e2);


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
A = A(:,M+1+mm:2*M-mm);

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
plot((mm:M-1-mm)'./M, A(4,:) ./ A(1,:))
plot((mm:M-1-mm)'./M, A(5,:) ./ A(1,:))
plot((mm:M-1-mm)'./M, A(6,:) ./ A(1,:))
plot((mm:M-1-mm)'./M, real(A(7,:)))
ylim([-100, 100])
legend(["$(4h^2 \partial_x^2 u) / u$", "$(4h^2 \partial_y^2 u) / u$", "$\Im (4h^2 \partial_x\partial_y u) / u$", "V"], 'Interpreter', 'Latex')


% disp( (imag(dyv)' * dxv) / (norm(dxv) * norm(dyv)) )
% disp( -norm(dxv) / norm(dyv) )



%% (related to previous) calculate c (coefficient in ode) for some values of alpha

alphas=5:0.1:20;
xx_yy = alphas * 0;
xx_xy = alphas * 0;
x_y = alphas * 0;

for id=1:length(alphas)
alpha = alphas(id);



[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2 * Up * Um, 1, 'smallest');
[~,v] = K2X(V,0,300,e1,e2);
phase = exp(1i*angle(v(1,1)));
V = V ./ phase;

Pot = Up*Um*const_fun(N);
DxV = Dx(N,0,f1,f2) * V;
DyV = Dy(N,0,f1,f2) * V;

DxxV = Dx(N,0,f1,f2) * DxV;
DxyV = Dy(N,0,f1,f2) * DxV;
DyyV = Dy(N,0,f1,f2) * DyV;


M = 500;
[~,v] = K2X(V,0,3*M,e1,e2);
[~,dxv] = K2X(DxV,0,3*M,e1,e2);
[~,dyv] = K2X(DyV,0,3*M,e1,e2);
[~,dxxv] = K2X(DxxV,0,3*M,e1,e2);
[~,dxyv] = K2X(DxyV,0,3*M,e1,e2);
[~,dyyv] = K2X(DyyV,0,3*M,e1,e2);
[~,pot] = K2X(Pot,0,3*M,e1,e2);


v = diag(v); dxv = diag(dxv).* (2 / alpha); dyv = diag(dyv).* (2 / alpha);
dxxv = diag(dxxv).* (2 / alpha)^2; dyyv = diag(dyyv).* (2 / alpha)^2; 
dxyv = diag(dxyv).* (2 / alpha)^2; pot = diag(pot);


A = [real(v), real(dxv), imag(dyv), real(dxxv), real(dyyv), imag(dxyv), pot].';
names = ["u", "dxu", "dyu", "dxxu", "dyyu", "dxyu", "potential V"];

mm = 20;
A = A(:,M+1+mm:2*M-mm);



ind = [2 3];
[c,~,~] = svd(A(ind,:), "econ");
x_y(id) = c(3)/c(4);

ind = [4 5];
[c,~,~] = svd(A(ind,:), "econ");
xx_yy(id) = c(3)/c(4);

ind = [4 6];
[c,~,~] = svd(A(ind,:), "econ");
xx_xy(id) = c(3)/c(4);

end

figure
hold on
plot(alphas, xx_yy)
plot(alphas, xx_xy)
plot(alphas, x_y)
legend(["$\sqrt{\partial_{yy} u / \partial_{xx} u}$",...
        "$(Im \partial_{xy} u) / \partial_{xx} u$", "$\partial_y u / Im(\partial_x u)$"],...
       'Interpreter', 'Latex')

%% Animation of potential changing (between f_1 and f_-2) and corresponding protected state

vid=VideoWriter('.\results\varying_potential_2_1.mp4','MPEG-4'); open(vid);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

Up2 = sym_potential(N,1,-1); Um2 = Up2.';

alpha = @(id) 10;
theta = (0:0.01:1) * pi;

for id=1:length(theta)
    thet = theta(id);
    
    clf(fig);
    tl=tiledlayout(1,3,'TileSpacing','compact');
    V = (cos(thet)^2 * Up - sin(thet)^2 * Up2)* (cos(thet)^2 * Um - sin(thet)^2 * Um2);
    [~,~,U] = svds((4*Dbar(N,0,f1,f2)^2 -alpha(id)^2 * V), 1, 'smallest');
    [z,u0] = K2X3(U,0,300,e1,e2);

    [~,v] = K2X3( V * const_fun(N), 0, 300, e1,e2 );
    v = v ./ exp(1i*angle(v(150,150)));

    title(tl, ['$V(z) = U(z)U(-z), U(z) = \cos^2(\theta)f_1 - \sin^2(\theta) f_{-2}, \theta = '...
               sprintf('%.2f', id / length(theta)) '\pi$'], 'Interpreter', "latex");

    nexttile    % plot highest points? nah for now, just wanna see
    hold on
    title("$h \log |u_0|$", 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.25,40);
    contourf(real(z), imag(z), max(log(abs(u0)) ./ alpha(id), minlevel), levels);
    hex(zS)
    axis equal
    colorbar
    xlim([-0.63,0.63]), ylim([-0.63,0.63])
    clim([minlevel 0.25]); 
    
    nexttile
    hold on
    title("$|V|$", 'Interpreter', 'Latex')
    levels=linspace(0,150,70);
    contourf(real(z), imag(z), min(abs(v), 150), levels)
    hex(zS,10)
    axis equal
    xlim([-0.63,0.63])
    ylim([-0.63,0.63])
    clim([0 60])
    colorbar

    nexttile
    hold on
    title("$angle(V)$", 'Interpreter', 'latex')
    surf(real(z), imag(z), angle(v), 'EdgeColor', 'none')
    clim([-pi,pi])
    colormap(gca,wheelmap)
    hex(zS,10)
    view(2)
    axis equal
    xlim([-0.63,0.63])
    ylim([-0.63,0.63])
    colorbar



    drawnow;
    frame=getframe(fig);
    writeVideo(vid,frame);

end

close(vid);


%% Comparison with Airy function in w coordinates

% prologue -------------


vid=VideoWriter('.\results\airy_holo_6.mp4','MPEG-4'); open(vid);
fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);


alphas = 0.1:0.1:19.5;

for id=1:length(alphas)
alpha = alphas(id);
% alpha = 17;

% ----------------------
disp(id)
clf(fig);
tl=tiledlayout(2,3,'TileSpacing','compact');
title(tl, ['$\alpha = ' sprintf('%.1f', real(alpha / 1i)) 'i$'], 'Interpreter', 'latex')

A = -1i*(8 * pi^3 * (alpha)^2 / 3)^(-1/3);

% Range of w
[x,y] = meshgrid(-1:0.005:1, -1:0.005:1);
w = x + 1i*y; w= 5*w;
z = zS + A * w;

% Get protected state
M=4000;
[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - (alpha)^2*Up*Um, 1, 'smallest');
[~,v] = K2X(V,0,M,e1,e2);


apply = @(f) arrayfun(@(z) interpolate(e1,e2,f,M,real(z),imag(z)), z);
phase = exp(1i*angle(interpolate(e1,e2,v,M,real(zS), imag(zS) ) ) );


v = apply(v);


[~,dv] = K2X(1i*Dbar(N,0,f1,f2)*V, 0,M,e1,e2);
dv = apply(A*dv);

dv = dv ./ phase;v = v ./ phase; 



% Airy function and derivatives
c = (2/3)*(16/3)^(1/3);
% c=0;
ww = conj(w) - c * alpha^(-2/3) * w.^2;
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
% wronsF = F0 .* dF2 - F2 .* dF0; % 0.8270
wronsF = F0 * 0 -3*sqrt(3)/(2*pi);
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
title("$|\alpha|^{-1}\log |u_1(w)|$", 'Interpreter', 'latex')
contourf(real(w), imag(w), max(log(abs(v)) ./ abs(alpha), -1.2), linspace(-1.2, 0.2, 50))
contourf(real(w), imag(w), abs(v), 32)
axis equal
clim([-1.2, 0])
colorbar

nexttile(4)
hold on
title("$angle(v(w))$", 'Interpreter', 'latex')
surf(real(w), imag(w), angle(v), 'EdgeColor', 'none')
clim([-pi,pi]), colormap(gca,wheelmap)
view(2)
axis equal
xlim([-5, 5]), ylim([-5 5])
colorbar

nexttile(2)
hold on
title("$|\alpha|^{-1}\log |g_0(w^3)|$", 'Interpreter', 'latex')
contourf(real(w), imag(w), min(max(log(abs(g0)) ./ abs(alpha), -1.2), 1),  linspace(-1.2, 1, 50))
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
title("$|\alpha|^{-1}\log |g_2(w^3)|$", 'Interpreter', 'latex')
contourf(real(w), imag(w), min(max(log(abs(g2)) ./ abs(alpha), -1.2), 1),  linspace(-1.2, 1, 50))
axis equal
clim([-1.2 1])
colorbar

nexttile(6)
hold on
title("$angle(g_2(w^3))$", 'Interpreter', 'latex')
surf(real(w), imag(w), angle(g2), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
view(2)
axis equal
xlim([-5, 5]), ylim([-5 5])
colorbar

% epilogue -------------
drawnow
frame=getframe(fig);
writeVideo(vid,frame)
end

close(vid)

% ----------------------


%% for some set of alpha, compute g_0(0) and g_0(2) from wronskian with Airy function


% prologue -------------

% [xx,yy] = meshgrid(0.1:0.15:15, 0.1:0.15:13);
% [xx,yy] = meshgrid(0.1:0.2:5, 0.1:0.5:5);
% alphas = xx + 1i * yy; alphas = alphas(:);
alphas = 5:0.03:18;
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
[~,v] = K2X(V,0,M,e1,e2);
[~,dv] = K2X(1i*Dbar(N,0,f1,f2)*V, 0,M,e1,e2);
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

% plot(alphas, log(g00))
% plot(alphas, log(g20))
% 
% Ak = Inv(4*Dbar(N,K,f1,f2)^2) *  Up * Um;
% Alphas = 1./sqrt(eigs(Ak, 500));
% Alphas = Alphas(abs(imag(Alphas)) < 0.001 & 0 < real(Alphas) & real(Alphas) < 22);
% for id=1:length(Alphas)
%     xline(real(Alphas(id)));
% end
% 
% legend(["$\log(g_0(0))$", "$\log(g_2(0))$", "magic $\alpha$"], 'Interpreter', 'Latex')
% saveas(fig,'./results/g_00_2.png')

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



%% Logarithmic derivative of protected state u

% V(x)

alpha = 10;
[~,~,V] = svds(D + alpha*U + K*speye(size(D)), 1, 'smallest');
V = V(1:length(V)/2); V = V / norm(V);
[z,v] = K2X(V,K,300,e1,e2);


rad = 20;
chi = ones(rad,rad); chi(300,300) = 0; %chi=circshift(chi, [-rad/2, -rad/2]);
v_conv =  ifft2(fft2(abs(v)) .* fft2(chi)) ./ rad^2;
v_conv = v_conv .* exp(1i*angle(v));

V_conv =X2K(v_conv,0,N,e1,e2);
[~,Dbar_v_conv]=K2X(Dbar(N,0,f1,f2) * V_conv, 0, 300, e1,e2);

vv = (Dbar_v_conv ./ v_conv) .^ 2;

figure
til=tiledlayout(1,4,'TileSpacing','compact');


% nexttile;
% hold on
% title('$\chi$ where $\int \chi = 1$', 'Interpreter', 'latex');
% contourf(real(z), imag(z), -chi ./ 2, levels);
% axis equal; clim([minlevel 0.25]); 


nexttile;
hold on
contourf(real(z), imag(z), max(min(real(vv), 40000), -40000), 32);
axis equal
colorbar

nexttile;
hold on
contourf(real(z), imag(z), max(log(abs(v)), -20), 32);
axis equal; colorbar;

%% Phase function branch cut WKB (from AM's Phys Rev Letter paper on WKB)

M = 300;

V = Up*Um*const_fun(N);

[~,v] = K2X(V,0,M,e1,e2);
v = sqrt(-v);



V2 = X2K(v,0,N,e1,e2);

Db = Dbar(N,0,f1,f2);
ww = spdiags([ 0*(1:N) 1 0*(1:N)]', 0, 2*N+1, 2*N+1);
Db = Db + kron(ww, ww);

V2 = Inv(Db) * V2;

val = V2((N)*(2*N+1) + N+1) *sqrt(sqrt(3)/2);
V2((N)*(2*N+1) + N+1) = 0;

[z,v] = K2X3(V2,0,M,e1,e2,'nophase');
% [~,v1] = K2X3(ROT(N,0)*V2,0,M,e1,e2);
% [~,v2] = K2X3(ROT(N,0)^2*V2,0,M,e1,e2);


V = Up*Um*const_fun(N);
[~,v0] = K2X3(V,0,M,e1,e2,'nophase');
v0 = v0 .^ (-0.25);
v0(imag(z) > 0) = -1i*v0(imag(z)>0);

v = v + 0.5 * (val*z'+val'*z); % (no diff? lol)


% vid=VideoWriter('.\results\physics_wkb_3.mp4','MPEG-4'); open(vid);
% fig=figure;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% 
% for alpha=0:0.1:10

    % clf(fig);
    tl=tiledlayout(1,2,'TileSpacing','compact');
    % title(tl, ['$\alpha = ' sprintf('%.1f', alpha) '$'], 'Interpreter', 'latex');


% vv = (exp(1i*v*alpha)  + exp(-1i*v*alpha)) .* v0;

nexttile
hold on
title("real part")
contourf(real(z), imag(z), real(v), 30)
axis equal; colorbar
% xlim([-0.63,0.63]), ylim([-0.63,0.63])

nexttile
hold on
title("imaginary part")
contourf(real(z), imag(z), imag(v), 30)
axis equal; colorbar
% xlim([-0.63,0.63]), ylim([-0.63,0.63])

% nexttile
% surf(real(z), imag(z), angle(vv), 'EdgeColor', 'none');
% clim([-pi,pi]); colormap(gca,wheelmap);view(2);
% hex(zS,10); axis equal; colorbar
% xlim([-0.63,0.63]), ylim([-0.63,0.63])
% 
% nexttile
% hold on
% contourf(real(z), imag(z), max(log(abs(vv)) ./ alpha, -1), -1:0.05:2)
% hex(zS); axis equal; colorbar
% xlim([-0.63,0.63]), ylim([-0.63,0.63])


%     drawnow;
%     frame=getframe(fig);
%     writeVideo(vid,frame);
% end
% 
% close(vid);


%% Log of protected state 2x2 matrix (non-scalar)

alpha = 17.6;
[~,s1,V1] = svds(D + K*speye(size(D)) + alpha*U, 1, 'smallest');
[~,s2,V2] = svds(D - K*speye(size(D)) + alpha*U, 1, 'smallest');

M = 600;
[z,v11,v12] = K2X2(V1,K,M,e1,e2);
[~,v21,v22] = K2X2(V2,-K,M,e1,e2);


% figure
% hold on
% title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex')
% plot((0:M-1)'./M, log(abs(diag(v11))) ./ alpha, 'color', 'black')
% plot((0:M-1)'./M, log(abs(diag(v12))) ./ alpha, 'color', 'red')
% plot((0:M-1)'./M, log(abs(diag(v21))) ./ alpha, 'color', 'blue')
% plot((0:M-1)'./M, log(abs(diag(v22))) ./ alpha, 'color', 'green')

% xlim([0 1]), ylim([-1 0.2])
% xline(1/3), xline(2/3), yline(0, 'Color', 'red')
% ax = gca; % Get current axes
% ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]

e = 0*z; a = 0*z; b = 0*z; c = 0*z; E1 = 0*z; E2 = 0*z;
for ii=1:size(z,1)
    for jj=1:size(z,2)
        v = logm([v11(ii,jj), v21(ii,jj); v12(ii,jj), v22(ii,jj)]) / alpha;
        e(ii,jj) = v(1,1) + v(2,2);
        a(ii,jj) = v(1,1) - v(2,2);
        b(ii,jj) = v(2,1) + v(1,2);
        c(ii,jj) = v(2,1) - v(1,2);

        EE = eigs([v11(ii,jj), v21(ii,jj); v12(ii,jj), v22(ii,jj)]);
        E1(ii,jj) = EE(1); E2(ii,jj) = EE(2);
    end
end

%

figure
hold on
contourf(real(z), imag(z), max(real(log(E1)), -20), 30);
hex(zS); axis equal; colorbar
xlim([-0.63,0.63]), ylim([-0.63,0.63])

%% For non-symmetric potential, magic angles changing for k in a loop

Up = sym_potential(N,0,0) + 0.001*fourier_shift(N,0,0);
Um = Up.';
%for f_{-2}, use   Up2 = sym_potential(N,1,-1); Um2 = Up2.';

Dk = @(k) D + k*speye(size(D));
U = [0*Up, Up; Um, 0*Up];


ks = [-0.12:0.000701:0.12, 0.15];
% ks = exp(1i*2*pi*(0:0.01:1));
ks = 0.1*ks - K;
colors = spring(length(ks));

figure
hold on

for id=1:length(ks)
    k=ks(id);
    Alphas = 1./(eigs(Inv(Dk(ks(id)))*U, 200));
    plot(real(Alphas), imag(Alphas), '.', 'color', colors(id,:))
    drawnow
    disp(id)
end

Alphas = 1./(eigs(Inv(Dk(0))*U, 200));
plot(real(Alphas), imag(Alphas), 'x', 'Color', 'black')


%% P(alpha) representation using position space basis

V = sym_potential(10,0,0) * sym_potential(10,0,0).' * const_fun(10);
M=300;
[z,v_pot] = K2X(V, 0, M, e1, e2);
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

%% D(alpha) using position space basis

k = -K;

M = 200;
[z,up0] = K2X(sym_potential(10,0,0)*const_fun(10), -2*K, M, e1, e2);
[~,um0] = K2X(sym_potential(10,0,0).'*const_fun(10), 2*K, M, e1, e2);

Db_maker = @(k) FE_Dbar2(M, f1, f2, e1, e2, ...
    @(z) 0*z+exp(1i*(k*e1'+k'*e1)/2), @(z) 0*z+exp(1i*(k*e2'+k'*e2)/2));
Db_1 = Db_maker(k-K);
Db_2 = Db_maker(k+K);

alpha = 5;
D = [Db_1, alpha*Mult(up0(:),M); alpha*Mult(um0(:),M), Db_2];
[~,s,v] = svds(D, 1, 'smallest');
disp(s)
v1 = v(1:length(v)/2); v2 = v(1+length(v)/2:end);
v1 = reshape(v1, [M M]); v2 = reshape(v2, [M M]);

z = z((1:2:M), (1:2:M));
v1 = v1((1:2:M), (1:2:M));
v2 = v2((1:2:M), (1:2:M));
v1 = v1 / norm(v1); v2 = v2 / norm(v2);


figure

tiledlayout(1,2,'TileSpacing','compact')
minlevel=-16;
levels=linspace(minlevel,1,28);

% Plot log(abs(u_0))

nexttile
hold on
contourf(real(z), imag(z), max(log(abs(v1)), minlevel), levels)
hex(zS); axis equal; colorbar


% Plot angle(u_0)

nexttile
hold on
contourf(real(z), imag(z), max(log(abs(v2)), minlevel), levels)
hex(zS); axis equal; colorbar


%% Section of -1 degree line bundle at magic angle




%% FUNCTIONS -----------------------------------------------------

% Below are functions related to the space L^2_k(C/(Z e1 + Z e2); C) 
% (or after fourier transform,  L^2(Z f1 + Z f2 + k; C) ).

% We represent our operators (Dbar, V, etc) in the fourier basis.
% To go from fourier to position space, use (K2X, K2X2, K2X3, K2X4).
% To go from position to fourier space, use X2K.
% We also have some symmetry operators.

% About fourier basis: a vector is a length (2N+1)^2 column vector V, represent 
% v(z) = \sum_{a,b\in{-N,...,N}} V((2N+1)*(b+N)+a+N+1) e^{i<z, a f1 + b f2+ k>} / sqrt(area_of_(e1,e2))

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



% Fourier transforms -----------

% Get position-space vector v corresponding to k-space vector V:
% v(z) = sum V(n1,n2) e^(i<n1 f1 + n2 f2 + k, z>)
%
% For given e1,e2,M,   
% z  will be MxM matrix  with values   z(a,b) = e1*(a-1)/M + e2*(b-1)/M
function [z,v] = K2X(V,k,M,e1,e2)
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
function V=X2K(v,k,N,e1,e2)
    M = size(v, 1);
    [y2,y1]=meshgrid(0:M-1, 0:M-1);
    z=(e1*y1+e2*y2)/M;

    v = v ./ exp(0.5i*(k'*z+k*conj(z)));
    
    V = fft2(v) * sqrt(abs(imag(e1' * e2))) / M^2;
    
    V = circshift(V, [N,N]);
    V = V(1:2*N+1, 1:2*N+1);
    V = V(:);
end


% Same as K2X but do for two layers (D(alpha) is 2x2 matrix)
function [z,v1,v2] = K2X2(V,k,M,e1,e2)
    K = 4*pi/3;

    [z,v1] = K2X(V(1:length(V)/2),k-K,M,e1,e2);
    [~,v2] = K2X(V(1+length(V)/2:end),k+K,M,e1,e2);

    phase=exp(1i * angle(v1(1)) );
    v1=v1 ./ phase; v2=v2./phase;
end

% Same as K2X, K2X2 but repeat over some copies of fundamental domain

function [z,v] = K2X3(V,k,M,e1,e2,nophase)
    [z0,v0]=K2X(V,k,M,e1,e2);

    z = [z0-e1-e2, z0-e1; z0-e2, z0];
    
    p2 = exp(-0.5i*(k'*e2+k*e2')); p1 = exp(-0.5i*(k'*e1+k*e1'));
    v = [v0 *p1 *p2, v0 *p1; v0 *p2, v0];
    if ~exist('nophase', 'var')
        phase = exp(1i * angle(v0(1)));
        v = v ./  phase;
    end
end

function [z,v1,v2] = K2X4(V,k,M,e1,e2,nophase)
    K = 4*pi/3;

    [z,v1] = K2X3(V(1:length(V)/2),k-K,M,e1,e2);
    [~,v2] = K2X3(V(1+length(V)/2:end),k+K,M,e1,e2);

    if ~exist('nophase', 'var')
        phase=exp(1i*angle(v1(length(v1)/2)));
        v1=v1 ./ phase; v2=v2./phase;
    end
end



% Symmetries of TBG ---------------

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


% General operator permuting fourier basis

% (FV)(F1(m1,m2),F2(m1,m2)) = V(m1,m2)

function F = PERMK(N,F1,F2)
    indx =  @(m1,m2) (2*N+1)*(m2+N) + m1+N + 1;
    [m1,m2]=meshgrid(-N:N,-N:N); m1=m1(:); m2=m2(:);

    fm1 = F1(m1,m2);
    fm2 = F2(m1,m2);

    mask = fm1 <= N & fm1 >= -N & fm2 <= N & fm2 >= -N;
    F = sparse(indx(fm1(mask),fm2(mask)), indx(m1(mask),m2(mask)), ones(length(find(mask)),1), (2*N+1)^2, (2*N+1)^2);
end




% Some position space operators --------------

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
% H1 and H2 are holomorphic "multipliers" defining a holomorphic line bundle:
%  - functions u satisfying  u(z + e1) = H1(z) u(z),  u(z + e2) = H2(z) u(z),
% to be multipliers, H1, H2 must satisfy compatibility condition
function Db=FE_Dbar2(M, f1, f2, e1, e2, H1, H2)

    A = spdiags(ones(M,1), 1, M, M); B = sparse(M,M); B(M,1) = 1;
    D1_front = kron(speye(M,M), A) + kron(spdiags(H1(e2*(0:M-1)' ./ M) ,0,M,M),   B);
    D1_back = kron(speye(M,M), A') + kron(spdiags(1./H1(e1*(M-1)/M + e2*(0:M-1)' ./ M) ,0,M,M),   B');

    D2_front = kron(A, speye(M,M)) + kron(B, spdiags(H2(e1*(0:M-1)' ./ M) ,0,M,M));
    D2_back = kron(A', speye(M,M)) + kron(B', spdiags(1./H2(e2*(M-1)/M + e1*(0:M-1)' ./ M) ,0,M,M));

    D1 = (D1_front - D1_back) * M/(2i);
    D2 = (D2_front - D2_back) * M/(2i);

    Db = (f1 * D1 + f2 * D2) ./ (4*pi);
end


% Multiply pointwise by V
function Mat=Mult(V, M)
    Mat = spdiags(V, 0, M^2, M^2);
end


% misc ----------------------------

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
    set1=C(count==1,:); set2=C(count==2,:); set3=C(count>2,:);
    
    scatter(set1(:,1),set1(:,2),dotsize,'blue','filled');
    scatter(set2(:,1),set2(:,2),dotsize,'red','filled');
    scatter(set3(:,1),set3(:,2),dotsize,'black','filled');

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
