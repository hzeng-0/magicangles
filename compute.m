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



%% Retrieve magic angles
save_angles=load('.\angles\1_-2_0.0002_circle.mat', 'save_angles').('save_angles');


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

[~,~,V] = svds(4*Dbar(N,0,f1,f2)^2 - (2+6i)^2*Up*Um, 1, 'smallest');
[z,v1] = K2X3(V,0,200,e1,e2);

fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
til=tiledlayout(1,2,'TileSpacing','compact');
minlevel=-16;
levels=linspace(minlevel,1,28);


nexttile; hold on;
contourf(real(z), imag(z), max(log(abs((v1))),minlevel), levels);
hex(zS); axis equal; colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off;

% nexttile; hold on;
% contourf(real(z), imag(z), max(log(abs((v2))),minlevel), levels);
% hex(zS); axis equal; colorbar;
% xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off;

nexttile; hold on;
surf(real(z), imag(z), angle(v1), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
hex(zS,10);
view(2); axis equal;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off;
colorbar;


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

% inp = @(x,y) (x'*y+x*y')/2;
% U_fun = @(x) exp(1i*inp(x,K)) + w*exp(1i*inp(x,w*K)) +  w^2*exp(1i*inp(x,w^2*K));
% 
% V_fun = @(x) U_fun(x) * U_fun(-x);
% 
% V = @(x) abs(exp(1i*x*K) - exp(-1i*x*K/2));
% xx = 0:0.001:1;
% VV = arrayfun(V, xx);
% disp(sum(VV) / length(VV));

 V = const_fun(N);
 V = Up*Um*V;
[z,v] = K2X3(V,K,300,e1,e2);

fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
til=tiledlayout(1,2,'TileSpacing','compact');


nexttile; hold on;
contourf(real(z), imag(z), abs(v), 32);
hex(zS); axis equal; colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]);



nexttile; hold on;
surf(real(z), imag(z), angle(v), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
hex(zS,10);
view(2); axis equal;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off;
colorbar;


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

vid=VideoWriter('.\results\different_k_state_x.mp4','MPEG-4'); open(vid);
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
    sizes = '$0:' +  string(norm(V10) >0.001) + ', 1:' +  string(norm(V11) >0.001) + ', 2:' +  string(norm(V12) >0.001 + '$');
    if norm(V10)> 0.001
        V1 = V10;
    elseif norm(V11) > 0.001
        V1 = V11;
    else
        V1 = V12;
    end
    V1 = V1 / norm(V1);
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
    title(sizes, 'Interpreter', 'Latex');
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


%% Animation for log |u_K| / alpha as alpha grows


vid=VideoWriter('.\results\log_norm_asymptotics_4.mp4','MPEG-4'); open(vid);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

alphas = 3:0.1:20;

for id=1:length(alphas)
    clf(fig);
    tl=tiledlayout(2,2,'TileSpacing','compact');
    title(tl, ['$\alpha = ' sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex');

    [~,~,V] = svds(D + alphas(id)*U + K*speye(size(D)), 1, 'smallest');
    [z,v1,~] = K2X4(V,K,300,e1,e2);


    nexttile(1, [2 1]);
    hold on;
    title('$\alpha^{-1} \log |u_K|$', 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.25,40);
    contourf(real(z), imag(z), max(log(abs(v1)) ./ alphas(id), minlevel), levels);
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]); 


    M = 600;
    [~,v1,~] = K2X2(V,K,M,e1,e2);
    
    nexttile;   
    hold on;
    title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex');
    plot((0:M-1)'./M, log(abs(diag(v1))) ./ alphas(id));
    xlim([0 1]);
    ylim([-1 0.2]);
    xline(1/3); xline(2/3); yline(0, 'Color', 'red');
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]

    nexttile;
    hold on;
    title('restricted to line from $0$ to $\omega$', 'Interpreter', 'latex');
    plot((M-1:-1:0)'./M, log(abs(v1(1,:))) ./ alphas(id));
    xlim([0 1]);
    ylim([-1 1]);
    xline(1/2); yline(0, 'Color', 'red');
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]

    drawnow;
    frame=getframe(fig);
    writeVideo(vid,frame);
end

close(vid);


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
    hold on;
    title('protected state on line from $-\frac{\sqrt{3}i}{3}$ to $-\frac{2\sqrt{3}i}{3}$', 'Interpreter', 'latex');
    plot((0:(MM/13)-1)' ./(MM/13) .*(2/sqrt(3)), v0(1:MM/13));
    xlim([0 2/sqrt(3)]);
    ylim([-3 3]);
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 1 1]; % [width, height, depth]


    nexttile;
    hold on;
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
    
    clf(fig);
    tl=tiledlayout(3,2,'TileSpacing','compact');
    title(tl, ['$U = \cos\theta f_1 + \sin\theta f_{-2},\,\,\theta = ' sprintf('%.3f', theta0(id)) '\pi$'], 'Interpreter', 'latex');
    
    [~,~,V] = svds(D + alphas(id)*(cos(theta(id))*U + sin(theta(id))*U2) + K*speye(size(D)), 1, 'smallest');
    [z,v1,~] = K2X4(V,K,300,e1,e2);
    v1_max = max(abs(v1), [], 'all');
    max_inds = find(abs(abs(v1) - v1_max) < 0.0001);
    
    
    nexttile(1, [3 1]);
    hold on;
    title(['$\alpha^{-1} \log |u_K|,\quad\quad \alpha = '  sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.25,40);
    contourf(real(z), imag(z), max(log(abs(v1)) ./ alphas(id), minlevel), levels);
    hex(zS); axis equal; colorbar;
    plot(real(z(max_inds)), imag(z(max_inds)), '.', 'color', 'red', 'MarkerSize', 10);
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]);

    
    M = 600;
    [~,v1,~] = K2X2(V,K,M,e1,e2);
    
    nexttile;   
    hold on;
    title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex');
    plot((0:M-1)'./M, log(abs(diag(v1))) ./ alphas(id));
    xlim([0 1]);
    ylim([-1 0.2]);
    xline(1/3); xline(2/3); yline(0, 'Color', 'red');
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]
    
    nexttile;
    hold on;
    title('restricted to line from $0$ to $\omega$', 'Interpreter', 'latex');
    plot((M-1:-1:0)'./M, log(abs(v1(1,:))) ./ alphas(id));
    xlim([0 1]);
    ylim([-1 0.5]);
    xline(1/2); yline(0, 'Color', 'red');
    ax = gca; % Get current axes
    ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]
    


    nexttile;
    Alphas = save_angles(:, 1+10*(id-1));
    
    hold on;
    scattermult([real(Alphas), imag(Alphas)], 13);
    xlim([-0.1 11]); ylim([-8 8]); yline(0);

    drawnow;
    frame=getframe(fig);
    writeVideo(vid,frame);

end

close(vid);

%% Bracket and asymptotic (log |u_K|)/alpha ?

% lambda = 1

 V = const_fun(N);
 V = Up*Um*V;

% V = (w*fourier_shift(N,0,1) + w^2*fourier_shift(N,-1,0)+ ...
%     fourier_shift(N,-1,-1) + w*fourier_shift(N,0,-1)+ ...
%     w^2*fourier_shift(N,1,0)+fourier_shift(N,1,1)) * const_fun(N);

dzV = 1i* Dbar(N,0,f1,f2)' * V;

M=600;
[z,v] = K2X(V,0,M,e1,e2);
[~,dzv] = K2X(dzV,0,M,e1,e2);
bracket = imag( sqrt(conj(v)) .* dzv );
bracket = abs(bracket); 
bracket = bracket - max(bracket, [], 'all');



factor = 3000;
alpha=20;

% [~,~,V1] = svds(D + alpha*U + K*speye(size(D)), 1, 'smallest');
% [~,v1,~] = K2X2(V1,K,M,e1,e2);

[~,~,V1] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2*Up*Um, 1, 'smallest');
[~,v1] = K2X(V1,K,M,e1,e2);

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

tl=tiledlayout(2,2,'TileSpacing','compact');
nexttile;
hold on;
title('$\alpha^{-1} \log |u_K|$', 'Interpreter', 'latex');
minlevel=-1;
levels=linspace(minlevel,0.25,40);
contourf(real(z), imag(z), max(log(abs(v1)) ./ alpha, minlevel), levels);
hex(zS); axis equal; colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]);

nexttile;
hold on;
title('$|\{q, \overline q\}|_{q=0}$', 'Interpreter', 'latex');
contourf(real(z), imag(z), max(bracket ./ factor, minlevel), levels);
hex(zS); axis equal; colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]);


nexttile;   
hold on;
title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex');
plot((0:M-1)'./M, log(abs(diag(v1))) ./ alpha);
plot((0:M-1)'./M, diag(bracket)./factor, 'color', 'green');
xlim([0 1]);
ylim([-1.3 0.2]);
xline(1/3); xline(2/3); yline(0, 'Color', 'red');
ax = gca; % Get current axes
ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]

nexttile;
hold on;
title('restricted to line from $0$ to $\omega$', 'Interpreter', 'latex');
plot((M-1:-1:0)'./M, log(abs(v1(1,:))) ./ alpha);
plot((M-1:-1:0)'./M, bracket(1,:) ./ factor, 'color', 'green');
xlim([0 1]);
ylim([-1.5 0.5]);
xline(1/2); yline(0, 'Color', 'red');
ax = gca; % Get current axes
ax.PlotBoxAspectRatio = [1 0.5 1]; % [width, height, depth]



% doesn't seem right.

%% Plot log bands / alpha    as alpha grows

alphas=3:0.05:20;

E = zeros(20, length(alphas));

for id=1:length(alphas)

    disp(alphas(id))

    E(:,id) = log( svds(D + alphas(id)*U, 20, 'smallest') ) ./ alphas(id);

end

fig=figure; hold on;
tl=tiledlayout(2,1,'TileSpacing', 'compact');
title(tl, '$k=0$', 'Interpreter', 'latex');

nexttile; hold on;
title('$\alpha^{-1}\log E_j$ for $j=1,2,3$', 'Interpreter', 'latex');
plot(alphas, E(1:3, :));
xlabel('$\alpha$', 'Interpreter', 'latex');
ylim([-1, 0.5]); xlim([3, 20]);

nexttile; hold on;
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
    [~,v1_scalar] = K2X3(V_scalar,0,300,e1,e2);


    nexttile(1);
    hold on;
    title('$\alpha^{-1} \log |u_{K,1}|$, where $D(\alpha)u_K = 0$', 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.25,40);
    contourf(real(z), imag(z), max(log(abs(v1)) ./ alphas(id), minlevel), levels);
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]); 

    nexttile(3);
    hold on;
    title('$\alpha^{-1} \log |v_0|$, where $v_0$ is state for smallest singular value of $Q(\alpha,0)-K^2$', 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.25,40);
    contourf(real(z), imag(z), max(log(abs(v1_scalar)) ./ alphas(id), minlevel), levels);
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]); 


    M = 600;
    [~,v1] = K2X(V,K,M,e1,e2);
    [~,v1_scalar] = K2X(V_scalar,0,M,e1,e2);
    
    nexttile
    hold on
    title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex');
    plot((0:M-1)'./M, log(abs(diag(v1))) ./ alphas(id), 'Color', 'blue');
    plot((0:M-1)'./M, log(abs(diag(v1_scalar))) ./ alphas(id), 'Color', 'red');
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
    plot((M-1:-1:0)'./M, log(abs(v1_scalar(1,:))) ./ alphas(id), 'Color', 'red');
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



%% Alpha changing through complex numbers

vid=VideoWriter('.\results\test0.mp4','MPEG-4'); open(vid);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);


for id=-1.2:0.02:1.2
    clf(fig)
alpha = 5+ 2i*id;

[~,s,V] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2*Up*Um, 1, 'smallest');
[z,v1] = K2X3(V,K,200,e1,e2);



til=tiledlayout(1,2,'TileSpacing','compact');
title(til, string(alpha))
minlevel=-16;
levels=linspace(minlevel,1,28);


nexttile; hold on;
contourf(real(z), imag(z), max(log(abs((v1))),minlevel), levels);
hex(zS); axis equal; colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([-6 1]); hold off;

% nexttile; hold on;
% contourf(real(z), imag(z), max(log(abs((v2))),minlevel), levels);
% hex(zS); axis equal; colorbar;
% xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off;

nexttile; hold on;
surf(real(z), imag(z), angle(v1), 'EdgeColor', 'none');
clim([-pi,pi]); colormap(gca,wheelmap);
hex(zS,10);
view(2); axis equal;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); hold off;
colorbar;

% nexttile; hold on;
% plot(real(alphas), imag(alphas), '.');
% plot(real(alpha), imag(alpha), 'x'); axis equal;
% xlim([4 7]), ylim([-1,1]); 

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

% get f
V = sqrt(sqrt(3)/2) * const_fun(N);
V = Up*Um*V;

dzV = 1i* Dbar(N,0,f1,f2)' * V;

M=1000;
[z,v] = K2X(V,0,M,e1,e2);
[~,dzv] = K2X(dzV,0,M,e1,e2);
bracket = imag( sqrt(conj(v)) .* dzv );
f = abs(bracket) ./ abs(v);
f(1,1) = f(2,1); % lol

% %%
% figure;hold on;
% surf(real(z),imag(z),abs(bracket) ./ 1500, 'EdgeColor', 'none')
% view(2); axis equal; colorbar;
% 

figure; hold on;
contourf(real(z),imag(z),f, 18)
view(2); axis equal; colorbar;
% 

% interpolate(e1,e2,f,M,0.2165,-0.8669)
% % interpolate(e1,e2,f,M,0,-0.8)
% [z,~]=K2X3(V,0,200,e1,e2);
% fun = @(z) interpolate(e1,e2,f,M,real(z),imag(z));
% figure; hold on;
% surf(real(z),imag(z), arrayfun(fun,z), 'EdgeColor', 'none');
% view(2); axis equal; colorbar;


ffun = @(z,etc) interpolate(e1,e2,f,M,z.x,z.y);
specifyCoefficients(model, "m",0,"d",0,"c",1,"a",0,"f",ffun);

generateMesh(model,"Hmax",0.01);
result=solvepde(model);

figure; hold on;
pdeplot(result.Mesh, XYData=result.NodalSolution, Contour="on", ColorMap="default", Levels=18); axis equal;
title("$\Phi$", 'Interpreter', 'latex')
% saveas(gcf,'./results/Phi_solution.png')


%% First 10 eigenfunctions

num=10;
[~,s,V] = svds(D + 17*U + K*speye(size(D)), num, 'smallest');
disp([(1:num)', diag(s)])

fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
til=tiledlayout(2,5,'TileSpacing','tight');


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


nexttile; hold on;
contourf(real(z), imag(z), max(log(abs(v1)) ./ 17,minlevel), levels);
hex(zS); axis equal; colorbar;
xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]);

 

end

%% Log eigenfunction / alpha   smoothened, laplacian

alpha = 20;
[~,~,V] = svds(D + alpha*U + K*speye(size(D)), 1, 'smallest');
V = V(1:length(V)/2); V = V / norm(V);
[z,v] = K2X(V,K,300,e1,e2);


rad = 8;
chi = ones(rad,rad); chi(300,300) = 0; %chi=circshift(chi, [-rad/2, -rad/2]);
v_conv =  ifft2(fft2(abs(v)) .* fft2(chi)) ./ rad^2;
v_conv = log(abs(v_conv)) ./ alpha;
% v_conv = max(log(abs(v)) ./ alpha, -0.5);

V_conv =X2K(v_conv,z,0,N);
[~,Lap_v_conv]=K2X(Dbar(N,0,f1,f2) * Dbar(N,0,f1,f2)' * V_conv, 0, 300, e1,e2);



figure;
til=tiledlayout(1,4,'TileSpacing','compact');

minlevel=-1;
levels=linspace(minlevel,0.25,40);

nexttile;
hold on;
title('$\chi$ where $\int \chi = 1$', 'Interpreter', 'latex');
contourf(real(z), imag(z), -chi ./ 2, levels);
axis equal; clim([minlevel 0.25]); 


nexttile;
hold on;
title('$\alpha^{-1} \log |u_K|$', 'Interpreter', 'latex');
contourf(real(z), imag(z), max(log(abs(v)) ./ alpha, minlevel), levels);
clim([minlevel 0.25]); axis equal;

nexttile;
hold on;
title('$\alpha^{-1} \log (|u_K| \ast \chi)$', 'Interpreter', 'latex');
contourf(real(z), imag(z), max(v_conv, minlevel), levels);
axis equal; clim([minlevel 0.25]); colorbar;

% nexttile;
% hold on;
% title('$\Delta(\alpha^{-1} \log (|u_K| \ast \chi))$', 'Interpreter', 'latex');
% surf(real(z), imag(z), real(Lap_v_conv),'EdgeColor', 'None');
% axis equal; 
% view(2);
% % hex(zS); axis equal; colorbar;
% ylim([-1.721,0]), xlim([-0.5,0.5]); 
%  colorbar;

 nexttile;
hold on;
title('max$(-3,\Delta(\alpha^{-1} \log (|u_K| \ast \chi)) )$', 'Interpreter', 'latex');
contourf(real(z), imag(z), min(max(real(Lap_v_conv),-3), 18),15);
axis equal; 
view(2);
% hex(zS); axis equal; colorbar;
ylim([-1.721,0]), xlim([-0.5,0.5]); 
 colorbar;

% saveas(gcf,'./results/Laplacian_of_phi_1.png')

%% Conjugate by multiplication with e^(phi/h)

M = 300;

% choice of phi

% {q, bar q}



% eigenfunction


%% Solve laplacian phi = {q, bar q}/|xi|^2
V = Up * Um * const_fun(N);
dzV = 1i* Dbar(N,0,f1,f2)' * V;
M=600;
[Z,v1] = K2X(V,0,M,e1,e2);
[~,dzv] = K2X(dzV,0,M,e1,e2);
bracket = imag( sqrt(conj(v1)) .* dzv );
bracket = abs(bracket) ./ abs(v1);

% figure;
% contourf(real(z), imag(z), bracket, 20); axis equal;
% %%

Bracket = X2K(bracket, z,0,N);
% bracket = bracket - max(bracket, [], 'all');

Lap = Dbar(N,0,f1,f2) * Dbar(N,0,f1,f2)'; Lap(2*N^2 + 2*N+1,2*N^2 + 2*N+1) = 1;

[z,v] = K2X(Inv(Lap) * Bracket,0, M,e1,e2);

figure;
contourf(real(z), imag(z), abs(v), 20);
axis equal;

%% Animation overlaying log |u|/alpha and log |2hDbar u|/alpha   as alpha grows



vid=VideoWriter('.\results\log_norm_scalar_overlay_1.mp4','MPEG-4'); open(vid);
fig=figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

alphas = 5:0.1:20;

for id=1:length(alphas)
    clf(fig);
    tl=tiledlayout(2,2,'TileSpacing','compact');
    title(tl, ['$\alpha = ' sprintf('%.1f', alphas(id)) '$'], 'Interpreter', 'latex');

    [~,~,V] = svds((4*Dbar(N,0,f1,f2)^2 -alphas(id)^2 * Up*Um), 1, 'smallest');
    hdbarV = 2*Dbar(N,0,f1,f2) * V ./ alphas(id);
    [z,v] = K2X3(V,0,300,e1,e2);
    [~,dv] = K2X3(hdbarV,0,300,e1,e2);

    nexttile(1);
    hold on;
    title('$h\log |u_{0}|$, where $P(\alpha)u_0 = 0$', 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.5,40);
    contourf(real(z), imag(z), max(log(abs(v)) ./ alphas(id), minlevel), levels);
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]); 

    nexttile(3);
    hold on;
    title('$h\log |2hD_{\bar z} u_{0}|$', 'Interpreter', 'latex');
    minlevel=-1;
    levels=linspace(minlevel,0.5,40);
    contourf(real(z), imag(z), max(log(abs(dv)) ./ alphas(id), minlevel), levels);
    hex(zS); axis equal; colorbar;
    xlim([-0.63,0.63]), ylim([-0.63,0.63]); clim([minlevel 0.25]); 


    M = 600;
    [~,v] = K2X(V,K,M,e1,e2);
    [~,dv] = K2X(hdbarV,0,M,e1,e2);
    
    nexttile
    hold on
    title('restricted to line from $0$ to $-\sqrt 3 i$', 'Interpreter', 'latex');
    plot((0:M-1)'./M, log(abs(diag(v))) ./ alphas(id), 'Color', 'blue');
    plot((0:M-1)'./M, log(abs(diag(dv))) ./ alphas(id), 'Color', 'red');
    xlim([0 1])
    ylim([-1 0.2])
    xline(1/3); xline(2/3); yline(0, 'Color', 'black')
    legend('u', '2hDbar u')
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


%% Misc functions

function fQ=interpolate(e1,e2,f,M,xQ,yQ)
    T = [real(e1), real(e2); imag(e1), imag(e2)];
    v = T^(-1) * [xQ; yQ];
    v = v - floor(v); v = ceil(v*M); v(v==0) = 1;
    % disp(v)
    ind = v(1,:) + M*(v(2,:)-1);
    fQ = f(ind);
end

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

% (2K+) Potential with TBG rotation symmetry containing K + n1 f1 + n2 f2

function U=sym_potential(N,n1,n2)

    w = exp(2i*pi/3);

    U = (-4i*pi/3)* (fourier_shift(N,n1-1,n2+1) + w*fourier_shift(N,-n2-1,n1-n2) + w^2*fourier_shift(N,n2-n1, -n1+1));
end

function V=const_fun(N)
    V = kron((-N:N).'==0, (-N:N).'==0);
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
    v = v ./ sqrt(sqrt(3)/2); % normalize for fundamental domain size

end

% Inverse

function V=X2K(v,z,k,N)

    v = v ./ exp(0.5i*(k'*z+k*conj(z)));
    v = v .* sqrt(sqrt(3)/2);
    
    V = ifft2(v);
    
    V = circshift(V, [N,N]);
    V = V(1:2*N+1, 1:2*N+1);
    V = V.'; V = V(:);

end


% Variations

function [z,v1,v2] = K2X2(V,k,M,e1,e2)

    K = 4*pi/3;

    [z,v1] = K2X(V(1:length(V)/2),k-K,M,e1,e2);
    [~,v2] = K2X(V(1+length(V)/2:end),k+K,M,e1,e2);

    phase=v1(1) / abs(v1(1));
    v1=v1 ./ phase; v2=v2./phase;

end

% repeat over some copies of fundamental domain

function [z,v] = K2X3(V,k,M,e1,e2)

    [z0,v0]=K2X(V,k,M,e1,e2);

    z = [z0-e1-e2, z0-e2; z0-e1, z0];
    
    p2 = exp(-0.5i*(k'*e2+k*e2')); p1 = exp(-0.5i*(k'*e1+k*e1'));
    v = [v0 *p1 *p2, v0 *p2; v0 *p1, v0];
    v = v ./ (v0(1) / abs(v0(1)));

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

function [V0,V1,V2] = PROJC(N,V,s)

    w = exp(2i*pi/3); R = ROT2(N,s);

    V0= (R^0 + w'^0 .* R + w^0 .* R^2)/3;
    V1= (R^0 + w'^1 .* R + w^1 .* R^2)/3;
    V2= (R^0 + w'^2 .* R + w^2 .* R^2)/3;
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



