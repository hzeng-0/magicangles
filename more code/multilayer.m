
om = exp(2i*pi/3);

e1 = 1i*om'*4*pi/3;     % basis of Lambda
e2 = 1i*(-om)*4*pi/3;

f1 = sqrt(3)*(-om);     % dual basis of Lambda^* (<e_i, f_j> = 2pi 1_{i==j})
f2 = sqrt(3)*(-om');


N = 32;        % determines size of "truncated fourier space" we use


Up = fourier_shift(N,0,0) + om*fourier_shift(N,0,-1) + om^2 * fourier_shift(N,1,0);
Um = Up.';


% L^2_{0,0}(C;C^4)
% = L^2_{i,0} + L^2_{0,0} + L^2_{i,2} + L^2_{0,1}

A = ROT_INC(N,1,f1,f2);
B = ROT_INC(N,0,f1,f2);

R0 = ROT(N,0);
R1 = ROT(N,1);
Extend0 = @(p) (R0^0 + om'^p*R0^1 + om^p * R0^2);
Extend1 = @(p) (R1^0 + om'^p*R1^1 + om^p * R1^2);

D0 = 2*Dbar(N,0,f1,f2);
D1 = 2*Dbar(N,1i,f1,f2);


% [C,C0] = ROT_INC2(N,0,f1,f2);
% m1 = size(A, 2);
% m2 = size(B, 2);
% Pr = spdiags(ones(m2,1), 1, m2, m2+1);
% 
% V = [sparse(m1,m1),    A' * Up * (3*Extend0(0) * C + C0);
%     B' * Um * (3*Extend1(0) * A),   sparse(m2, m2+1)];
% V(2*m1+2*m2+1,2*m1+2*m2+1) = 0;
% D_inv = [ Inv(A'*D1*A),      sparse(m1,m2),    -t1*Inv(A'*D1^2 * A),   sparse(m1, m2+1);
%           sparse(m2+1, m1),  Pr'*Inv(B'*D0*B), sparse(m2+1, m1),  (1/t1)*(C0' * C0);
%           sparse(m1,m1),    sparse(m1,m2),       Inv(A'*D1*A),   sparse(m1,m2+1);
%           sparse(m2,m1),   -t1*Inv(B'*D0^2*B),  sparse(m2, m1),   Inv(B'*D0*B)*Pr];
% alphas = 1./eigs(D_inv * V, 200);


T =  B'*Um*(Extend1(0)*A)  *  Inv(A'*D1*A)   *  A'*Up*Extend0(0)*B   *   Inv(B'*D0*B);

alphas = 1./sqrt(eigs(T, 200));    % alphas in (2.27)

magicalphas = 1./sqrt(eigs( Um*Inv(2*Dbar(N,2i,f1,f2))*Up*Inv(2*Dbar(N,1i,f1,f2)) , 400));   % magic alphas for tbg

figure
hold on
title('crosses are $\alpha$ with $-1/\alpha \in Spec_{L^2_{0,0}} D_2(0)^{-1} V$, dots are magic alpha for TBG', 'Interpreter', 'Latex')
plot(real(alphas), imag(alphas), 'x')
scattermult([real(magicalphas), imag(magicalphas)], 10)



%%


D = @(alpha, k) [2*Dbar(N,k+1i,f1,f2), alpha*Up;
                 alpha*Um,  2*Dbar(N,k,f1,f2)];
ss = size(Up);
Z = sparse(ss(1), ss(2));
Tp = [speye(ss), Z;
      Z, Z];
Tm = [Z, Z;
      Z, speye(ss)];


% alpha = alphas(abs(alphas - 4.60233 - 3.399i) < 0.01);
alphas = 0.3:0.01:1.7;


k_range = (-0.5:0.05:2.5)*1i;
k_range2 = 0 + (-0.1:0.01:0.1)*1i;
numbands=2;
numbands2=1;
bands = zeros(numbands, length(k_range));
bands2 = zeros(numbands, length(k_range2));

vid=VideoWriter('.\bands_multilayer_6.mp4','MPEG-4'); open(vid);
fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);


for alpha=alphas

    for id=1:length(k_range)
        D2 = [D(alpha, k_range(id)), Tp;
            Tm,  D(0, k_range(id))];
        s = svds(D2, numbands, 'smallest');
        bands(:,id) = s;

    end

    for id=1:length(k_range2)
        D2 = [D(alpha, k_range2(id)), Tp;
            Tm,  D(0, k_range2(id))];
        s2 = svds(D2, numbands2, 'smallest');
        bands2(:,id) = s2;

    end

    clf(fig);
    hold on
    tl=tiledlayout(1,2,'TileSpacing','compact');
    title(tl,['$\alpha = $' sprintf('%.2f', alpha)], 'Interpreter', 'latex')

    nexttile
    hold on
    plot(imag(k_range), bands);
    plot(imag(k_range), -bands);
    xline(0)
    xline(1)
    xline(2)
    xticklabels(arrayfun(@(k) sprintf('%gi', k) ,  xticks   ,  'UniformOutput' ,  false))
    xlim([min(imag(k_range)), max(imag(k_range))])
    ym = max(bands(numbands,:));
    ylim([-ym*1.1,ym*1.1]);

    nexttile
    hold on
    plot(imag(k_range2), bands2);
    plot(imag(k_range2), -bands2);
    xline(0)
    xticks(imag(k_range2(1:5:end)));
    xticklabels(arrayfun(@(k) sprintf('%gi', k) ,  xticks   ,  'UniformOutput' ,  false))
    xlim([min(imag(k_range2)), max(imag(k_range2))])
    ym = max(bands2(numbands2,:));
    ylim([-ym*1.1,ym*1.1]);


    drawnow
    frame=getframe(fig);
    writeVideo(vid,frame)
end

close(vid)







% FUNCTIONS


% Derivatives --------------

function Dbar = Dbar(N,k,f1,f2)
    D0 = spdiags((-N:1:N)', 0, 2*N+1, 2*N+1); E = speye(2*N+1, 2*N+1);
    Dbar = 0.5 * (kron(E, f1 * D0) + kron(f2 * D0, E) + k*kron(E,E));
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




% Symmetries ----------

% Inclusion of V, L^2_{s*i}(C;C^1) = V + RV + R^2V (+C  if s = 0).
function Inc = ROT_INC(N,s,f1,f2)
    [yy,xx] = meshgrid(-N:N, -N:N);
    k = xx * f1 + yy * f2 + s*1i;
    k = k(:);
    vv = find( (abs(k) > 0.0001) & (angle(k) >= -pi/3-0.00001) & (angle(k) < pi/3 - 0.000001) );
    m = length(vv);
    Inc = sparse(vv, 1:m, ones(m,1), (2*N+1)^2, m);
end

function [Inc, Inc0] = ROT_INC2(N,s,f1,f2)
    [yy,xx] = meshgrid(-N:N, -N:N);
    k = xx * f1 + yy * f2 + s*1i;
    k = k(:);
    vv = find( (abs(k) > 0.0001) & (angle(k) >= -pi/3-0.00001) & (angle(k) < pi/3 - 0.000001) );
    m = length(vv);
    Inc = sparse(vv, 2:m+1, ones(m,1), (2*N+1)^2, m+1);
    Inc0 = sparse(find(abs(k) < 0.0001), 1, 1, (2*N+1)^2, m+1);
end


% For k in {0,i,-i}, R[u](x) = u(wx) acts on L^2_k(C;C^1)
function R = ROT(N,s)         % k = si
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
