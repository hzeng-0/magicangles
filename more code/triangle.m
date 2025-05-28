function [z,tr,logu]=triangle(alpha,ep)

% Constants

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




M = 600;
MM = 2*M;

[~,~,U1] = svds(4*Dbar(N,0,f1,f2)^2 - alpha^2*(Up*Um), 1, 'smallest');
[z,u1] = F2X3(U1,0,M,e1,e2);

[~,u2] = F2X3(Dx(N,0,f1,f2)*U1,0,M,e1,e2);
u1 = u1 + u2;

u1 = u1 ./ exp(1i*angle(u1(1,1)));



% Triangle and phase function for triangle

% ep = 0.03;
tr = real(z) < 0.5-ep & real(om*z) < 0-ep & real(om' * z) < 0-ep;
% [~,v] = F2X3(Up*Um*const_fun(N), 0,M, e1, e2);

logu = log_continue(tr, log(u1));



end

function v_new=log_continue(triangle, v)
z_ind = find(triangle);
MM = sqrt(length(triangle(:)));
v_new = v;

import java.util.LinkedList
q = LinkedList();
add(q, z_ind(1));
add(q, imag(v(z_ind(1)))); % this thing doesnt handle complex numbers

visited = false(size(triangle));
visited(z_ind(1)) = true;

while ~isEmpty(q)

    ind = remove(q);
    ref = remove(q);

    diff = ref - imag(v_new(ind));
    v_new(ind) = v_new(ind) + 2i*pi*round(diff/(2*pi));


    ind2 = floor((ind-1) / MM)+1; 
    ind1 = ind - MM*(ind2-1);
    for ii=-1:1
        for jj=-1:1
            if ~visited(ind1+ii, ind2+jj) && triangle(ind1+ii, ind2+jj)
                add(q, (ind2+jj-1)*MM+ind1+ii);
                add(q, imag(v_new(ind1, ind2)));
                visited(ind1+ii, ind2+jj) = true;
            end
        end
    end
end
end



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