function [boundary,unsure,areas,counts,z,bra]=bracket_area(pot,M,M2,num_iter)

[F,DF,F0,DF0,x,y]=BRA2(pot,M);


w = exp(2i*pi/3);
z = (w^2*x - w*y) * 3/(2*pi);
bra = imag(sqrt(F0));



x = x(:); y=y(:);
Epsilon = y(2) - y(1);

M = 2*M+1;
indx = reshape(1:M^2, M,M);
up = circshift(indx, [1 0]);
left = circshift(indx, [0 1]);
down = circshift(indx, [-1 0]);
right = circshift(indx, [0 -1]);
Neigh = [up(:), left(:), down(:), right(:)];

v = F0(:); dv = DF0(:);
v1 = v(up(:));
v2 = v(left(:));
v3 = v(up(left(:)));

N = 0; areas = [];
boundary = [];

for iter=1:num_iter

    disp('iteration'); disp(iter);

% A - blocks not intersecting zero set of bracket
% B - blocks intersecting zero set
% C - blocks for which not sure

% At each iteration put some blocks from C in A or B, then subdivide remaining blocks of C


% In the end we return

% boundary - all blocks in B
% unsure - blocks in C at the end
% areas, counts - areas of connected components of union of all blocks in A

% Could make the boundary less jagged looking by not putting blocks in B but
% keeping them in C for several iterations.

% Can be made more accurate by being more selective about
% what is put in B.
% However, this has not been an issue with the zero sets so far.

indx = (1:length(x))';

dv = max([dv, dv(Neigh(:,1)), dv(Neigh(:,2)), dv(Neigh(:,3)), dv(Neigh(:,4))], [],  2);
C = indx > N  & ...
        abs(v) <= dv * Epsilon * sqrt(2);     % a little sus %%%%%%%%%%%%%%%%%%%

% Use Matlab graph library to get connected components

E1 = indx > N & ~EDGEP(v(indx), v1(indx));
E2 = indx > N & ~EDGEP(v(indx), v2(indx));
E3 = indx > N & ~EDGEP(v1(indx), v3(indx));
E4 = indx > N & ~EDGEP(v2(indx), v3(indx));

E = E1 | E2 | E3 | E4;

A = ~C & ~E;
B = ~C & E;

boundary = [boundary; (w^2*x(B) - w*y(B)) * 3/(2*pi)];


A1 = A & A(Neigh(indx,1)); 
A2 = A & A(Neigh(indx,2));
A3 = A & A(Neigh(indx,3)); 
A4 = A & A(Neigh(indx,4));

G = graph([find(C); find(A1);    find(A2);     find(A3);      find(A4)]', ...
          [find(C); Neigh(A1,1); Neigh(A2,2);  Neigh(A3,3);   Neigh(A4,4)]' );

[bins,binsizes] = conncomp(G); 
bins = bins';

for k=1:N
    binsizes(bins(k)) =  binsizes(bins(k)) + areas(k) - 1;
end


%
order = [find(A); find(C)];
order_invert(order) = 1:length(order);

[c,~,bins] = unique(bins(order), 'stable');
binsizes = binsizes(c);

N = length(binsizes) - length(find(C));
n = length(find(C));

areas =  binsizes(1:N);



Neigh2=zeros(N+n, 4);

for j=1:4
    Neigh2(1:N,j) = 1:N;
    next = Neigh(C,j);
    nextB = B(next);

    Neigh2(N+find(nextB), j) = N+find(nextB);
    Neigh2(N+find(~nextB), j) = bins(order_invert(next(~nextB)));
end

x = x(C); y=y(C);


% Subdivide blocks in C
if iter ~= num_iter
    Epsilon=Epsilon/M2;
    
    xx = zeros(n*M2^2,1);
    yy = zeros(n*M2^2,1);
    Neigh = zeros(N+n*M2^2,4);
    
    for j=1:4
    Neigh(1:N,j)=1:N;
    end
    
    for k1=0:M2-1
        for k2=0:M2-1
            k = k1*M2+k2;
            xx(k*n+1:(k+1)*n) = x-Epsilon*k1;
            yy(k*n+1:(k+1)*n) = y-Epsilon*k2;
    
            for j=1:4
                kk1 = k1; kk2 = k2;
                if j==1
                    kk2 = kk2 + 1;
                elseif j==2
                    kk1 = kk1 + 1;
                elseif j==3
                    kk2 = kk2 - 1;
                elseif j==4
                    kk1 = kk1 - 1;
                end
                kk1 = mod(kk1,M2); kk2 = mod(kk2,M2);
                kk = kk1*M2+kk2;
                
                if j==1 && k2 == M2-1 || j==2 && k1 == M2-1 || j==3 && k2 == 0 || j==4 && k1 == 0
                    next=true;
                else
                    next=false;
                end
    
                if next
                    X= Neigh2(N+1:end,j) <= N;
                    Y= Neigh2(N+1:end,j) == (N+1:N+n)';

                    Neigh(N+k*n+find(X),j) = Neigh2(N+find(X),j);
                    Neigh(N+k*n+find(Y),j) = N+k*n+find(Y);
                    
                    Neigh(N+k*n+find(~X & ~Y),j) = kk*n+Neigh2(N+find(~X & ~Y),j);
                else
                    Neigh(N+k*n+1:N+(k+1)*n ,j) = N+kk*n+1: N+(kk+1)*n;
                end
            end
        end
    end
    x = [zeros(N, 1); xx]; y = [zeros(N, 1); yy];
    disp('calculating bracket for num values:')
    disp(length(x));

    v = F(x(:), y(:));
    v1 = F(x(:), y(:)-Epsilon);
    v2 = F(x(:)-Epsilon, y(:));
    v3 = F(x(:)-Epsilon, y(:)-Epsilon);
    dv = DF(x(:), y(:));

    areas = areas * M2^2;
end

end

unsure = (w^2*x - w*y) *3/(2*pi);


% process areas
areas = areas / sum(areas);
areas = areas(areas > 0.00001);

disp(areas);

% identify "same" areas

[areas,~,ic] = uniquetol(sqrt(areas),15/M, 'DataScale', 1); %%%%%%%%%%%%%%%%%%%%%%%%
areas = areas' .^ 2;
counts=accumarray(ic,1);

areas = areas(length(areas):-1:1); counts = counts(length(areas):-1:1);

end



function dbound=dbound   %%%%%%%%%%%%%%%%%%%%%%%%%


end

function edge=EDGEP(z1,z2)        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = real(z1); b = imag(z1); c = real(z2); d = imag(z2);

    edge = b.*d > 0 | abs(d) .* a + abs(b) .* c <= 0;

    % edge = b .* d > 0 |  (b .* d == 0 & abs(d) .* a + abs(b) .* c < 0) | (b .* d < 0  &  abs(d) .* a + abs(b) .* c <= 0);

end


function [F,DF,F0,DF0,x,y]=BRA2(pot,M)
    
    U = pot2U(pot,12); % for more general potentials make sure N big enough
    
    flipU=@(U)flip(flip(U,1),2);
    conjU=@(U)conj(flipU(U));
    
    V=conv2(U,flipU(U));
    zV=DZ(V);
    xV=DX(V);
    yV=DY(V);
    xzV=DX(zV);
    yzV=DY(zV);
    
    F = @(x,y) EVF(conjU(V),x,y) .* EVF(zV,x,y).^2;
    
    DF= @(x,y) sqrt(  abs( EVF(conjU(xV),x,y) .* EVF(zV,x,y).^2 ...
        + 2*EVF(conjU(V),x,y).*EVF(xzV,x,y).*EVF(zV,x,y) ).^2 ...
        +  abs( EVF(conjU(yV),x,y).*EVF(zV,x,y).^2 ...
        + 2*EVF(conjU(V),x,y).*EVF(yzV,x,y).*EVF(zV,x,y) ).^2   );
    
    
    B = conv2(conjU(V), conv2(zV,zV));
    [F0,x,y] = F2P(B,M); 
    
    xB = conv2(conjU(xV), conv2(zV,zV)) + 2*conv2(conjU(V), conv2(xzV,zV));
    yB = conv2(conjU(yV), conv2(zV,zV)) + 2*conv2(conjU(V), conv2(yzV,zV));
    F1 = F2P(xB,M); F2 = F2P(yB,M);
    
    DF0 = sqrt( abs(F1).^2 + abs(F2).^2  );

end


function U=pot2U(pot,N)
    U = U_mat(0,0,N) * 0;

    for id=1:3:length(pot)
        U = U + U_mat(pot(id), pot(id+1), N) * pot(id+2);
    end
    U = 1i*U;
end

% m1 e1 + m2 e2 + K = e1/3 (3m1 - 1) + e2/3 (3m2 + 1)

function U=U_mat(m1,m2,N)
    w = exp(2i*pi/3);
    U = zeros(2*N+1);
    a1 = 3*m1-1; a2 = 3*m2+1;
    U(N+1+a1,N+1+a2)=1; 
    U(N+1-a2,N+1+a1-a2)=w; 
    U(N+1+a2-a1,N+1-a1)=w^2;
end


function DX=DX(V)
    N=round((size(V,1)-1)/2);
    DX = V .* kron((-N:N),ones(2*N+1,1));
end

function DY=DY(V)
    N=round((size(V,1)-1)/2);
    DY = V .* kron((-N:N)',ones(2*N+1,1)');
end

function DZ=DZ(V)
    w=exp(2i*pi/3);
    K=4*pi/3;
    DZ = (w*K*DX(V) + w'*K*DY(V)) / 2 ;
end


function v=EVF(V,x,y)
    N = round((size(V,1)-1)/2);
    v = 0*x;
    [j1,j2]=find(V);
    for k=1:length(j1)
        a1=j1(k)-N-1; a2=j2(k)-N-1;
        v = v + V(j1(k),j2(k)) .* exp(1i*(a1*x+a2*y));
    end
end


function [v,x,y] = F2P(V0,M)

    N = (size(V0,1)-1)/2;

    V = zeros(2*M+1);
    V(1:2*N+1, 1:2*N+1)=V0;
    V = circshift(V, [-N,-N]);

    v = fft2(V);

    v = circshift(v,[M,M]);
    [x,y]=meshgrid(-M:M, -M:M);
    x = 2*pi*x/(2*M+1);  y = 2*pi*y/(2*M+1);
end