function comp = comp(k,F,N)


% k is a complex number, F is a vector of length 5, N is the number of
% n


  
% computes the eigenvalue of ( D_x + k)^(-2)( 2 + i V ( x ) ) on R/2pi Z 
% using 2N+1 Fourier modes; k \notin Z 
% where V ( x) = F0 + F1 cos x + F2 sin x + F3 cos 2x + F4 sin 2x 
  
  if (nargin < 2)
    F = [2,1i,0,0,0];
end 


  
  if (nargin < 3)
    N = 20;
end
  
n = -N:1:N;


y = linspace(0,2*pi,10^5);


Fy = F(1) + F(2)*cos(y) + F(3)*sin(y) + F(4)*cos(2*y) + F(5)*sin(2*y);


W0 = 10^-5*sum(sqrt(Fy));


comp = W0^-1*(n+k)+10^-12*1i; 


comp = [ comp, W0^-1*(n-k)+10^-12*1i];


comp = sort(comp,'ComparisonMethod','real');


	 
end
