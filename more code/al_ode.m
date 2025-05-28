function al_ode = al_ode(k,F,N)


% k is a complex number, F is a vector of length 5, N is the number of Fourier
% modes
  
% computes the eigenvalue of ( D_x + k)^(-2) V ( x )  on R/2pi Z 
% using 2N+1 Fourier modes; k \notin Z 
% where V ( x) = F0 + F1 cos x + F2 sin x + F3 cos 2x + F4 sin 2x 





  


  n = -N:1:N;


  N = 2*N+1;  


  A = spdiags(ones(N,1),1,N,N);


  C1 = (A + A')/2;


  C2 = 2*C1^2 - speye(N);


  S1 = (A - A')/(2*1i);


  S2 = 2*S1*C1;
  
  A = F(1)*speye(N)+ F(2)*C1 + F(3)*S1 + F(4)*C2 + F(5)*S2;


P = sparse(diag((n + k).^-2))*A;


ee = eigs(P,200);


al_ode = 1./sqrt(ee);


al_ode = [-al_ode, al_ode];


al_ode = sort(al_ode,'ComparisonMethod','real');


end
