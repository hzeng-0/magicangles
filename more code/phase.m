function phase = phase(x,k,F,al)


% k is a complex number, F is a vector of length 5, N is the size of the grid


  
% computes the phause as a function of x:
% xk + \alpha \int_0^x sqrt V (y)dy 
% where V ( x) = (F0 + F1 cos x + F2 sin x + F3 cos 2x + F4 sin 2x)^2; 
  
  if (nargin < 3)
    F = [2,1i,0,0,0];
end 


  
  if (nargin < 3)
    al = 1;
    end


phase=F(1)*x+F(2)*sin(x)-F(3)*(cos(x)-1)+0.5*F(4)*(cos(2*x)-1)-0.5*F(5)*sin(2*x);


  phase = k*x + al*phase;


	 
end
