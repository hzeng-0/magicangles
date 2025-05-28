function [A00,A01,A10,A11] = monodromy_mat(z,alpha,gamma)
om = exp(2i*pi/3);
K = 4*pi/3;

f1 = (4i*pi/sqrt(3))*om;
f2 = (4i*pi/sqrt(3))*om^2;


f = @(z,w,a) cos(0.5*(z * a' + w *a));
% V = @(z,w) (2i*K^2*(f(z,w,f1+f2) + om*f(z,w,f1) + om'*f(z,w,f2))).^2;

% for BMH potential:
V = @(z,w) -2*K^2 * (f(z,w,f1+f2) + om'*f(z,w,f1) + om*f(z,w,f2));

sha = size(z);
z = z(:);

w0 = @(z) 0*z;
sol = zeros(length(z),2,2);
sol(:,1,1) = 1; sol(:,2,2) = 1;

ep = 1/1000;
for t=0:ep:0.999999
    sol(:,1,:) = sol(:,1,:) + ep*0.5i*alpha*gamma'*sol(:,2,:);
    VV = V(z+gamma,w0(z)+t*gamma'); W(:,1,1) = VV; W(:,1,2) = VV;
    sol(:,2,:) = sol(:,2,:) + ep*0.5i*alpha*gamma'* W .* sol(:,1,:);
end

A00 = reshape(sol(:,1,1), sha);
A01 = reshape(sol(:,2,1), sha);
A10 = reshape(sol(:,1,2), sha);
A11 = reshape(sol(:,2,2), sha);

end