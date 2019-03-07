function D2=RBFkernel(X,Y,gamma)

D2=exp(-gamma*(X*X'+Y*Y'-2*X*Y'));


% function kxy=RBFkernel(x,y,sigma)
% 
% % xt=transpose of x vector, yt= transpose of y vector
% 
% xsq = sum(x.^2);
% ysq = sum(y.^2);
%    
% Dsq = xsq + ysq -2*x'*y;
% kxy=exp(-sigma*Dsq);
