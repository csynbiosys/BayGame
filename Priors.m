figure; 
% plot a uniform distribution in [0,1]
 A = 0.1*4e-2;
 B = 10*4e-2;
 meas = B-A;
 X = A:0.01:B;
%  y = 1/meas*ones(size(X));
%  plot(X,y); hold on
 xlim([A-0.5,B+0.5])
 
 % trying to plot the corresponding normal distribution
m = mean([X(1),X(end)])
s = meas/4
norm = normpdf(X,m,s);
plot(X,norm); hold on; 


