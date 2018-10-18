figure; 
% plot a uniform distribution in [0,1]
 A = 1;
 B = 5;
 meas = B-A;
 X = A:0.01:B;
 y = 1/meas*ones(size(X));
 plot(X,y); hold on
 xlim([A-0.5,B+0.5])
 
 % trying to plot the corresponding normal distribution
m = mean([X(1),X(end)])
s = meas/6;
norm = normpdf(X,m,s);
plot(X,norm); hold on; 

norm2 = normpdf(X,m,s+1);
plot(X,norm2); hold on; 