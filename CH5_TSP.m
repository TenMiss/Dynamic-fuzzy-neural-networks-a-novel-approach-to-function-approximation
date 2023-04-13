% Mackey—Glass time-series prediction in chapter 5 by D-FNN
% Revised 11-3-2006
% Copyright Wu Shiqian.
clc
clear
close all

x=ones(1,4000);
x(1)=1.2;
for t=18:4017
   x(t+1)=0.9*x(t)+0.2*x(t-17)/(1+x(t-17).^10);
end
x1=x(136:635);
x2=x(130:629);
x3=x(124:623);
x4=x(118:617);
p=[x1;x2;x3;x4];
t=x(142:641);
[s2,q]=size(t);
% Setting initial values
kdmax=2;
kdmin=0.25;
gama=0.98;
beta=0.95;
width0=1;
emax=1.1;
emin=0.02;
k=1.2;
kw=1.1;
kerr=0.00025;

parameters(1)= kdmax;      
parameters(2)= kdmin;     
parameters(3)= gama;     
parameters(4)= emax;      
parameters(5)= emin;     
parameters(6)= beta;      
parameters(7)= width0;    
parameters(8)= k;        
parameters(9)= kw;       
parameters(10)= kerr;  
[w1,w2,width,rule,e,RMSE] = DFNN(p,t,parameters);
TA=RBF(dist(w1,p),1./width');
TA0=sum(TA);
[u,v]=size(w1);
TA1=TA./(ones(u,1)*TA0);
TA2=transf(TA1,p);
outTA2=w2*TA2;

figure
plot(rule,'r');
title('Fuzzy rule eneration');
xlabel('Input sample patterns');

figure
plot(e,'r');
title('Actual output error e(i)');
xlabel('Sample patterns');

figure
plot(RMSE,'r');
title('Root mean squared error');
xlabel('Input sample patterns');

figure; 
k=1:500;  
plot(k,t,'r-',k,outTA2,'r--');
title('Comparison between desired and actual outputs in training');
xlabel('Time t');

% Testing result
x1=x(636:1135);
x2=x(630:1129);
x3=x(624:1123);
x4=x(618:1117);
ALLIN=[x1;x2;x3;x4];
tt=x(642:1141);
A=RBF(dist(w1,ALLIN),1./width');
SA=sum(A);
[u,v]=size(w1);
A1=A./(ones(u,1)*SA);
A2=transf(A1,ALLIN);
A3=w2*A2;
sse=sumsqr(tt-A3)/(s2*500);
rmse=sqrt(sse);
k=501:1000;
figure;   
plot(k,tt,'r-',k,A3,'r--');
title('Comparison between desired and predicted outputs');
xlabel('Time t');

figure;   
plot(k,tt-A3,'r');
title('Prediction error');
xlabel('Time t');
