% Static function approximation by D-FNN
% The underlying function is Hermite function shown in chapter 5
% Revised 11-3-2006
% Copyright Wu Shiqian.
clc
clear
close all

rand('seed',5);
p=rand(1,200)*8-4;
t=1.1*(1-p+2*p.^2).*exp(-p.^2./2);

[r,q]=size(p);
[s2,q]=size(t);

% Setting initial values
kdmax=2; kdmin=0.2; gama=0.977; beta=0.9; width0=2;
emax=1.1; emin=0.02; k=1.1; kw=1.1; kerr=0.0015;
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
ylabel('No of rules');

figure
plot(e,'r');
title('Actual output error');
xlabel('Input sample patterns');

figure
plot(RMSE,'r');
title('Root mean squared error (RMSE)');
xlabel('Input sample patterns');

figure;   
plot(p,t,'r+',p,outTA2,'bo');
title('Desired and actual outputs');
xlabel('Input data');