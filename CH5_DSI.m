% Nonlinear Dynamic System Identification in chapter 5 by D-FNN
% Revised 11-3-2006
% Copyright Wu Shiqian.
clc
clear
close all
y=zeros(1,200);
y(1)=0;
y(2)=1;
a=sin(2*pi/25);
p(:,1)=[1;0;a];
t(1)=sin(2*pi*2/25);
for k=3:200
    y(k+1)=(y(k)*y(k-1)*(y(k)+2.5))./(1+y(k).^2+y(k-1).^2)+sin(2*pi*k/25);
    x1(k-1)=y(k);
    x2(k-1)=y(k-1);
    x3(k-1)=sin(2*pi*k./25);
    p=[x1;x2;x3];
    t(k-1)=y(k+1);
end

% Setting initial values
kdmax=4;
kdmin=0.2;
gama=0.97;
beta=0.9;
width0=2;
emax=1;
emin=0.02;
k=1.1;
kw=1.1;
kerr=0.002;
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
k=1:199;
plot(k,t,'r-',k,outTA2,'ro');
title('Comparison between desired and actual outputs');
xlabel('Time t');