function [w1, w2, width, rule, e, RMSE] = DFNN(p, t, parameters)
% This is N-FNN training program.
% Input:
%       p is the input data, which is r by q matrix. r is the No. of input
%       t is the output data,which is s2 by q matrix. q is the No. of sample data.
%       parameters is a vector which defines the predefined parameters 
%       parameters(1)= kdmax     (max of accommodation criterion)
%       parameters(2)= kdmin     (min of accommodation criterion)
%       parameters(3)= gama      (decay constant)
%       parameters(4)= emax      (max of output error)
%       parameters(5)= emin      (min of output error)
%       parameters(6)= beta      (convergence constant)
%       parameters(7)= width0    (the width of the first rule)
%       parameters(8)= k         (overlap factor of RBF units)
%       parameters(9)= kw        (width updating factor)
%       parameters(10)= kerr     (significance of a rule)
% Output:
%       w1 is the centers of the RBF units, which is a u by r matrix
%       w2 is the weights, which is s2 by u(1+r) matrix
%       width is the widths of RBF units, which is a 1 by u matrix
% Revised 11-3-2006
% Copyright Wu Shiqian.
if nargin<2
    error('Not enough input arguments')
end
if size(p,2)~=size(t,2)
    error('The input data are not correct')
end
[r,q]=size(p);
[s2,q]=size(t);
% Setting predefined parameters
kdmax=parameters(1);
kdmin=parameters(2);
gama=parameters(3);
emax=parameters(4);
emin=parameters(5);
beta=parameters(6);
width0=parameters(7);
k=parameters(8);
kw=parameters(9);
kerr=parameters(10);
ALLIN=[];
ALLOUT=[];
CRBF=[];
%When first sample data coming
ALLIN=p(:,1);
ALLOUT=t(:,1);
% Seeting up the initial DFNN
CRBF=p(:,1);
w1=CRBF';
width(1)=width0;
rule(1)=1;
% Caculating the first out error
a0=RBF(dist(w1,ALLIN),1./width');
a01=[a0 p(:,1)'];
w2=ALLOUT/a01';
a02=w2*a01';
RMSE(1)=sqrt(sumsqr(ALLOUT-a02)/s2);
% When other sample data coming
for i=2:q
   IN=p(:,i);
   OUT=t(:,i);
   ALLIN=[ALLIN IN];
   ALLOUT=[ALLOUT OUT];
   [r,N]=size(ALLIN);
   [s r]=size(w1);
   dd=dist(w1,IN);
   [d,ind]=min(dd);
   kd=max(kdmax*gama.^(i-1),kdmin);
   % Caculating the actual output of ith sample data
   ai=RBF(dist(w1,IN),1./width');
   ai=ai/sum(ai);
   ai1=transf(ai,IN);
   ai2=w2*ai1;
   errout=t(:,i)-ai2;
   errout=sum(errout.*errout)/s2;
   e(i)=sqrt(errout);
   ke=max(emax*beta.^(i-1),emin);
   if d > kd
       if e(i) > ke
           CRBF=[CRBF IN];    % Add a new rule
           wb=k*d;
           width=[width wb];
           w1=CRBF';
           [u,v]=size(w1);
           % Caculating outputs of RBF after growing for all coming data
           A=RBF(dist(w1,ALLIN),1./width');
           A0=sum(A);
           A1=A./(ones(u,1)*A0);
           A2=transf(A1,ALLIN);
           if u*(r+1)<=N         
               % caculating error reduction rate
               tT=ALLOUT';
               PAT=A2';
               [W,AW]=orthogonalize(PAT);
               SSW=sum(W.*W)';
               SStT=sum(tT.*tT)';
               err=((W'*tT)'.^2)./(SStT*SSW');  %err is s2*u(r+1), which corresponds to w2
               errT=err';
               err1=zeros(u,s2*(r+1));
               err1(:)=errT;  % err1 corresponds to ww2,which is u*s2(r+1)
               err21=err1';   
               % err21 ,s2(r+1)*u,corresponds to w21, in which every element
               % means the importance of the correspond weight coefficient in
               % w21. The bigger, the more important.
               err22=sum(err21.*err21)/(s2*(r+1));
               err23=sqrt(err22);
               No=find(err23<kerr);
               if ~isempty(No)
                   CRBF(:,No)=[];
                   w1(No,:)=[];
                   width(:,No)=[];
                   err21(:,No)=[];
                   [uu,vv]=size(w1);
                   AA=RBF(dist(w1,ALLIN),1./width');
                   AA0=sum(AA);
                   AA1=AA./(ones(uu,1)*AA0);
                   AA2=transf(AA1,ALLIN);
                   w2=ALLOUT/AA2;         
                   outAA2=w2*AA2;
                   sse0=sumsqr(ALLOUT-outAA2)/(i*s2);
                   RMSE(i)=sqrt(sse0);
                   rule(i)=uu;
                   w2T=w2';
                   ww2=zeros(uu,s2*(r+1));
                   ww2(:)=w2T;
                   w21=ww2';  
               else
                   w2=ALLOUT/A2;   % w2 is s2*u(r+1)
                   outA2=w2*A2;
                   sse0=sumsqr(ALLOUT-outA2)/(s2*i);
                   RMSE(i)=sqrt(sse0);
                   rule(i)=u;
                   w2T=w2';
                   ww2=zeros(u,s2*(r+1));
                   ww2(:)=w2T;
                   w21=ww2';    %w21 (size: s2(r+1)*u) is very important,in which every 
                   %i-th column means the all weights to all the outputs
               end
           else 
               w2=ALLOUT/A2;   
               outA2=w2*A2;
               sse0=sumsqr(ALLOUT-outA2)/(s2*i);
               RMSE(i)=sqrt(sse0);
               rule(i)=u;
               w2T=w2';
               ww2=zeros(u,s2*(r+1));
               ww2(:)=w2T;
               w21=ww2';   
           end
       else         
           a=RBF(dist(w1,ALLIN),1./width');
           a0=sum(a);
           a1=a./(ones(s,1)*a0);
           a2=transf(a1,ALLIN);
           w2=ALLOUT/a2;
           outa2=w2*a2;
           sse1=sumsqr(ALLOUT-outa2)/(s2*i);
           RMSE(i)=sqrt(sse1);
           rule(i)=s;
       end
   else
       if e(i) > ke
           width(ind)=kw*width(ind);
           aa=RBF(dist(w1,ALLIN),1./width');
           aa0=sum(aa);
           aa1=aa./(ones(s,1)*aa0);
           aa2=transf(aa1,ALLIN);
           w2=ALLOUT/aa2;
           outaa2=w2*aa2;
           sse2=sumsqr(ALLOUT-outaa2)/(i*s2);
           RMSE(i)=sqrt(sse2);
           rule(i)=s;
       else
           aa1=RBF(dist(w1,ALLIN),1./width');
           aa01=sum(aa1);
           aa11=aa1./(ones(s,1)*aa01);
           aa21=transf(aa11,ALLIN);
           w2=ALLOUT/aa21;
           outaa21=w2*aa21;
           sse3=sumsqr(ALLOUT-outaa21)/(s2*i);
           RMSE(i)=sqrt(sse3);
           rule(i)=s;
       end
   end
end