%--------calculation of 




empirical prob------------
% Generating sample from Trinomial(n,p1,p2)
%-------------------------------------------------
clc; clear all;close all;
n=10;       m=500;
p1=1/3;     q1=1-p1;
p2=1/4;     q2=1-p2;
p3=1-p1-p2; q3=1-p3;
%-------------Simulation---------------------------
x=binornd(n,p1,m,1);
y=zeros(m,1);
for i=1:m;
    y(i)=binornd(n-x(i),p2/(1-p1),1,1);
end
s=[x y];
%-------------Estimation--------------------------
smean=mean(s);  pmean=[n*p1 n*p2];
svcm=cov(s);    pvcm=[  n*p1*q1     -n*p1*p2; 
                        -n*p1*p2    n*p2*q2];
fprintf('\n\t Sample mean       :');disp(smean);
fprintf('\n\t population mean   :');disp(pmean);
fprintf('\n\t Sample and population V-C-M:\n'); 
disp([svcm pvcm]);
%-------------Calculation of empirical probabilities--
N=zeros(n+1,n+1);
for i=1:n+1
    for j=1:n+1;
        N(i,j)=sum(x==i-1 & y==j-1);
    end
end
p=N/m;
set(gcf,'color',[1 1 1])
bar3(p)
%------Simulation and estimation for X|Y=y0---------
y0=2;
cs=binornd(n-y0,p1/(1-p2),m,1);
scm=mean(cs);    pcm=(n-y0)*p1/(1-p2);
scv=var(cs);     pcv=(n-y0)*p1*p3/(1-p2)^2;
fprintf('\n\t Sample mean of X|Y=%d     : %f',y0, scm);
fprintf('\n\t popula mean of X|Y=%d     : %f',y0, pcm);
fprintf('\n\t Sample variance of X|Y=%d : %f',y0, scv);
fprintf('\n\t popula variance of X|Y=%d : %f\n\n',y0, pcv);
%------Simulation and estimation for Y|X=x0---------
clear cs;
x0=1;
cs=binornd(n-x0,p2/(1-p1),m,1);
scm=mean(cs);    pcm=(n-x0)*p2/(1-p1);
scv=var(cs);     pcv=(n-x0)*p2*p3/(1-p1)^2;
fprintf('\n\t Sample mean of Y|X=%d     : %f',x0, scm);
fprintf('\n\t popula mean of Y|X=%d     : %f',x0, pcm);
fprintf('\n\t Sample variance of Y|X=%d : %f',x0, scv);
fprintf('\n\t popula variance of Y|X=%d : %f',x0, pcv);
%------------END PROGRAM------------------------------
%-----------------------------------------------------
%-----------------------------------------------------
%-----------------------------------------------------
%-----------------------------------------------------
%-----------------------------------------------------
%-------------------SAMPLING FROM BPOISSON----------------
clc;clear all;close all;
l1=2;l2=1;l12=1;
n=1000;
u1=poissrnd(l1,n,1);    u2=poissrnd(l2,n,1);
u12=poissrnd(l12,n,1);
x1=u1+u12;  x2=u2+u12;
s=[x1 x2];
for i=0:max(x1);
    for j=0:max(x2);
        F(i+1,j+1)=sum(x1==i&x2==j);
    end
end
fprintf('\n\t empirical counts are :\n');
disp(F);
ep=F/n;
bar3(ep);
sm=mean(s); pm=[l1+l12 l2+l12];
svcm=cov(s);  pvcm=[l1+l12 l12;l12 l2+l12];
fprintf('\n\t sam and pop means are : \n');
disp([sm pm]);
fprintf('\n\t sam and poop vcm are : \n');
disp([svcm pvcm]);
%-----conditional sampling x1|x2----------------
gx2=2;
U1=poissrnd(l1,n,1);
V=binornd(gx2,l12/(l2+l12),n,1);
sc=U1+V;
for i=0:max(sc);
    cF(i+1)=sum(sc==i);
end
bar(cF/n);
smean=mean(sc);    pmean=gx2*l12/(l2+l12)+l1;
svar=var(sc);   pvar=gx2*l12/(l2+l12)^2+l1;
fprintf('\n\t smean=%d :%f\n',gx2,smean);
fprintf('\n\t pmean=%d :%f\n',gx2,pmean);
fprintf('\n\t svar=%d :%f\n',gx2,svar);
fprintf('\n\t pvar=%d :%f\n',gx2,pvar)
%----------------------------------------------------
%----------------------------------------------------
%----------------------------------------------------
%----------------------------------------------------
%-------------SAMPLING BEXP--------------------------
clc;clear all;close all;
n=1000;
l1=1;l2=2;l12=5;
x1=exprnd(1/l1,n,1);x2=exprnd(1/l2,n,1);
x12=exprnd(1/l12,n,1);
y1=min(x1,x12);    y2=min(x2,x12);
y=[y1 y2];
sm=mean(y);
pm=[1/l1+l12 1/l2+l12];
svcm=diag(cov(y))';
pvcm=[1/(l1+l12)^2 1/(l2+l12)^2];
fprintf('\n\t bivariate exponential :\n')
fprintf('\n\t sample mean \t population mean \n\n');
disp([sm' pm']);
fprintf('\n\t sample vcm \t population vcm \n\n')
disp([svcm pvcm]);
f=1-exp(-((y1/l1+y2/l2+max(y1,y2)/l12)));
scatter3(y1,y2,f)
%-----------------------------------------------------
%-----------------------------------------------------
%-----------------------------------------------------
clc;clear all;close all;
m1=0;m2=0;s1=1;s2=1;r=0;
syms x y
z1=(x-m1)/s1;   z2=(y-m2)/s2;
k=(2*pi*s1*s2*sqrt(1-r*r))^-1;
f=k*exp(-(z1^2+z2^2-2*r*z1*z2)/(2*(1-r*r)));
set(gcf,'color',[1 1 1]);
subplot(1,1,1);
ezsurf(f);
%---------------------------------------------------
%---------------------------------------------------
%---------------------------------------------------
%-----------MVND(model sampling)-----------------------
clc;clear all;close all;
n=100;   p=5;
m=[4.32 14.01 1.95 2.17 2.45]';
s=[ 4.308  1.683  1.803  2.155  -0.253;
    1.683  1.786  0.588  0.177  0.176;
    1.803   0.588  0.81  1.065   -0.158;
    2.155  0.177  1.065  1.970  -0.357;
    -0.253  0.176  -0.158  -0.357  0.504];
%-------Sample of Np(m,s) ---------------------------
c=chol(s)';         z=normrnd(0,1,p,n);
for i=1:n
    x(:,i)=m+c*z(:,i);
end%---------parameter estimation----------------------
m_hat=mean(x,2);    A=x*(eye(n,n)-ones(n,n)/n)*x';
s_ue=A/(n-1);   s_mle=A/n;
fprintf('\n\t The ue of mean and vcm are is:\n');
disp([m_hat s_ue]);
%------correlation matrix--------------------------
D_hat=diag(diag(s_mle)); D=diag(diag(s));
R=D^(-0.5)*s*D^(-0.5);
R_hat=D_hat^(-0.5)*s_mle*D_hat^(-0.5);
fprintf('\n\t The pop and sample corr matrix are :\n');
disp(R);
disp(R_hat);
%---------------------------------------------------
%---------------------------------------------------
%---------------------------------------------------
%---------------conditional sampling ---------------
clc;clear all;close all;
% sampling from X1|X2=x2
%-------------------------------------------------
p=6;    n=100;      q=3;
m=[9.47,25.56,13.25,1.44,27.29,8.80]';
s=[02.57 00.85 01.56 01.79 01.33 00.42
    00.85 37.00 03.34 13.47 07.59 00.52
    01.56 03.34 08.44 05.77 02.00 00.50 
    01.79 13.47 05.77 34.01 10.50 01.77
    01.13 07.59 02.00 10.50 23.01 03.43
    00.42 00.52 00.50 01.77 03.43 04.59];
x2=[2.15 30 9]';
%-------------manipulation--------------------
m1=m(1:q);              m2=m(q+1:p);
s11=s(1:q,1:q);         s12=s(1:q,q+1:p);
s21=s12';               s22=s(q+1:p,q+1:p);
%--------------sample from X1|X2=x2------------
cm=m1+s12*inv(s22)*(x2-m2);
cs=s11-s12*inv(s22)*s21;
x=mvnrnd(cm,cs,n)';
cm_hat=mean(x,2);
cs_hat=cov(x');
%----------o\p---------------------------------
fprintf('\n\t  population and sample mean of X1|X2=x2 :\n');
disp([cm cm_hat]);
fprintf('\n\t  population V-C-M of X1|X2=x2:\n');
disp(cs);
fprintf('\n\t  Sample V-C-M of X1|X2=x2:\n');
disp(cs_hat);
%-------------------------------------------------
% sampling from X2|X1=x1
%-------------------------------------------------
x1=[9 24 12]';
cm_1=m2+s21*inv(s11)*(x1-m1);
cs_1=s22-s21*inv(s11)*s12;
x=mvnrnd(cm_1,cs_1,n)';
%------------estimation------------------------
cm_hat_1=mean(x,2);
cs_hat_1=cov(x');
%----------o\p---------------------------------
fprintf('\n\t  population and sample mean of X2|X1=x1 :\n');
disp([cm_1 cm_hat_1]);
fprintf('\n\t  population V-C-M of X2|X1=x1:\n');
disp(cs_1);
fprintf('\n\t  Sample V-C-M of X2|X1=x1:\n');
disp(cs_hat_1);
%-------------------------------------------------
% sampling from X2 X4 x6)|(X1 X3 X5)=(x1 x3 x5)
%-------------------------------------------------
x2=[10 11 30]';
%-------------manipulation--------------------
m1=m([2 4 6]);              m2=m([1 3 5]);
s11=s([2 4 6],[2 4 6]);     s12=s([2 4 6],[1 3 5]);
s21=s12';                   s22=s([1 3 5],[1 3 5]);
%--------------sample from X1|X2=x2------------
cm_2=m1+s12*inv(s22)*(x2-m2);
cs_2=s11-s12*inv(s22)*s21;
x=mvnrnd(cm,cs,n)';
%------------estimation------------------------
cm_hat_2=mean(x,2);
cs_hat_2=cov(x');
%----------o\p---------------------------------
fprintf('\n\t  population and sample mean of (X2 X4 x6)|(X1 X3 X5)=(10 11 30) :\n');
disp([cm_2 cm_hat_2]);
fprintf('\n\t  population V-C-M of (X2 X4 x6)|(X1 X3 X5)=(10 11 30):\n');
disp(cs_2);
fprintf('\n\t  Sample V-C-M of (X2 X4 x6)|(X1 X3 X5)=(10 11 30):\n');
disp(cs_hat_2);
%----------------------------------------------------
%----------------------------------------------------
%---------------------------------------------------
%------HOTLLING T^2 STATISTICS----------------------
%-----------one sample problem----------------------
clc;clear all;close all;
p=5;n=100;
m=[4.32,14.01,1.95,2.17,2.45]';
mo=[5.52,15.47,2.51,2.96,3.30]';
s=[ 4.308  1.683  1.803  2.155  -0.253;
    1.683  1.786  0.588  0.177  0.176;
    1.803   0.588  0.81  1.065   -0.158;
    2.155  0.177  1.065  1.970  -0.357;
    -0.253  0.176  -0.158  -0.357  0.504];
%----------manupulation-----------------------------
x=mvnrnd(m,s,n)';
m_hat=mean(x,2);  s_hat=cov(x');
%--------- testing of hypothesis -----------------
t2=n*(m_hat-mo)'*inv(s_hat)*(m_hat-mo);
f_c=((n-p)/(n-1)*p)*t2;
p_value=1-fcdf(p,n-p,f_c);
%--------------- o/p-------------------------------
fprintf('\n\t sample mean :\n');
disp(m_hat);
fprintf('\n\t hotelling T^2 statistics :\n');
disp(t2);
fprintf('\n\t F calculated value :\n');
disp(f_c);
fprintf('\n\t sample v-c-m :\n');
disp(s_hat);
fprintf('\n\t p_value :\n');
disp(p_value);
if p_value<0.05
    fprintf('\n\t concusion: reject H0 :\n')
else
    fprintf('\n\t conclusion: fail to reject H0 or accept H0 :\n')
end
%---------------------------------------------------------
%--------------------------------------------------------
%-------------------------------------------------------
%-----hotelling one sample when data given---------------
clc;clear all;close all;
load data_sweat;
[p n]=size(x);
alpha=0.05; mo=[40 50 10]';
%--------manupulation---------------------------------
m_hat=mean(x,1)';    s_hat=cov(x);    % some mistake due to data store so reverse procedure than other
t2=n*(m_hat-mo)'*inv(s_hat)*(m_hat-mo);
f_c=(n-p)*t2/((n-1)*p);
p_value=1-fcdf(p,n-p,f_c);
%----------------o/p---------------------------------
fprintf('\n\t sample mean :\n');
disp(m_hat);
fprintf('\n\t hotelling T^2 statistics :\n');
disp(t2);
fprintf('\n\t F calculated value :\n');
disp(f_c);
fprintf('\n\t sample v-c-m :\n');
disp(s_hat);
fprintf('\n\t p_value :\n');
disp(p_value);
if p_value<0.05
    fprintf('\n\t concusion: reject H0 :\n')
else
    fprintf('\n\t conclusion: fail to reject H0 or accept H0 :\n')
end
%--------------------------------------------------------
%--------------------------------------------------------
%-----hotelling T^2 statistics base on two sample data---
clc;clear all;close all;
load data_iris;
[p n]=size(x);
n1=50;n2=50;n=n1+n2;
X=x(:,1:50);   Y=x(:,51:100);
xbar=mean(X,2); ybar=mean(Y,2);
s1=cov(X'); s2=cov(Y');
A=(n1-1)*s1+(n2-1)*s2;
s_pooled=A/n-2;
D2=(xbar-ybar)'*inv(s_pooled)*(xbar-ybar);
T2=n1*n2*D2/(n1+n2);    f_c=(n-p-1)*T2/((n-2)*p);
p_value=1-fcdf(p,n-p-1,f_c);
%---------------------o/p------------------------------
fprintf('\n\t Mahalnobis square distance :%f',D2);
fprintf('\n\t Hotelling T2 =%f',T2);
fprintf('\n\t calculated F=%f',f_c);
fprintf('\n\t p_value =%f',p_value);
if p_value<0.05
    fprintf('\n\t concusion: reject H0 :\n')
else
    fprintf('\n\t conclusion: fail to reject H0:\n')
end
%-------------------------------------------------------
%-------------------------------------------------------
%---------CANONICAL CORRELATION-------------------------
clc;clear all;close all;
p=4; q=2;
s=[8 2 3 1
    2 5 -1 3
    3 -1 6 -2
    1 3 -2 7];
s11=s(1:q,1:q);     s12=s(1:q,q+1:p);
s21=s12';   s22=s(q+1:p,q+1:p);
r1=s11^(-0.5)*s12*inv(s22)*s21*s11^(-0.5);
r2=s22^(-0.5)*s21*inv(s11)*s12*s22^(-0.5);
[p1 l1]=eig(r1);    [p2 l2]=eig(r2);
U=s11^(-0.5)*p1;    V=s22^(-0.5)*p2;
r0=sqrt(diag(l1));
%-----------o/p-----------------------------------------
fprintf('\n\t First cc coef=%f',r0(1));
fprintf('\n\t Coefficient of canonical variable in col :\n');
disp([U(:,1) V(:,1)]);
fprintf('\n\t Second cc coef=%f',r0(2));
fprintf('\n\t Coefficient of canonical variable in col :\n');
disp([U(:,2) V(:,2)]);
% verification
ccmat=U'*s12*V;
disp(ccmat)
%---------------------------------------------------------
%---------------------------------------------------------
%-------CANONICAL CORRELATION TO IRIS DATA----------------
clc;clear all;close all;
load data_iris;
s=cov(x');  xstd=diag(diag(s))^(-0.5);
[p n]=size(x);  q=2;
%------ manupulation -----------------------------------
s11=s(1:q,1:q);     s12=s(1:q,q+1:p);
s21=s12';   s22=s(q+1:p,q+1:p);
r1=s11^(-0.5)*s12*inv(s22)*s21*s11^(-0.5);
r2=s22^(-0.5)*s21*inv(s11)*s12*s22^(-0.5);
[p1 l1]=eig(r1);    [p2 l2]=eig(r2);
U=s11^(-0.5)*p1;    V=s22^(-0.5)*p2;
r0=sqrt(diag(l1));
%-----------o/p-----------------------------------------
fprintf('\n\t First cc coef=%f',r0(1));
fprintf('\n\t Coefficient of canonical variable in col :\n');
disp([U(:,1) V(:,1)]);
fprintf('\n\t Second cc coef=%f',r0(2));
fprintf('\n\t Coefficient of canonical variable in col :\n');
disp([U(:,2) V(:,2)]);
%--------------------------------------------------------
%-------------------------------------------------------
%------ discriminant analysis----------------------------
clc;clear all;close all;
x1=[90 115.5 94.8 91.5 117 140.1 138 112.8 99 123 81 111
    18.4 16.8 21.6 20.8 23.6 19.2 17.6 22.4 20 20.8 22 20];
x2=[105 82.8 94.8 73.2 114 79.2 89.4 96 77.4 63 81 93
    19.6 20.8 17.2 20.4 17.6 17.6 16 18.4 16.4 18.8 14 14.8];
x=[x1 x2];
x0=[85,18.8]';
n1=12;n2=12;
m1=mean(x1,2);  m2=mean(x2,2);
s1=cov(x1');    s2=cov(x2');
A=(n1-1)*s1+(n2-1)*s2;
s=A/(n1+n2-2);
%-----------------FDF--------------------------------------------------------------
d=inv(s)*(m1-m2);  %since for sample,so d instead of del.
c=(m1-m2)'*inv(s)*(m1+m2)/2;
for i=1:n1+n2
   if( d'*x(:,i)>c)
       pop(i)=1;
   else
       pop(i)=2;
       
   end
end
n2_1=sum(pop(1:n1)==2);
n1_2=sum(pop(n1+1:n1+n2)==1);
n=n2_1+n1_2;
p2_1=n2_1/n1;
p1_2=n1_2/n2;
fprintf('\n\t coefficint of FDF is:');
disp(d');
fprintf('\n\t constantc=%f\n',c);
fprintf('\n\t n(2/1)=%f\n',n2_1);
fprintf('\n\t n(1/2)=%f\n',n1_2);
fprintf('\n\t total missclassified =%f\n',n);
fprintf('\n\t p(2/1)=%f\n',p2_1);
fprintf('\n\t p(1/2)=%f\n',p1_2);
fprintf('\n\t missclassification prob=%f\n',n/(n1+n2));
%--Classification of x0---------------------------------------------------
if (d'*x0>c)
    fprintf('\n\t Given individual belogs to Population 1');
else
    fprintf('\n\t Given individual belogs to Population 2');
end
%--------------------------------------------------------
%--------------------------------------------------------
%------MULTIPLE AND PARTIAL CORRELATION------------------
clc;clear all;close all;
p=5;n=300;q=1;
s=[ 4.308 1.683 1.803 2.155 -0.253
    1.683 1.768 0.588 0.177 0.176
    1.803 0.588 0.81 1.065 -0.158
    2.155 0.177 1.065 1.970 -0.357
    -0.253 0.176 -0.158 -0.357 0.504]
m=zeros(p,1);
%------------manupulation--mcc---------------------
x=mvnrnd(m,s,n);    
s_hat=cov(x)*(n-1)/n;
s11=s(1:q,1:q); s12=s(1:q,q+1:p);
s21=s12';   s22=s(q+1:p,q+1:p);
s11_hat=s_hat(1:q,1:q); s12_hat=s_hat(1:q,q+1:p);
s21_hat=s12';   s22_hat=s_hat(q+1:p,q+1:p);
r1_2345=sqrt(1-(det(s)/(s11*det(s22))));
sr1_2345=sqrt(1-(det(s_hat)/(s11_hat*det(s22_hat))));
fprintf('\n\t         MCC r1.2345=%f \n',r1_2345);
fprintf('\n\t  Sample MCC r1.2345=%f \n',sr1_2345);
%-------------manupulation--pcc-------------------------
q=2;
s11=s(1:q,1:q);         s12=s(1:q,q+1:p);
s21=s12';               s22=s(q+1:p,q+1:p);
s11_2=s11-s12*inv(s22)*s21;
r12_345=s11_2(1,2)/sqrt(s11_2(1,1)*s11_2(2,2));
fprintf('\n\t Partial corr. coeff r12.345=%f \n',r12_345);

s11_h=s_hat(1:q,1:q);         s12_h=s_hat(1:q,q+1:p);
s21_h=s12_h';               s22_h=s_hat(q+1:p,q+1:p);
s11_2_h=s11_h-s12_h*inv(s22_h)*s21_h;
r12_345_h=s11_2_h(1,2)/sqrt(s11_2_h(1,1)*s11_2_h(2,2));
fprintf('\n\t  Sample Partial corr. coeff r12.345=%f \n',r12_345_h);
%-------------------------------------------------------
%-------------------------------------------------------
%------PRINCIPLE COMPONENT --------------------------
clc;clear all;close all;
%load pop_data
x=[5.935	14.2	2.265	2.27	2.91
1.523	13.1	0.597	0.75	2.62
2.599	12.7	1.237	1.11	1.72
4.009	15.2	1.649	0.81	3.02
4.687	14.7	2.312	2.50	2.22
8.044	15.6	3.641	4.51	2.36
2.766	13.3	1.244	1.03	1.97
6.538	17.0	2.618	2.39	1.85
6.451	12.9	3.147	5.52	2.01
3.314	12.2	1.606	2.18	1.82
3.777	13.0	2.119	2.83	1.80
1.530	13.8	0.798	0.84	4.25
2.768	13.6	1.336	1.75	2.64
6.585	14.9	2.763	1.91	3.17
]';
[p n]=size(x);
%-----------Manupulation-------------------------
s=cov(x');      
R=(diag(diag(s))^-0.5)*s*(diag(diag(s))^-0.5);
[P L]=eig(R);   l=diag(L);
P=P(:,[1 2 5 4 3]);
l=l([1 2 5 4 3]);
TSV=sum(l);     pv=l/TSV;       cpv=cumsum(pv*100);
%-----------Output-----------------------------------
fprintf('\n\t Coefficients of %d  Principle Components (in column) are :\n',p);
disp(P);
fprintf('\n\t Variance of PC                :\n'); disp(l');
fprintf('\n\t proportion of Variance of PC  :\n'); disp(pv');
fprintf('\n\t Cumulative percent of Var. of PC:\n'); disp(cpv');
set(gcf,'color',[1 1 1]);   
plot(cpv,'-bx')
%------------------------------------------------------
%------------------------------------------------------
%--------MANOVA----------------------------------------
clc;clear all;close all;
n1=3;n2=2;n3=3;k=3;
x1=[9 6 9
    3 2 2];
x2=[2 0
    0 4];
x3=[3 1 2
    8 9 7];
x=[x1 x2 x3];    n=n1+n2+n3;
m1=mean(x1,2);   s1=cov(x1');   A1=(n1-1)*s1;
m2=mean(x2,2);   s2=cov(x2');   A2=(n2-1)*s2;
m3=mean(x3,2);   s3=cov(x3');   A3=(n3-1)*s3;
m=mean(x,2);     s=cov(x');     A=(n-1)*s;
T=A;             W=A1+A2+A3;    B=T-W;
L=det(W)/det(T);
f_c=(n-k-1)*(1-sqrt(L))/(sqrt(L)*(k-1));
pvalue=1-fcdf(f_c,2*(k-1),2*(n-k-1));
%----------------Ouptut---------------------------
fprintf('\n\t Matrix of Beteen SS and SCP: \n ');disp(B);
fprintf('\t Matrix of Within SS and SCP: \n ');disp(W);
fprintf('\t Matrix of Total  SS and SCP: \n ');disp(T);
fprintf('\n\t Wilks Lambda= %f \n',L);
fprintf('\n\t F-calculated= %f \n',f_c);
fprintf('\n\t p-value     = %f \n',pvalue);
if pvalue<0.05
fprintf('\n\t Decision   : Reject H0:m1=m2=mk (at %f LOS)',0.05); 
else
fprintf('\n\t Decision   : Do not reject H0:m1=m2=mk'); 
end

%-----------------END OF MATALAB PRACTICAL---------------