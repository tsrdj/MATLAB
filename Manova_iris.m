%------- MANOVA----------------------------------
% Testing H0:m1=m2=...=mk for Iris data
%-------------------------------------------------
clc;clear all;close all;
%----------------input----------------------------
load iris_data
n1=50;   n2=50;   n3=50;   k=3;    a=0.05;  p=4;
x1=x(:,1:n1);
x2=x(:,n1+1:n1+n2);
x3=x(:,n1+n2+1:n1+n2+n3);
%----------------Computational Part---------------
n=n1+n2+n3;
m1=mean(x1,2);   s1=cov(x1');   A1=(n1-1)*s1;
m2=mean(x2,2);   s2=cov(x2');   A2=(n2-1)*s2;
m3=mean(x3,2);   s3=cov(x3');   A3=(n3-1)*s3;
m=mean(x,2);     s=cov(x');     A=(n-1)*s;
T=A;             W=A1+A2+A3;    B=T-W;
L=det(W)/det(T);
f_c=(n-p-2)*(1-sqrt(L))/(sqrt(L)*p);
pvalue=1-fcdf(f_c,2*p,2*(n-p-2));
%----------------Ouptut---------------------------
fprintf('\n\t Matrix of Beteen SS and SCP: \n ');disp(B);
fprintf('\t Matrix of Within SS and SCP: \n ');disp(W);
fprintf('\t Matrix of Total  SS and SCP: \n ');disp(T);
fprintf('\n\t Wilks Lambda= %f \n',L);
fprintf('\n\t F-calculated= %f \n',f_c);
fprintf('\n\t p-value     = %f (at %f LOS) \n',pvalue,a);
if pvalue<a
fprintf('\n\t Decision   : Reject H0:m1=m2=...=mk'); 
else
fprintf('\n\t Decision   : Do not reject H0:m1=m2=...=mk'); 
end
fprintf('\n\t Mean vectors (in column) are:\n');
disp([m1 m2 m3]);






