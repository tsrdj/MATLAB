%------- MANOVA                             -----
% Testing H0:m1=m2=...=mk
%-------------------------------------------------
clc;clear all;close all;
% %----------------Data set 1: input----------------------------
n1=3;   n2=2;   n3=3;   k=3;    a=0.010;
x1=[9 6 9
    3 2 2];
x2=[0 2
    4 0];
x3=[3 1 2
    8 9 7];
%----------------Data set 2: input----------------------------
% n1=5;   n2=3;   n3=4;   k=3;    a=0.010;
% x1=[6 5 8 4 7
%     7 9 6 9 9];
% x2=[3 1 2
%     3 6 3];
% x3=[2 5 3 2
%     3 1 1 3];
%----------------Computational Part---------------

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
if pvalue<a
fprintf('\n\t Decision   : Reject H0:m1=m2=...=mk (at %f LOS)',a); 
else
fprintf('\n\t Decision   : Do not reject H0:m1=m2=...=mk'); 
end






