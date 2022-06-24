%-----------------------------------------------------
% Testing H0:m1=m2 (Two sample problem)
% One Sample data given in X matrix of size p x n1 
% Second Sample data given in Y matrix of size p X n2
%------------------------input-------------------------
clc;clear all;      load data_iris
[p n]=size(x);      alpha=0.05;
n1=50;n2=50;        n=n1+n2;
%---------------------manipulation and calculation--------------------
X=x(:,51:100);            Y=x(:,101:150);
xbar=mean(X,2);         ybar=mean(Y,2);   
s1=cov(X') ;            s2=cov(Y') ;        
A=(n1-1)*s1+(n2-1)*s2;  s_pooled=A/(n-2);
D2=(xbar-ybar)'*inv(s_pooled)*(xbar-ybar);
T2=n1*n2*D2/n;          F_c=(n-p-1)*T2/((n-2)*p);
p_value=1-fcdf(p,n-p-1,F_c);
%--------------output------------------------
fprintf('\n\t MSD D2        =%f',D2);
fprintf('\n\t Hotellings T2 =%f',T2);
fprintf('\n\t Calculated F  =%f',F_c);
fprintf('\n\t p_value       =%f',p_value);
if (p_value<alpha)
    fprintf('\n\t conclusion: Ho is rejected \n\t')
else
    fprintf('\n\t  conclusion: Fail to rejecte H0 \n\t')
end
