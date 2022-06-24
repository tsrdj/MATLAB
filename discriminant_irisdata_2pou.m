%-------------------------------------------------
% Discriminant analysis on iris data
%-------------------------------------------------
clc; clear all;
%-----------Iris data--------------------------------
load iris
x1=c(:,1:50);       x2=c(:,51:100);
n1=50;              n2=50; 
%-----------Manupulation-------------------------
m1=mean(x1,2);  s1=cov(x1');
m2=mean(x2,2);  s2=cov(x2');
A=(n1-1)*s1+(n2-1)*s2;  s=A/(n1+n2-2);
%-----------FDF-----------------------------------
d=inv(s)*(m1-m2);     k=d'*(m1+m2)/2;
%-----------Classification among P1 and P2---------
x=[x1 x2];
for i=1:n1+n2;
    if d'*x(:,i)>k
         pop(i)=1;
     else 
         pop(i)=2;
     end
 end
 n2_1=sum(pop(1:n1)==2);     n1_2=sum(pop(n1+1:n1+n2)==1);
n=n2_1+n1_2;
p2_1=n2_1/n1;               p1_2=n1_2/n2;
fprintf('\n\t Classification in P1 and P2 :\n')
fprintf('\n\t Coefficient of FDF is :\n');       disp(d');
fprintf('\n\t Constant c=%f \n',k);
fprintf('\n\t n(2|1)=%f \n',n2_1)
fprintf('\n\t n(1|2)=%f \n',n1_2)
fprintf('\n\t Total misclassified=%d \n',n)
fprintf('\n\t P(2|1)=%f \n',p2_1)
fprintf('\n\t P(1|2)=%f \n',p1_2)

