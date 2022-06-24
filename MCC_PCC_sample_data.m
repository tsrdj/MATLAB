clc;clear all;
x=[65 72 54 68 55 59 78 58 57 51
56 58 48 61 50 51 55 48 52 42
9 11 8 13 10 8 11 10 11 7];
[p n]=size(x);
s=cov(x)*(n-1)/n;  % MLE of S
%----------Calculation and output-------
sr1_23=sqrt(1-(det(s)/(s(1,1)*det(s(2:p,2:p)))));     
fprintf('\n\t  Sample MCC r1.23=%f \n',sr1_23);
% %----------PCC Calculation and output-------
q=2;
s11=s(1:q,1:q);         s12=s(1:q,q+1:p);
s21=s12';               s22=s(q+1:p,q+1:p);
s11_2=s11-s12*inv(s22)*s21;
sr12_3=s11_2(1,2)/sqrt(s11_2(1,1)*s11_2(2,2));
fprintf('\n\t  Sample Partial corr. coeff r12.3=%f \n',sr12_3);



% y=[65 72 54 68 55 59 78 58 57 51
% 56 58 48 61 50 51 55 48 52 42
% 9 11 8 13 10 8 11 10 11 7];