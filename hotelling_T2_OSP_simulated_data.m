%-----------------------------------------------------
% Testing H0:m=mo (One sample problem) based on simualation
% Sample data is to be simulated in X matrix of size n X p
%------------------------input-------------------------
clc;clear all;close all;
p=5; n=100;
m=[4.32,14.01,1.95,2.17,2.45]';
mo=[5.52,15.47,2.51,2.96,3.30]';
s=[4.308 1.683 1.803 2.155 -0.253;
    1.683 1.768 0.588 0.177 0.176;
    0.803 0.588 0.81 1.065 -0.158;
    2.155 0.177 1.065  1.970 -0.357;
    -0.253 0.176 -0.158 -0.357  0.504];
%---------------------manipulation---------------------
x=mvnrnd(m,s,n)';
m_hat=mean(x,2);    s_hat=cov(x') ;
%--------------------testing of hypothesis-------------
t2=n*(m_hat-mo)'*inv(s_hat)*(m_hat-mo);
F_c=((n-p)/(n-1)*p)*t2;
p_value=1-fcdf(p,n-p,F_c);
%-----------------------Ouptut---------------------------
fprintf('\n\t sample mean:');       disp(m_hat');
fprintf('\n\t sample v_c_m is\n');  disp(s_hat);
fprintf('\n\t The value of hotelling t2 =%f\n',t2);
fprintf('\n\t The value of calculated F =%f\n',F_c)
fprintf('\n\t The p_value for testing Ho:=%f\n',p_value);
if (p_value<0.05)
    fprintf('\n\t conclusion: Ho is rejected \n\t')
else
    fprintf('\n\t conclusion: Fail to reject Ho. \n\t')
end
