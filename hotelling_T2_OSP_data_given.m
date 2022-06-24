%-----------------------------------------------------
% Testing H0:m=mo (One sample problem)
% Sample data given in X matrix of size n X p
% TOH Ho:m=m0 where m0=[4
%                       50
%                       10]
%------------------------input-------------------------
clc;clear all;
load data_sweat
[p n]=size(x);
alpha=0.05;
mo=[40 50 10]';
%---------------------manipulation and calculation--------------------
m_hat=mean(x,2);        s_hat=cov(x') ;
t2=n*(m_hat-mo)'*inv(s_hat)*(m_hat-mo);
F_c=(n-p)*t2/((n-1)*p);
p_value=1-fcdf(p,n-p,F_c);
%--------------output------------------------
fprintf('\n\t sample mean:');       disp(m_hat');
fprintf('\n\t sample v_c_m is\n');  disp(s_hat);
fprintf('\n\t The value of hotelling t2 is=%f',t2);
fprintf('\n\t The value of calculated F  is=%f',F_c);
fprintf('\n\t The p_value is=%f',p_value);
if (p_value<alpha)
    fprintf('\n\t conclusion: Ho is rejected \n\t')
else
    fprintf('\n\t  conclusion: Fail to rejecte H0 \n\t')
end
