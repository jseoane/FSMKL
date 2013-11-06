% function minimizing
function [ps A B resnorm] = procesprob(yapp,ypred)
% Calculate the parameters for probability estimation using Platt model
%
% Input
% yapp		Training output
% ypred		predicted train output
%
% Output
% ps		Vector of probabilities
% A		Parameter A
% B		Parameter B
% resnorm	optimization information
%
% FSMKL
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license


yt = (yapp+1)/2;
ypred = (ypred+1)/2;

np = sum(yt>0);
nn = sum(yt<=0);
%x0 = [-0.012 0.001];
x0 = [0 log((np+1)/(nn+1))];

hitarget = (np+1)/(np+2);
lotarget = 1/(nn+2);

ytemp = ones(length(yt),1)*lotarget;

ytix = find(yt==1);
ytemp(ytix) = hitarget;
yt = ytemp;




[x, resnorm]=lsqnonlin(@myfun,x0);


ps = pdi(ypred,x(1),x(2));

A = x(1);
B = x(2);

 %plot(-2:0.1:2,pdi(-2:0.1:2,x(1),x(2)))
 %hold all
 %plot(ypred,pdi(ypred,x(1),x(2)),'*r')



function F = myfun(x)
  
  pi = 1./(1+(exp(x(1).*yt+x(2))));
  t1 = ypred.*log(pi)+(1-ypred).*log(1-pi);  
  F = sum(t1);
end

    function pi = pdi(x,A,B)
        pi = 1./(1+(exp(A.*x+B)));
    end

end
