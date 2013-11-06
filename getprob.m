function [ps ] = getprob(ypred,A, B)
% Return a probability meassure following Platt approach
% Input
% ypred	output vector of the SVM predicted values
% A	adjusted curve parameter
% B	adjusted vurve parameter
%
% Output
% ps	probability of positive output
%
% FS-MKL
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license



        ps = 1./(1+(exp(A.*ypred+B)));
        
        
end

