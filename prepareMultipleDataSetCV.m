function [xapp,xtest,yapp,ytest,IDX,dim,offst] = prepareMultipleDataSet(x,y,nbtrain,classcode, fold,c)
% Build cross validation dataset for each dataset
%
% Input
% x 		cell with one dataset in each cell
% y 		output
% nbtrain 	number of samples for training
% classcode 	codification of class ouput [-1 1]
% fold 		fold of the nfold cross validation to be used
% c 		cross validation object from "cvpartition"
% 
% Output
% xapp		Train input data
% xtest		Test input data
% yapp		Train output data
% ytest		Test ouput data
% IDX		original index of best ranked features in the original datasets
% dim		cell with dimension of each dataset
% offst		offset of feature in each dataset
%
% FSMKL
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license

numdatasets = length(x);

offset = 0;
xapp = [];
xtest = [];
for(i=1:numdatasets)
    xi = x{i};

    if(i==1)
	indice1.app = find(training(c,fold) );
	indice1.test = find(test(c,fold) );
 	xapp1 = xi(indice1.app,:);
	yapp1 = y(indice1.app);
	xtest1 = xi(indice1.test,:);
	ytest1 = y(indice1.test);
        
        yapp = yapp1;
        ytest = ytest1;
    else
        xapp1 = x{i}(indice1.app,:);
        yapp1 = y(indice1.app);
        xtest1 = x{i}(indice1.test,:);
        ytest1 = y(indice1.test);
    end
            
        
    [xapp1,xtest1]=normalizemeanstd(xapp1,xtest1);
    dim{i} = size(xapp1,2);
    X{i}=xapp1;
    Xt{i}=xtest1;
    
    %%feature selection over each individual dataset
    %for each datasource, make different rank feature
    [IDX1, Z1] = rankfeatures(xapp1', yapp1);
    IDX{i}=IDX1+offset;
    offst{i}=offset;
    offset = offset + size(xapp1,2);
    
    xapp = [xapp xapp1];
    xtest = [xtest xtest1];
end



