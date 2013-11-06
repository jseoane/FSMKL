function [kernelcellaux,kerneloptioncellaux,variablecellaux,ds]=CreateKernelListWithVariable(variablecell,dim,kernelcell,kerneloptioncell,threVar,rankedIDX)
% create list of kernels, from original SimpleMKL "CreateKernelListWithVariable.m"
%
% Input
% variablecell	vector with types of feature-based kernel construction 
%		all: include all features in one kernel
%		single: include one kernel per feature
%		random: include random features in each kernel
%		ranked: include combinations of first n ranked features in each kernel
%		rankedoptim: include just the firnst n ranked features per kernel
% 				(use this option when high memory usage)
%
% dim		dimension cell from "prepareMultipleDataSetCV"
% kernelcell	cell of kernel types to include (see svmkernel for more info)
% kerneloptioncell	parameters of each kernel
% threVar	max number of features to be included in a ranked o rankedoptim kernel
% rankedIDX	original index of best ranked features in the original datasets
%
% Output
% kernelcellaux	auxiliar structure to store kernel information
% kerneloptioncellaux	auxiliar structure to store kernel option information
% variablecellaux	auxiliar structure to store kernel features information
% ds		index of the datasets
%
% This file is a modified version of CreateKernellistWithVariable to include rankedfeatures
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license
% This code is based in SimpleMKL by:
% A. Rakotomamonjy, F. Bach, Y. Grandvalet, S. Canu
% SimpleMKL,  Journal of Machine Learning Research, Vol. 9, pp 2491-2521, 2008


j=1;
for i=1:length(variablecell)
    switch variablecell{i}
        case 'all'
            kernelcellaux{j}=kernelcell{i};
            kerneloptioncellaux{j}=kerneloptioncell{i};
            variablecellaux{j}=1:dim{1};
            ds(j)=1;
            j=j+1;
        case 'single'
            for k=1:dim{1}
                kernelcellaux{j}=kernelcell{i};
                kerneloptioncellaux{j}=kerneloptioncell{i};
                variablecellaux{j}=k;
                ds(j)=1;
                j=j+1;
            end;
        case 'random'
            kernelcellaux{j}=kernelcell{i};
            kerneloptioncellaux{j}=kerneloptioncell{i};
            indicerand=randperm(dim{1});
            ds(j)=1;
            nbvarrand=floor(rand*dim{1})+1;
            variablecellaux{j}=indicerand(1:nbvarrand);
            j=j+1;
            
        case 'ranked'
            numDS = length(rankedIDX);
            for(l=1:numDS)
                if(dim{l}==1)
                    kernelcellaux{j}=kernelcell{i};
                    kerneloptioncellaux{j}=kerneloptioncell{i};
                    variablecellaux{j}=rankedIDX{l};
                    ds(j)=l;
                    j = j+1;
                else
                    if(dim{l}<threVar(l))
                        lim = dim{l};
                    else
                        lim=threVar(l);
                    end
                    for(k=2:lim)
                        kernelcellaux{j}=kernelcell{i};
                        kerneloptioncellaux{j}=kerneloptioncell{i};
                        
                        ds(j) = l;
                        %ds{j}=l;
                        if(k<=dim{l})
                            variablecellaux{j}=rankedIDX{l}(1:k);
                            
                        else
                            variablecellaux{j}=rankedIDX{l}(1:dim{l});
                            
                        end
                        j=j+1;
                    end
                end
            end
            
        case 'rankedoptim'
            numDS = length(rankedIDX);
            for(l=1:numDS)
                if(dim{l}==1)
                    kernelcellaux{j}=kernelcell{i};
                    kerneloptioncellaux{j}=kerneloptioncell{i};
                    variablecellaux{j}=rankedIDX{l};
                    ds(j)=l;
                    j = j+1;
                else
                    if(dim{l}<threVar(l))
                        lim = dim{l};
                    else
                        lim=threVar(l);
                    end
                    if(lim>0)
                        kernelcellaux{j}=kernelcell{i};
                        kerneloptioncellaux{j}=kerneloptioncell{i};
                        ds(j) = l;
                        variablecellaux{j}=rankedIDX{l}(1:lim);
                        j=j+1;
                    end
                    
                end
            end
            
            
    end

end;
