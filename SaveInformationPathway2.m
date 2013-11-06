function SaveInformationPathway(kernelsInfoA,filename, dsnames,offst)
% store the data related with FSMKL in a excel file
% foreach k, order by weigth, get:
% number dataset
% number of kernel
% weigth
% patway
% kernel-related variables (features selected in Xapp) with offset included
% kernel-related offset
% index of the variables of x{i} in names vector
%
% Input
% kernelsInfoA	kernel info structure from GetKernelInformationPathway
% filename	filename fo the new file
% dsnames	vector with name of datasets
% offst		offset of each feature in each dataset
%
%
% FSMKL
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license

numkernels = length(kernelsInfoA);

fid = fopen(filename,'w');
fprintf(fid,'DS number, Kernel number, W, importance, votes, pathway, genes \n')

for (i = 1:numkernels)
    k = kernelsInfoA{i};
    
    %ds = idivide(int16(k.kerNum),int16(3),'floor')+1;
    
    
    excel.nkern = k.kernNum;
    excel.nds = k.ds;
    excel.w = k.sigma;
    excel.importance = k.importance;
    excel.votes = k.votes;
    excel.pathway = k.pathway;
    
    
    ds_names = dsnames(k.dsor(1):k.dsor(2)); 
    varIDX = k.varIDX-offst{k.ds};
    %idx = k.indx{k.ds}(varIDX);
    idx = k.indx(varIDX);
    excel.genenames = ds_names(idx);
    
    fprintf(fid,'%d,%d,%f,%f,%d,%s, ',excel.nds,excel.nkern, excel.w ,excel.importance,excel.votes, excel.pathway);
    for(j=1:49)
        try
            fprintf(fid,'%s ,', excel.genenames{j});
        catch
            fprintf(fid,'0 ,');
        end
    end
    fprintf(fid,'0 \n');
    
    
end

fclose(fid);
