function [a]=GetKernelInformationPathway(W,InfoKernel,ds2,infok,n,importance,votes)
% 
% Extract information from selected kernels
%
% Input
% W		kernel parameter
% InfoKernel	auxiliar structure with information about kernels
% ds2		Index of original datasets
% infok		auxiliary structure with information about datasets
% n		Number of kernels to include inforamtion
% importance	Sumatory of kernel parameters
% votes		Number of times some kernel was selected
%
% Output
% a		Structure with information of kernels
%
% FSMKL
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license
% This code is based in SimpleMKL by:
% A. Rakotomamonjy, F. Bach, Y. Grandvalet, S. Canu
% SimpleMKL,  Journal of Machine Learning Research, Vol. 9, pp 2491-2521, 2008




[B,idx] = sort(W,2,'descend');

%for(i=1:length(W))
for(i=1:n)
j = idx(i);

   k = idivide(int16(j),int16(6),'floor')+1;   % 3??? cada dataset tiene 3 kernels? por?
   if(W(j)>0)
       
       fprintf('Kernel %d, sigma=%f, kernel %s,option %f variables %s \n ',j,W(j),InfoKernel(j).kernel,InfoKernel(j).kerneloption, int2str(InfoKernel(j).variable'))
       ds1 = int16(ds2(j));
       ptwy = char(infok{ds1}.pathway);
       fprintf('Dataset %d, pathway %s, variables  \n ',ds1,ptwy)
       a{i}.sigma = W(j);
       a{i}.varIDX = InfoKernel(j).variable;
       a{i}.pathway = ptwy;
       a{i}.genes = infok{ds1}.genes;
       a{i}.indx = infok{ds1}.indx;
       a{i}.kernNum = j;
       a{i}.ds = ds1;
       a{i}.dsor = infok{ds1}.dsor;
       a{i}.importance = importance(j);
       a{i}.votes = votes(j);

       

       %get the names of the selected variables
       %cnv_names(infok{ds1}.indx(InfoKernel(j).variable))
   end
    
    
end
