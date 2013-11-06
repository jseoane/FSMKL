function K=mklkernelEff(xapp,InfoKernel,Weight,options,xsup,beta)
% Build a disk cached version of the kernel. Use it if your memory needs are high.
%
% Input
% xapp		train input data
% infoKernel	structure with kernel information
% Weigth	original weigth of kernel
% options	option structure for efficient kernel (not used)
% xsup		support vectors
% beta		MKL beta parameter
%
% Output
% K		Kernel
%
%
% FSMKL
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license
% This code is based in SimpleMKL by:
% A. Rakotomamonjy, F. Bach, Y. Grandvalet, S. Canu
% SimpleMKL,  Journal of Machine Learning Research, Vol. 9, pp 2491-2521, 2008




if nargin <5
    xsup=xapp;
    beta=[];

     symbols = ['a':'z' 'A':'Z' '0':'9'];
     MAX_ST_LENGTH = 5;
     stLength = randi(MAX_ST_LENGTH);
     nums = randi(numel(symbols),[1 stLength]);
     filename = strcat('kern_temp_',symbols (nums),'.tmp');
    
    
    lxapp = size(xapp,1);
    lweigh = size(Weight,2);
    %build file
    tic
    fid = fopen(filename,'w');
    
    fwrite(fid, zeros(lweigh,lxapp,lxapp), 'double');
    fclose(fid);
    toc
    keff = memmapfile(filename,'Format',{'double' [lxapp,lxapp] 'minik'}, 'Writable',true);
    tic
    for k=1:length(Weight)

        Kr=svmkernel(xapp(:,InfoKernel(k).variable),InfoKernel(k).kernel,InfoKernel(k).kerneloption, xsup(:,InfoKernel(k).variable));
 
        Kr=Kr*Weight(k);
%         if options.efficientkernel
%             Kr=build_efficientK(Kr);
%         end;

        keff.data(k).minik=Kr;
        

    end;
    toc
    tic
    K = build_efficientKEFF(keff.data);
    clear keff;
    delete(filename);
    toc
else
    ind=find(beta);
    K=zeros(size(xapp,1),size(xsup,1));
    for i=1:length(ind);
        k=ind(i); 
        Kr=svmkernel(xapp(:,InfoKernel(k).variable),InfoKernel(k).kernel,InfoKernel(k).kerneloption, xsup(:,InfoKernel(k).variable));
        Kr=Kr*Weight(k);
        K=K+ Kr*beta(k);
    end;

end;
