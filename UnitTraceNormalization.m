function [Weigth,InfoKernel,ds2]=UnitTraceNormalization(x,kernelvec,kerneloptionvec,variablevec,ds)
%	Perform kernel normalization (unit trace normalization)
%
% Input
% x		Training input data
% kernelvec	Kernel structure
% kerneloptionvec	Kernel option structure
% variablevec	Features information structure
% ds		Information about the dataset
%
% Output
% Weigth	Kernel weigth
% InfoKernel	Structure with information of kernels
% ds2		Information about the dataset
%
% FSMKL
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license
% This code is based in SimpleMKL by:
% A. Rakotomamonjy, F. Bach, Y. Grandvalet, S. Canu
% SimpleMKL,  Journal of Machine Learning Research, Vol. 9, pp 2491-2521, 2008



chunksize=200;
N=size(x,1);
nbk=1;
for i=1:length(kernelvec);
    % i
    for k=1:length(kerneloptionvec{i})

        somme=0;

        chunks1=ceil(N/chunksize);

        for ch1=1:chunks1
            ind1=(1+(ch1-1)*chunksize) : min( N, ch1*chunksize);
            somme=somme+sum(diag(svmkernel(x(ind1,variablevec{i}),kernelvec{i},kerneloptionvec{i}(k))));
        end;
        %         for j=1:N
        %             somme=somme+svmkernel(x(j,variablevec{i}),kernelvec{i},kerneloptionvec{i}(k));
        %
        %         end
        if somme~=0
            Weigth(nbk)=1/somme;
	    ds2(nbk)=ds(i);
            InfoKernel(nbk).kernel=kernelvec{i};
            InfoKernel(nbk).kerneloption=kerneloptionvec{i}(k);
            InfoKernel(nbk).variable=variablevec{i};
            InfoKernel(nbk).Weigth=1/somme;
            nbk=nbk+1;
%         else
%             A
        end;
    end;
end;
