function Kse = build_efficient_KEFF(Ks);
% build efficient representation of the kernel matrices
% 
% Input
% Ks	Original kernel representation
% 
% Output
% Kse	New efficient kernel codificaation
% efficient_type = 0 -> not efficient
%                  1 -> efficient
%matlab7 =  str2num(version('-release')) >=14;
% modification from build_efficientK.m from SimpleMKL
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license
% This code is based in SimpleMKL by:
% A. Rakotomamonjy, F. Bach, Y. Grandvalet, S. Canu
% SimpleMKL,  Journal of Machine Learning Research, Vol. 9, pp 2491-2521, 2008
%


Kse.nbkernel = size(Ks,1);
Kse.n = size(Ks(1).minik,1);

nbkernel = size(Ks,1);
n = size(Ks(1).minik,1);
if isa(Ks(1).minik,'single');
    Kse.data = zeros(n*(n+1)/2,nbkernel,'single');
else
    Kse.data = zeros(n*(n+1)/2,nbkernel);
end


for j=1:nbkernel
    if isa(Ks(1).minik,'single');
        Kse.data(:,j) = vectorize_single(Ks(j).minik);
    else
        Kse.data(:,j) = vectorize(Ks(j).minik);
    end
end

