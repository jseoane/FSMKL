function [idxpos] = samplingEqual(yapp)
%% this function performs a sub sampling of the elements of the higher set in unbalanced datasets
% Input
% yapp	Training input vector
%
% Output
% index of new positions
%
%
% FSMKL
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license

	n = length(yapp);
	posidx = find(yapp==1);
	negidx = find(yapp==-1);
	posidxl = length(posidx);
	negidxl = length(negidx);

	if(posidxl >= negidxl)
	  posidx2 = posidx(unidrnd(posidxl,1,length(negidx)));
	  idxpos= [posidx2 ;negidx];
	else          
	  negidx2 = negidx(unidrnd(negidxl,1,length(posidx)));
	  idxpos = [posidx ;negidx2];
	end
end




