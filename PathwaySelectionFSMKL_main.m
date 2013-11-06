% Feature Selection Multiple Kernel Learning Algorithm for
% Pathway selection 
% CNV + EXP + pathw + clinic + ER = surv
% jseoane
% j.seoane@bristol.ac.uk
% This code is protected under GPL license
% This code is based in SimpleMKL by:
% A. Rakotomamonjy, F. Bach, Y. Grandvalet, S. Canu
% SimpleMKL,  Journal of Machine Learning Research, Vol. 9, pp 2491-2521, 2008


% load data 
% this algoritm was developed originally to be used with METABRIC dataset
% this example load artifically generated exp and cnv dataset

load('randomData.mat')

l_exp = size(exp,2);
l_cnv = size(cnv,2);


% MKL/SVM output must be in terms of 1/-1
y = output;
y1 = logical(y);
y = -1*ones(length(y1),1);
y(y1)= 1;


cnv_names = strrep(cnv_names(1:l_cnv),'"X','"');
cnv_names = strrep(cnv_names,'_eg','');

%load pathway mapp from KEGG
load('pathways.mat');

genelist_cnv = unique(pathway_cnv(:,1));
%------------------------------------------------------
% For each dataset, pathways and features need to be related (if there are some relation)
%------------------------------------------------------

% all data are stored grouped in pathways in cell "x", and related information are stored in "infok"

% CNV
% generate different datasets for each pathway in CNV

% first found which genes in CNV are related with a pathway and which not.
% note that most of genes are in more than one pathway
no_path_idx = find(strcmp(pathway_cnv(:,2),'"nil"'));	% genes that are not in any pathway
path_idx= find(~strcmp(pathway_cnv(:,2),'"nil"'));	% genes that are in, at least, one pathway
path = pathway_cnv(path_idx,2);				% pathways in this CNV dataset
uni_path = unique(path);				


% this parameters controls the number of genes that are not included in any
% pathway that are included in the FSMKL
max_unique_kernels = 25;

% in this loop we build a cell dataset with information about pathways in CNV (x and infok)
idx = 1;	% this index indicate the number of total datasets, grouped by pathways
%foreach pathway, cell dataset with samples
for(i=1:length(uni_path))
    pathname = uni_path(i);
    pathname_idx = find(strcmp(pathway_cnv(:,2),pathname));
    genesPath = pathway_cnv(pathname_idx,1);
    [int, ix1,ix2]=intersect(genesPath,cnv_names');   
    x{i}=cnv(:,ix2);			%save data in x
    infok{i}.pathway = pathname;	%save pathway name
    infok{i}.genes = genesPath;		%save gene names in this pathway
    infok{i}.indx = ix2; 		%save index of elements of this pathway
    infok{i}.dsor = [1,l_cnv];		%save position of CNV dataset

end

idx = idx + length(uni_path)-1;

% same that in previous loop but with genes that ARE NOT in any pathway

%rank most important CNV features outside pathway
genes_nopath = pathway_cnv(no_path_idx,1);
[int, ix1,ix2]=intersect(genes_nopath,cnv_names');   
cnv_nopath = cnv(:,ix2);
gene_name_nopath =cnv_names(ix2);

% rankfeatures from the bioinformatics toolbox
[rankedIDx] = rankfeatures(cnv_nopath',y);

if(length(no_path_idx)>max_unique_kernels)
    unique_kern = max_unique_kernels;
else
    unique_kern = length(no_path_idx);
end

for(i=1:unique_kern)
    x{i+idx}=cnv_nopath(:,rankedIDx(i));	% save data in x
    infok{i+idx}.pathway = 'nopathway';		% save pathway name (no pathway)
    infok{i+idx}.genes = gene_name_nopath(rankedIDx(i));	% save gene names
    infok{i+idx}.indx = rankedIDx(i); %ix2;			% save index of elements 
    infok{i+idx}.dsor = [1,l_cnv];		% save position of CNV dataset
end

idx = idx + unique_kern;



% EXP
%generate different datasets for each pathway in EXP
% first found which genes in EXP are related with a pathway and which not.
% note that most of genes are in more than one pathway
no_path_idx = find(strcmp(pathway_exp(:,2),'"nil"'));	% genes that are not in any pathway
path_idx= find(~strcmp(pathway_exp(:,2),'"nil"'));	% genes that are in, at least, one pathway
path = pathway_exp(path_idx,2);				% pathways in this EXP dataset
uni_path = unique(path);

%foreach pathway, cell dataset with samples
% in this loop we build a cell dataset with information about pathways in EXP (x and infok)
for(i=1:length(uni_path))
    pathname = uni_path(i);
    pathname_idx = find(strcmp(pathway_exp(:,2),pathname));
    genesPath = pathway_exp(pathname_idx,3);
    [int, ix1,ix2]=intersect(genesPath,exp_names');
    x{idx+i}=exp(:,ix2);			% save data in x
    infok{i+idx}.pathway = pathname;		% save pathway name
    infok{i+idx}.genes = genesPath;		% save gene names
    infok{i+idx}.indx = ix2;			% save index of features in this pathway
    infok{i+idx}.dsor = [l_cnv+1,l_cnv+l_exp];	% save position of EXP dataset
end
idx = idx + length(uni_path)-1;

% same that in previous loop but with genes that ARE NOT in any pathway

%rank most important EXP features outside pathway
genes_nopath = pathway_exp(no_path_idx,3);
[int,ix1,ix2]=intersect(genes_nopath,exp_names');
exp_nopath = exp(:,ix2);
gene_name_nopath =exp_names(ix2);
[rankedIDx] = rankfeatures(exp_nopath',y);

if(length(no_path_idx)>max_unique_kernels)
    unique_kern = max_unique_kernels;
else
    unique_kern = length(no_path_idx);
end

for(i=1:unique_kern)   
    x{i+idx}=exp_nopath(:,rankedIDx(i));	% save data in x
    infok{i+idx}.pathway = 'nopathway';		% save pathway name (nopathway)
    infok{i+idx}.genes = gene_name_nopath(rankedIDx(i));	% save gene names
    infok{i+idx}.indx = rankedIDx(i); %ix2;   	% save index of features
	 infok{i+idx}.dsor = [l_cnv+1,l_cnv+l_exp];	% save position of exp dataset
end
idx = idx + unique_kern;

indexDsor = l_cnv+l_exp +1;

numericClinicVar = numericNewClinic;

% we consider each diferent clinic variable as a new dataset

% Clinical variables
% numeric

% age
x{idx+1} = cell2mat(numericClinicVar(:,1));	% save data in x
infok{idx+1}.pathway = 'clinicalVar.age';	% save pathway name "age"
infok{1+idx}.genes ='age';			% save gene names "age"
infok{1+idx}.indx = 1;				% save index of feature (just one feature)
infok{1+idx}.dsor = [indexDsor,indexDsor];	% save position of age dataset
idx = idx +1;

indexDsor = indexDsor+1;

% size tumour
x{idx+1} = cell2mat(numericClinicVar(:,2));	% save data in x
infok{idx+1}.pathway = 'clinicalVar.size';	% save pathway name "size"
infok{1+idx}.genes ='size';			% save gene names "size"
infok{1+idx}.indx = 1;				% save index of feature (just one feature)
infok{1+idx}.dsor = [indexDsor,indexDsor];	% save position of "size" dataset
idx = idx +1;

indexDsor = indexDsor+1;


% NPI
x{idx+1} = cell2mat(numericClinicVar(:,3));	% save data in x
infok{idx+1}.pathway = 'clinicalVar.npi';	% save pathway name "NPI"
infok{1+idx}.genes ='npi';			% save gene names "NPI"
infok{1+idx}.indx = 1;				% save index of feature (just one feature)
infok{1+idx}.dsor = [indexDsor,indexDsor];	% save position of NPI dataset
idx = idx +1;

indexDsor = indexDsor+1;



clinicNumNames = {'age' 'size' 'npi'};

%categorical
% categorical variables are stored as different dataset, but in binary format
% the categorical variables has been binarized as dummy variable
binaryClinicVar = binaryNewClinic;
% group
x{idx+1} = binaryClinicVar(:,1:5);		% save data in x
infok{idx+1}.pathway = 'clinicalVar.group';	% save pathway name "group"
infok{1+idx}.genes ={'"1"'    '"2"'    '"3"'    '"4"'    '"other"'};	% save gene names "group"
infok{1+idx}.indx = 1:5;			% save index of feature (1:5)
infok{1+idx}.dsor = [indexDsor,indexDsor+4];	% save position of "group" dataset
idx = idx +1;

indexDsor = indexDsor +5;

% grade
x{idx+1} = binaryClinicVar(:,6:9);		% save data in x
infok{idx+1}.pathway = 'clinicalVar.grade';	% save pathway name "grade"
infok{1+idx}.genes ={ '"1"'    '"2"'    '"3"'    '"null"'};	% save gene names "grade"
infok{1+idx}.indx = 1:4;			% save index of feature (1:4)
infok{1+idx}.dsor = [indexDsor,indexDsor+3];	% save position of "grade" dataset
idx = idx +1;

indexDsor = indexDsor +4;

% stage
x{idx+1} = binaryClinicVar(:,10:15);		% save data in x
infok{idx+1}.pathway = 'clinicalVar.stage';	% save pathway name "stage"
infok{1+idx}.genes = { '"0"'    '"1"'    '"2"'    '"3"'    '"4"'    '"null"'};	% save gene names "stage"
infok{1+idx}.indx = 1:6;			% save index of feature (1:6)
infok{1+idx}.dsor = [indexDsor,indexDsor+5];	% save position of "stage" dataset
idx = idx +1;

indexDsor = indexDsor +6;
% histological type
x{idx+1} = binaryClinicVar(:,16:28);		% save data in x
infok{idx+1}.pathway = 'clinicalVar.histo_type';% save pathway name "hist type"
infok{1+idx}.genes ={'BENIGN' 'DCIS' 'IDC' 'IDC+ILC' 'IDC-MED' 'IDC-MUC' 'IDC-TUB' 'ILC' 'INVASIVE TUMOUR' 'MIXED NST AND A SPECIAL TYPE' 'OTHER INVASIVE' 'OTHER' 'PHYL'}; % save gene names "hist type"
infok{1+idx}.indx = 1:13;			% save index of feature (1:13)
infok{1+idx}.dsor = [indexDsor,indexDsor+12];	% save position of "hist type" dataset
idx = idx +1;
indexDsor = indexDsor +13;

% HER SNP 6
x{idx+1} = binaryClinicVar(:,29:32);		% save data in x
infok{idx+1}.pathway = 'clinicalVar.HER_SNP';	% save pathway name "HER SNP 6"
infok{1+idx}.genes ={ '"GAIN"' '"LOSS"'    '"NEUT"'    '"UNDEF"'}; % save gene names "HER SNP 6"
infok{1+idx}.indx = 1:4;			% save index of feature (1:4)
infok{1+idx}.dsor = [indexDsor,indexDsor+3];	% save position of "HER SNP 6" dataset
idx = idx +1;

indexDsor = indexDsor +4;

% cellularity
x{idx+1} = binaryClinicVar(:,33:36);		% save data in x
infok{idx+1}.pathway = 'clinicalVar.cellularity';	% save pathway name "cellularity"
infok{1+idx}.genes ={'"high"' '"low"' '"moderate"' '"undef"'};% save gene names "cellularity"
infok{1+idx}.indx = 1:4;			% save index of feature (1:4)
infok{1+idx}.dsor = [indexDsor,indexDsor+3];	% save position of "cellularity" dataset
idx = idx +1;

indexDsor = indexDsor +4;


% PAM 50 subtype
x{idx+1} = binaryClinicVar(:,37:42);		% save data in x
infok{idx+1}.pathway = 'clinicalVar.pam50';	% save pathway name "PAM 50"
infok{1+idx}.genes = {'"Basal"' '"Her2"' '"LumA"' '"LumB"' '"NC"'  '"Normal"'};% save gene names "PAM 50"
infok{1+idx}.indx = 1:6;			% save index of feature (1:6)
infok{1+idx}.dsor = [indexDsor,indexDsor+5];	% save position of "PAM 50" dataset
idx = idx +1;

indexDsor = indexDsor +6;


% ER HIS
x{idx+1} = ERstatus5(:);			% save data in x
infok{idx+1}.pathway = 'clinicalVar.ER';	% save pathway name "ER"
infok{1+idx}.genes = {'ER+'};			% save gene names "ER"
infok{1+idx}.indx = 1;				% save index of feature (1)
infok{1+idx}.dsor = [indexDsor,indexDsor];	% save position of "ER" dataset
idx = idx +1;

indexDsor = indexDsor +1;

clinicCatNames = {'"1"'    '"2"'    '"3"'    '"4"'    '"other"' '"1"'    '"2"'    '"3"'    '"null"' '"0"'    '"1"'    '"2"'    '"3"'    '"4"'    '"null"' 'BENIGN' 'DCIS' 'IDC' 'IDC+ILC' 'IDC-MED' 'IDC-MUC' 'IDC-TUB' 'ILC' 'INVASIVE TUMOUR' 'MIXED NST AND A SPECIAL TYPE' 'OTHER INVASIVE' 'OTHER' 'PHYL' '"GAIN"' '"LOSS"'    '"NEUT"'    '"UNDEF"' '"high"' '"low"' '"moderate"' '"undef"' '"Basal"' '"Her2"' '"LumA"' '"LumB"' '"NC"'  '"Normal"' 'ER+'};


% max number of genes for pathway used
threshold = 15;    
threeVar = threshold*ones(1,length(x));



%%algorithm parametrization

nbiter=5; %store folds in CV
ratio=(1-nbiter/nbiter);

C = [300];
verbose=1;

options.algo='svmclass'; % Choice of algorithm in mklsvm can be either
                         % 'svmclass' or 'svmreg'
%------------------------------------------------------
% choosing the stopping criterion
%------------------------------------------------------
options.stopvariation=0; % use variation of weights for stopping criterion 
options.stopKKT=0;       % set to 1 if you use KKTcondition for stopping criterion    
options.stopdualitygap=1; % set to 1 for using duality gap for stopping criterion

%------------------------------------------------------
% choosing the stopping criterion value
%------------------------------------------------------
options.seuildiffsigma=1e-2;        % stopping criterion for weight variation 
options.seuildiffconstraint=0.1;    % stopping criterion for KKT
options.seuildualitygap=0.01;       % stopping criterion for duality gap

%------------------------------------------------------
% Setting some numerical parameters 
%------------------------------------------------------
options.goldensearch_deltmax=1e-1; % initial precision of golden section search
options.numericalprecision=1e-8;   % numerical precision weights below this value
                                   % are set to zero 
options.lambdareg = 1e-8;          % ridge added to kernel matrix 

%------------------------------------------------------
% some algorithms paramaters
%------------------------------------------------------
options.firstbasevariable='first'; % tie breaking method for choosing the base 
                                   % variable in the reduced gradient method 
options.nbitermax=500;             % maximal number of iteration
options.seuil=0;                   % forcing to zero weights lower than this 
options.seuilitermax=10;           % vaILMN_1815666lue, for iterations lower than this one 

options.miniter=0;                 % minimal number of iterations 
options.verbosesvm=0;              % verbosity of inner svm algorithm 
options.efficientkernel=1;         % use efficient storage of kernels 


classcode=[1 -1];

[nbdata,dim]=size(x{1});

nbtrain=floor(nbdata*ratio);
%------------------------------------------------------------------------
%                   Building the kernels parameters
%------------------------------------------------------------------------
kernelt={'gaussian' 'poly'  };

kerneloptionvect={[0.5 1 2 ]  [1 2] };

variablevec={ 'rankedoptim' 'rankedoptim' }; %NEW ranked, rankedoptim


importance = zeros(1,79848);	%store importance (kernel parameter) for each kernel
votes = zeros(1,79848);		%store the number of times each kernel is selected during CV
testR = zeros(1,nbiter);	%store test accuracy
trainR = zeros(1,nbiter);	%store train accuracy
probTot = zeros(nbdata,1);	%store probability of each sample
predictions = zeros(nbdata,1);	%store predicion for each sample


% build cross validation folds
c = cvpartition(nbdata,  'kfold' , nbiter); % return the indexes on each fold

%for each CV interation
for i=1: nbiter
    
   %prepare data of each dataset in folds for CV  
   [xapp,xtest,yapp,ytest,IDX,dim,offst] = prepareMultipleDataSetCV(x,y,nbtrain,classcode,i,c); 
    
   % create kernels
    [kernel,kerneloptionvec,variableveccell,ds]=CreateKernelListWithVariable(variablevec,dim,kernelt,kerneloptionvect,threeVar,IDX );
    
    %normalize kernels
    [Weight,InfoKernel,ds2]=UnitTraceNormalization(xapp,kernel,kerneloptionvec,variableveccell,ds);
 
    if(i==1)
	importance = zeros(1,length(Weight));
	votes = zeros(1,length(Weight));
    end
	
    % create big MKL kernel (high memory consumtion depending of sample size 
    % and the number of datasets and pathways used)
    % final kernel size is samples*samples*kernels (# kernesl=length*Weigth)
    K=mklkernelEff(xapp,InfoKernel,Weight,options);

  
    
    
    %MKL minimization
    [beta,w,b,posw,story(i),obj(i)] = mklsvm(K,yapp,C,options,verbose);

    %eval train error
    ktrain = mklkernel(xapp,InfoKernel,Weight,options,xapp(posw,:),beta);
    ypredt = ktrain*w+b;

    % using a modified ROC function from G. Cardillo
    ROCout= roc3([(ypredt+1)/2, (yapp+1)/2],[],0.05,0);
    trainR(i) = ROCout.perf; % store train performance

    % sampling function to calculate probability measure parameters
    ibal = samplingEqual(yapp);
    yapp2 = yapp(ibal);
    ypred2 = ypredt(ibal);

    [prob A B] = procesprob(yapp2,ypred2);



    %eval test error
    Kt=mklkernel(xtest,InfoKernel,Weight,options,xapp(posw,:),beta);
    ypred=Kt*w+b;
    ROCout2= roc3([(ypred+1)/2, (ytest+1)/2],[],0.05,0);


    testR(i)=ROCout2.perf;  % store test performance
    pps(i)=ROCout2.pps; % store number of correct matches
    ppn(i)=ROCout2.ppn; % store number of incorrect matches
    auroc(i)=ROCout2.AUC; % store AUC
    cm{i} = ROCout2.cotable; % store confusion matrix


    probTest = getprob(ypred,A,B);


    indice1.test = find(test(c,i));
    probTot(indice1.test) = probTest;
    predictions(indice1.test)=ypred;


   % pvalue for probability measure
   threP = 0.05;
   posC = probTest>=(1-threP);
   negC = probTest<=(threP);
   
   indxConf = or(posC, negC);



   ROCout3= roc3([(ypred(indxConf)+1)/2, (ytest(indxConf)+1)/2],[],0.05,0)
   testRProb(i) =ROCout3.perf;
   aurocProb(i) =ROCout3.AUC;
   cmProb{i} = ROCout3.cotable;

   
    importance = importance + beta;
    votes = votes + (beta>0);
clear('K');

end;%
%------------------------------------------------------
%store pathway information
%------------------------------------------------------
a = GetKernelInformationPathway2(beta,InfoKernel,ds2,infok,200, importance,votes)
filesalida = strcat('outputPath_','clini_surv_005','.xls'); %output file 1
SaveInformationPathway2(a,filesalida,[cnv_names exp_names clinicNumNames clinicCatNames],offst);

%------------------------------------------------------
%remove outliers with prob measure less than 0.8
%------------------------------------------------------
outliers = find((probTot>0.2)&(probTot<0.8));
nooutliers = find((probTot<=0.2)|(probTot>=0.8));


for i = 1:length(x) 
 x{i}(outliers,:) = [];
end;

y(outliers,:)=[];



[nbdata,dim]=size(x{1});

nbtrain=floor(nbdata*ratio);


importance2 = zeros(1,230);
votes2 = zeros(1,230);
testR2 = zeros(1,nbiter);
trainR2 = zeros(1,nbiter);
%------------------------------------------------------
%% second CV and MKL minization, this time without outliers
%------------------------------------------------------

c = cvpartition(nbdata,  'kfold' , nbiter); % return the indexes on each fold

for i=1: nbiter
    
       
    % prepare data in folds for CV   
   [xapp,xtest,yapp,ytest,IDX,dim,offst] = prepareMultipleDataSetCV(x,y,nbtrain,classcode,i,c); 
    
     
   %create kernels
    [kernel,kerneloptionvec,variableveccell,ds]=CreateKernelListWithVariable(variablevec,dim,kernelt,kerneloptionvect,threeVar,IDX );

    % normalize kernels
    [Weight,InfoKernel,ds2]=UnitTraceNormalization(xapp,kernel,kerneloptionvec,variableveccell,ds);
 
    if(i==1)
	importance = zeros(1,length(Weight));
	votes = zeros(1,length(Weight));
    end

	% create big kernel	
    K=mklkernelEff(xapp,InfoKernel,Weight,options);


    
  
    
   % mkl minimization
    [beta,w,b,posw,story(i),obj(i)] = mklsvm(K,yapp,C,options,verbose);



%eval train error
    ktrain = mklkernel(xapp,InfoKernel,Weight,options,xapp(posw,:),beta);
    ypredt = ktrain*w+b;
   ROCout= roc3([(ypredt+1)/2, (yapp+1)/2],[],0.05,0);
    trainR_2(i) = ROCout.perf;


    ibal = samplingEqual(yapp);
    yapp2 = yapp(ibal);
    ypred2 = ypredt(ibal);

     [prob A B] = procesprob(yapp2,ypred2);




%eval test error
    Kt=mklkernel(xtest,InfoKernel,Weight,options,xapp(posw,:),beta);
    ypred=Kt*w+b;

    ROCout2= roc3([(ypred+1)/2, (ytest+1)/2],[],0.05,0);
testR_2(i)=ROCout2.perf;
pps_2(i)=ROCout2.pps;
ppn_2(i)=ROCout2.ppn;
auroc_2(i)=ROCout2.AUC;
cm_2{i} = ROCout2.cotable;


		%calculate probability
    probTest = getprob(ypred,A,B);


     indice1.test = find(test(c,i));
     probTot(indice1.test) = probTest;
     predictions(indice1.test)=ypred;



   threP = 0.05;
   posC = probTest>=(1-threP);
   negC = probTest<=(threP);
   
   indxConf = or(posC, negC);


ROCout3= roc3([(ypred(indxConf)+1)/2, (ytest(indxConf)+1)/2],[],0.05,0);
testRProb_2(i) =ROCout3.perf;
aurocProb_2(i) =ROCout3.AUC;
cmProb_2{i} = ROCout3.cotable;

    
    importance = importance + beta;
    votes = votes + (beta>0);
clear('K');

end;%


%------------------------------------------------------   
%store pathway information                                
%------------------------------------------------------   
a = GetKernelInformationPathway2(beta,InfoKernel,ds2,infok,200, importance,votes)
filesalida = strcat('outputPathRec_','clini_surv_005','.xls'); %output file pathways without outliers
SaveInformationPathway2(a,filesalida,[cnv_names exp_names clinicNumNames clinicCatNames],offst);


%------------------------------------------------------   
% print results                                 
%------------------------------------------------------ 

mean(trainR) %mean train error of CV
std(trainR) %std train error of CV
testR				%test error
mean(testR)
std(testR)
pps		%number of correct evaluated examples for each CV fold
ppn   %number of failed evaluated examples for each CV fold
auroc  % test auroc
mean(auroc)
std(auroc)
cm{:}	% test confusion matrix for each CV fold

testRProb			% test error with probability measure
mean(testRProb)
std(testRProb)
aurocProb			% test auroc with probability measure
mean(aurocProb)
std(aurocProb)
cmProb{:}			% test confusion matrix for each CV fold with probability measure

testR_2			% test error without outliers
mean(testR_2)
std(testR_2)
auroc_2			% test auroc without outliers
mean(auroc_2)
std(auroc_2)
cm_2{:}			% test confusion matrix for each CV fold without outliers

testRProb_2		% test error without outliers with probability measure
mean(testRProb_2)
std(testRProb_2)
aurocProb_2		% test auroc without outliers with probability measure
mean(aurocProb_2)
std(aurocProb_2)
cmProb_2{:}		% test confusion matrix without outliers with probability measure





% THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR 
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
% IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
