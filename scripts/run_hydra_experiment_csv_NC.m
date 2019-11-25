function [CIDX,ARI] = run_hydra_experiment_csv_NC(featureCSV,outputDir,varargin)

if nargin==0
    printhelp()
    return
end
    

% function returns estimated subgroups by hydra for clustering
% configurations ranging from K=1 to K=10, or another specified range of
% values. The function returns also the Adjusted Rand Index that was
% calculated across the cross-validation experiments and comparing
% respective clustering solutions.
%
% INPUT
%  
% REQUIRED 
% featureCSV : .csv file containing the input features. (REQUIRED)
%              every column of the file contains values for a feature, with 
%              the exception of the first and last columns. We assume that 
%              the first column contains subject identifying information 
%              while the last column contains label information. First line 
%              of the file should contain header information. Label
%              convention: -1 -> control group - 1 -> pathological group
%              that will be partioned to subgroups 
% outputDir : directory where the output from all folds will be saved (REQUIRED)
%
% OPTIONAL
%
% covCSV : .csv file containing values for different covariates, which
%           will be used to correct the data accordingly (OPTIONAL). Every
%           column of the file contains values for a covariate, with the
%           exception of the first column, which contains subject
%           identifying information. Correction is performed by solving a 
%           solving a least square problem to estimate the respective 
%           coefficients and then removing their effect from the data. The 
%           effect of ALL provided covariates is removed. If no file is 
%           specified, no correction is performed.
%
% NOTE: featureCSV and covCSV files are assumed to have the subjects given
%       in the same order in their rows
%
% C : regularization parameter (positive scalar). smaller values produce 
%     sparser models (OPTIONAL - Default 0.25)
% reg_type : determines regularization type. 1 -> promotes sparsity in the
%            estimated hyperplanes - 2 -> L2 norm (OPTIONAL - Default 1)
% balance : takes into account differences in the number between the two
%           classes. 1-> in case there is mismatch between the number of
%           controls and patient - 0-> otherwise (OPTIONAL - Default 1)
% init : initialization strategy. 0 : assignment by random hyperplanes 
%        (not supported for regression), 1 : pure random assignment, 2: 
%        k-means assignment, 3: assignment by DPP random 
%        hyperplanes (default)
% iter : number of iterations between estimating hyperplanes, and cluster
%        estimation. Default is 50. Increase if algorithms fails to
%        converge
% numconsensus : number of clustering consensus steps. Default is 20.
%                Increase if algorithm gives unstable clustering results.
% kmin : determines the range of clustering solutions to evaluate
%             (i.e., kmin to kmax). Default  value is 1.
% kmax : determines the range of clustering solutions to evaluate
    %             (i.e., kmin to kmax). Default  value is 10.
    % kstep: determines the range of clustering solutions to evaluate
    %             (i.e., kmin to kmax, with step kstep). Default  value is 1.
    % cvfold: number of folds for cross validation. Default value is 10.
    % vo : verbose output (i.e., also saves input data to verify that all were 
    %      read correctly. Default value is 0
    %
    % OUTPUT: 
    % CIDX: sub-clustering assignments of the disease population (positive
    %       class).
    % ARI: adjusted rand index measuring the overlap/reproducibility of
    %      clustering solutions across folds
    %
    % NOTE: to compile this function do
    % mcc -m  run_hydra_experiment_csv 
%
%
% EXAMPLE USE:
% run_hydra_experiment_csv('test.csv','.','covCSV','test_covar.csv','C',0.25,'kmax',3,'init',3,'cvfold',5)

% hydra parameters
% Set as 1 if data is kernelized otherwise 0
% data need to be kernelized for high-dimensional 
params.kernel=0; 

% input parser
p = inputParser;

% required parameters
p.addRequired('featureCSV',@(x)validateattributes(x,{'char'},...
        {'nonempty'}));
p.addRequired('outputDir',@(x)validateattributes(x,{'char'},...
        {'nonempty'}));

% optional parameters
p.addParameter('covCSV',[],@(x)validateattributes(x,{'char'},...
        {'nonempty'}));
p.addParameter('C',0.25); 
p.addParameter('reg_type',1);
p.addParameter('balance',1);
p.addParameter('init',3);
p.addParameter('iter',50);
p.addParameter('numconsensus',20);
p.addParameter('kmin',1);
p.addParameter('kmax',10);
p.addParameter('kstep',1);
p.addParameter('cvfold',10);
p.addParameter('vo',0);

% parse input 
disp('Parsing input...')
parse(p,featureCSV,outputDir,varargin{:});

% create output directory
if (~exist(outputDir,'dir'))
    [status,~,~] = mkdir(outputDir);
    if (status == 0)
        error('run_hydra_experiment_csv:argChk','Cannot create output directory!');
    end
end
    

% assign input parameters to params structure
if(isdeployed)
    if(all(strcmp('C',p.UsingDefaults)==0))
        params.C=str2double(p.Results.C) ; % Regularization parameter
    else
        params.C=p.Results.C ; 
    end
    if(all(strcmp('reg_type',p.UsingDefaults)==0))
        params.reg_type=str2double(p.Results.reg_type); % Regularization type
    else
        params.reg_type=p.Results.reg_type;
    end
    if(all(strcmp('balance',p.UsingDefaults)==0))
        params.balanceclasses=str2double(p.Results.balance); % account for class inbalance
    else
        params.balanceclasses=p.Results.balance;
    end
    if(all(strcmp('init',p.UsingDefaults)==0))
        params.init_type=str2double(p.Results.init); % Algorithm initialization type
    else
        params.init_type=p.Results.init; 
    end
    if(all(strcmp('iter',p.UsingDefaults)==0))
        params.numiter=str2double(p.Results.iter); % number of iterations
    else
        params.numiter=p.Results.iter;
    end
    if(all(strcmp('numconsensus',p.UsingDefaults)==0))
        params.numconsensus=str2double(p.Results.numconsensus); % number of clustering consensus steps
    else
        params.numconsensus=p.Results.numconsensus;
    end
    if(all(strcmp('kmin',p.UsingDefaults)==0))
        params.kmin = str2double(p.Results.kmin) ; % min number of clusters to consider
    else
        params.kmin = p.Results.kmin ; 
    end
    if(all(strcmp('kmax',p.UsingDefaults)==0))
        params.kmax = str2double(p.Results.kmax) ; % max number of clusters to consider
    else
        params.kmax = p.Results.kmax ;
    end
    if(all(strcmp('kstep',p.UsingDefaults)==0))
        params.kstep = str2double(p.Results.kstep) ; % step defining the range clustering solutions to consider
    else
        params.kstep = p.Results.kstep ;
    end
    if(all(strcmp('cvfold',p.UsingDefaults)==0))
        params.cvfold = str2double(p.Results.cvfold) ; % number of folds in k-fold cross validation
    else
        params.cvfold = p.Results.cvfold ; 
    end
    if(all(strcmp('vo',p.UsingDefaults)==0))
        params.vo = str2double(p.Results.vo) ; % verbose output or no
    else
        params.vo = p.Results.vo ; 
    end
else
    params.C=p.Results.C ; % Regularization parameter
    params.reg_type=p.Results.reg_type; % Regularization type
    params.balanceclasses=p.Results.balance; % account for class inbalance
    params.init_type=p.Results.init; % Algorithm initialization type
    params.numiter=p.Results.iter; % number of iterations
    params.numconsensus=p.Results.numconsensus; % number of clustering consensus steps
    params.kmin = p.Results.kmin ; % min number of clusters to consider
    params.kmax = p.Results.kmax ; % max number of clusters to consider
    params.kstep = p.Results.kstep ; % step defining the range clustering solutions to consider
    params.cvfold = p.Results.cvfold ; % number of folds in k-fold cross validation    
    params.vo = p.Results.vo ; % verbose output or no    
end

% confirm validity of optional input arguments
validateFcn_reg_type = @(x) (x==1) || (x == 2);
validateFcn_balance = @(x) (x==0) || (x == 1);
validateFcn_init = @(x) (x==0) || (x == 1) || (x==2) || (x == 3) || (x == 4);
validateFcn_iter = @(x) isscalar(x) && (x>0) && (mod(x,1)==0);
validateFcn_consensus = @(x) isscalar(x) && (x>0) && (mod(x,1)==0);
validateFcn_kmin = @(x) isscalar(x) && (x>0) && (mod(x,1)==0);
validateFcn_kmax = @(x,y) isscalar(x) && (x>0) && (mod(x,1)==0) && (x>y);
validateFcn_kstep = @(x,y,z) isscalar(x) && (x>0) && (mod(x,1)==0) && (x+y<z);
validateFcn_cvfold = @(x) isscalar(x) && (x>0) && (mod(x,1)==0);
validateFcn_vo = @(x) (x==0) || (x == 1);

if(~validateFcn_reg_type(params.reg_type))
    error('run_hydra_experiment_csv:argChk','Input regularization type (reg_type) should be either 1 or 2!');
end
if(~validateFcn_balance(params.balanceclasses))
    error('run_hydra_experiment_csv:argChk','Input balance classes (balance) should be either 1 or 2!');
end
if(~validateFcn_init(params.init_type))
    error('run_hydra_experiment_csv:argChk','Initialization type can be either 0, 1, 2, 3, or 4!');
end
if(~validateFcn_iter(params.numiter))
    error('run_hydra_experiment_csv:argChk','Number of iterations should be a positive integer!');
end
if(~validateFcn_consensus(params.numconsensus))
    error('run_hydra_experiment_csv:argChk','Number of clustering consensus steps should be a positive integer!');
end
if(~validateFcn_kmin(params.kmin))
    error('run_hydra_experiment_csv:argChk','Minimum number of clustering solutions to consider should be a positive integer!');
end
if(~validateFcn_kmax(params.kmax,params.kmin))
    error('run_hydra_experiment_csv:argChk','Maximum number of clustering solutions to consider should be a positive integer that is greater than the minimum number of clustering solutions!');
end
if(~validateFcn_kstep(params.kstep,params.kmin,params.kmax))
    error('run_hydra_experiment_csv:argChk','Step number of clustering solutions to consider should be a positive integer that is between the minimun and maximum number of clustering solutions!');
end
if(~validateFcn_cvfold(params.cvfold))
    error('run_hydra_experiment_csv:argChk','Number of folds for cross-validation should be a positive integer!');
end
if(~validateFcn_vo(params.vo))
    error('run_hydra_experiment_csv:argChk','VO parameter should be either 0 or 1!');
end

disp('Done');
disp('HYDRA runs with the following parameteres');
disp(['featureCSV: ' featureCSV]);
disp(['OutputDir: ' outputDir]);
disp(['covCSV:' p.Results.covCSV])
disp(['C: ' num2str(params.C)]);
disp(['reg_type: ' num2str(params.reg_type)]);
disp(['balanceclasses: ' num2str(params.balanceclasses)]);
disp(['init_type: ' num2str(params.init_type)]);
disp(['numiter: ' num2str(params.numiter)]);
disp(['numconsensus: ' num2str(params.numconsensus)]);
disp(['kmin: ' num2str(params.kmin)]);
disp(['kmax: ' num2str(params.kmax)]);
disp(['kstep: ' num2str(params.kstep)]);
disp(['cvfold: ' num2str(params.cvfold)]);
disp(['vo: ' num2str(params.vo)]);

% csv with features
fname=p.Results.featureCSV;
if (~exist(fname,'file'))
    error('run_hydra_experiment_csv:argChk','Input feature .csv file does not exist');
end

% csv with features
covfname=p.Results.covCSV;
if(~isempty(covfname))
    if(~exist(covfname,'file'))
        error('run_hydra_experiment_csv:argChk','Input covariate .csv file does not exist');
    end
end 

% input data
% assumption is that the first column contains IDs, and the last contains
% labels
disp('Loading features...');
input=readtable(fname,'Delimiter','comma');
ID=input{:,1};
XK=input{:,2:end-1};
Y=input{:,end};

%% Begin
Y = Y(Y==-1); XK = XK(Y==-1,:); %taking only controls -1

% input covariate information if necesary
if(~isempty(covfname))
    disp('Loading covariates...');
    covardata = readtable(covfname,'Delimiter','comma') ;
    IDcovar = covardata{:,1};
    covar = covardata{:,2:end};
    disp('Done');
end
covar = covar(Y==-1,:);%taking only controls -1

YNC = Y(1:ceil(size(Y,1)/2));
YP = -Y(ceil(size(Y,1)/2)+1:end);%making pseudo-patient group
Y = cat(1,YNC,YP);

n = length(Y) ;
idx = randperm(n) ;
YY = Y ;
YY(idx,1) = Y(:,1);  % colum element randomly 
Y = YY;

%% End 

% z-score imaging features
XK=zscore(XK);
covar = zscore(covar);
disp('Done');

        
% NOTE: we assume that the imaging data and the covariate data are given in
% the same order. No test is performed to check that. By choosing to have a
% verbose output, you can have access to the ID values are read by the
% software for both the imaging data and the covariates

% verify that we have covariate data and imaging data for the same number
% of subjects
if(~isempty(covfname))
    if(size(covar,1)~=size(XK,1))
        error('run_hydra_experiment_csv:argChk','The feature .csv and covariate .csv file contain data for different number of subjects');
    end
end

% residualize covariates if necessary
if(~isempty(covfname))
    disp('Residualize data...');
    [XK0,~]=GLMcorrection(XK,Y,covar,XK,covar);
    disp('Done');
else
    XK0=XK; 
end

% for each realization of cross-validation
clustering=params.kmin:params.kstep:params.kmax;
part=make_xval_partition(size(XK0,1),params.cvfold); %Partition data to 10 groups for cross validation
% for each fold of the k-fold cross-validation
disp('Run HYDRA...');
for f=1:params.cvfold
    % for each clustering solution
    for kh=1:length(clustering)
        params.k=clustering(kh);
        disp(['Applying HYDRA for ' num2str(params.k) ' clusters. Fold: ' num2str(f) '/' num2str(params.cvfold)]);
        model=hydra2018(XK0(part~=f,:),Y(part~=f,:),[],params);
        YK{kh}(part~=f,f)=model.Yhat;
    end
end
disp('Done');

disp('Estimating clustering stabilitiy...')
% estimate cluster stability for the cross-validation experiment
ARI = zeros(length(clustering),1);
for kh=1:length(clustering)
    tmp=cv_cluster_stability(YK{kh}(Y~=-1,:));
    ARI(kh)=tmp(1);
    %% GC ARI(kh,:)=tmp(1,:);
end
disp('Done')

disp('Estimating final consensus group memberships...')
% Computing final consensus group memberships
CIDX=-ones(size(Y,1),length(clustering)); %variable that stores subjects in rows, and cluster memberships for the different clustering solutions in columns
for kh=1:length(clustering)
    CIDX(Y==1,kh)=consensus_clustering(YK{kh}(Y==1,:),clustering(kh));
end
disp('Done')

disp('Saving results...')
if(params.vo==0)
    save([outputDir '/HYDRA_results.mat'],'ARI','CIDX','clustering','ID','model');
else
    save([outputDir '/HYDRA_results.mat'],'ARI','CIDX','clustering','ID','XK','Y','covar','IDcovar','model');
end
disp('Done')
end

function [score,stdscore]=cv_cluster_stability(S)
k=0;
for i=1:size(S,2)-1
    for j=i+1:size(S,2)
        k=k+1;
        zero_idx=any([S(:,i) S(:,j)]==0,2);
    [a(k),b(k),c(k),d(k)]=RandIndex(S(~zero_idx,i),S(~zero_idx,j));
    end
end
% score=[mean(a) mean(b) mean(c) mean(d)];
% stdscore=[std(a) std(b) std(c) std(d)];

%% :GC 
score=[a b c d];
stdscore=[a b c d];
end

function [AR,RI,MI,HI]=RandIndex(c1,c2)
%RANDINDEX - calculates Rand Indices to compare two partitions
% ARI=RANDINDEX(c1,c2), where c1,c2 are vectors listing the 
% class membership, returns the "Hubert & Arabie adjusted Rand index".
% [AR,RI,MI,HI]=RANDINDEX(c1,c2) returns the adjusted Rand index, 
% the unadjusted Rand index, "Mirkin's" index and "Hubert's" index.
%
% See L. Hubert and P. Arabie (1985) "Comparing Partitions" Journal of 
% Classification 2:193-218

%(C) David Corney (2000)   		D.Corney@cs.ucl.ac.uk

if nargin < 2 | min(size(c1)) > 1 | min(size(c2)) > 1
   error('RandIndex: Requires two vector arguments')
   return
end

C=Contingency(c1,c2);	%form contingency matrix

n=sum(sum(C));
nis=sum(sum(C,2).^2);		%sum of squares of sums of rows
njs=sum(sum(C,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(C.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
D=  -t2+t3;		%no. disagreements

if t1==nc
   AR=0;			%avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end

RI=A/t1;			%Rand 1971		%Probability of agreement
MI=D/t1;			%Mirkin 1970	%p(disagreement)
HI=(A-D)/t1;	%Hubert 1977	%p(agree)-p(disagree)

function Cont=Contingency(Mem1,Mem2)

if nargin < 2 | min(size(Mem1)) > 1 | min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end
end
end

function IDXfinal=consensus_clustering(IDX,k)
[n,~]=size(IDX);
cooc=zeros(n);
for i=1:n-1
    for j=i+1:n
        cooc(i,j)=sum(IDX(i,:)==IDX(j,:));
    end
    %cooc(i,i)=sum(IDX(i,:)==IDX(i,:))/2;
end
cooc=cooc+cooc';
L=diag(sum(cooc,2))-cooc;

Ln=eye(n)-diag(sum(cooc,2).^(-1/2))*cooc*diag(sum(cooc,2).^(-1/2));
Ln(isnan(Ln))=0;
[V,~]=eig(Ln);
try
    IDXfinal=kmeans(V(:,1:k),k,'emptyaction','drop','replicates',20);
catch
    disp('Complex Eigenvectors Found...Using Non-Normalized Laplacian');
    [V,~]=eig(L);
    IDXfinal=kmeans(V(:,1:k),k,'emptyaction','drop','replicates',20);
end

end

function [part] = make_xval_partition(n, n_folds)
% MAKE_XVAL_PARTITION - Randomly generate cross validation partition.
%
% Usage:
%
%  PART = MAKE_XVAL_PARTITION(N, N_FOLDS)
%
% Randomly generates a partitioning for N datapoints into N_FOLDS equally
% sized folds (or as close to equal as possible). PART is a 1 X N vector,
% where PART(i) is a number in (1...N_FOLDS) indicating the fold assignment
% of the i'th data point.

% YOUR CODE GOES HERE

s=mod(n,n_folds);r=n-s;
p1=ceil((1:r)/ceil(r/n_folds));
p2=randperm(n_folds);p2=p2(1:s);
p3=[p1 p2];
part=p3(randperm(size(p3,2)));
end

function [X0train,X0test]=GLMcorrection(Xtrain,Ytrain,covartrain,Xtest,covartest)

X1=Xtrain(:,:);
C1=covartrain(:,:);
B=[C1 ones(size(C1,1),1)];
Z=X1'*B*inv(B'*B);
X0train=(Xtrain'-Z(:,1:end-1)*covartrain')';
X0test=(Xtest'-Z(:,1:end-1)*covartest')';
end

function printhelp()
fprintf('function returns estimated subgroups by hydra for clustering\n');
fprintf(' configurations ranging from K=1 to K=10, or another specified range of\n');
fprintf(' values. The function returns also the Adjusted Rand Index that was\n');
fprintf(' calculated across the cross-validation experiments and comparing\n');
fprintf(' respective clustering solutions.\n');
fprintf('\n');
fprintf(' INPUT\n');
fprintf('  \n');
fprintf(' REQUIRED \n');
fprintf(' featureCSV : .csv file containing the input features. (REQUIRED)\n');
fprintf('              every column of the file contains values for a feature, with \n');
fprintf('              the exception of the first and last columns. We assume that \n');
fprintf('              the first column contains subject identifying information \n');
fprintf('              while the last column contains label information. First line \n');
fprintf('              of the file should contain header information. Label\n');
fprintf('              convention: -1 -> control group - 1 -> pathological group\n');
fprintf('              that will be partioned to subgroups \n');
fprintf(' outputDir : directory where the output from all folds will be saved (REQUIRED)\n');
fprintf('\n');
fprintf(' OPTIONAL\n');
fprintf('\n');
fprintf(' covCSV : .csv file containing values for different covariates, which\n');
fprintf('           will be used to correct the data accordingly (OPTIONAL). Every\n');
fprintf('           column of the file contains values for a covariate, with the\n');
fprintf('           exception of the first column, which contains subject\n');
fprintf('           identifying information. Correction is performed by solving a \n');
fprintf('           solving a least square problem to estimate the respective \n');
fprintf('           coefficients and then removing their effect from the data. The \n');
fprintf('           effect of ALL provided covariates is removed. If no file is \n');
fprintf('           specified, no correction is performed.\n');
fprintf('\n');
fprintf(' NOTE: featureCSV and covCSV files are assumed to have the subjects given\n');
fprintf('       in the same order in their rows\n');
fprintf('\n');
fprintf(' C : regularization parameter (positive scalar). smaller values produce \n');
fprintf('     sparser models (OPTIONAL - Default 0.25)\n');
fprintf(' reg_type : determines regularization type. 1 -> promotes sparsity in the\n');
fprintf('            estimated hyperplanes - 2 -> L2 norm (OPTIONAL - Default 1)\n');
fprintf(' balance : takes into account differences in the number between the two\n');
fprintf('           classes. 1-> in case there is mismatch between the number of\n');
fprintf('           controls and patient - 0-> otherwise (OPTIONAL - Default 1)\n');
fprintf(' init : initialization strategy. 0 : assignment by random hyperplanes \n');
fprintf('        (not supported for regression), 1 : pure random assignment, 2: \n');
fprintf('        k-means assignment, 3: assignment by DPP random \n');
fprintf('        hyperplanes (default)\n');
fprintf(' iter : number of iterations between estimating hyperplanes, and cluster\n');
fprintf('        estimation. Default is 50. Increase if algorithms fails to\n');
fprintf('        converge\n');
fprintf(' numconsensus : number of clustering consensus steps. Default is 20.\n');
fprintf('                Increase if algorithm gives unstable clustering results.\n');
fprintf(' kmin : determines the range of clustering solutions to evaluate\n');
fprintf('             (i.e., kmin to kmax). Default  value is 1.\n');
fprintf(' kmax : determines the range of clustering solutions to evaluate\n');
fprintf('             (i.e., kmin to kmax). Default  value is 10.\n');
fprintf(' kstep: determines the range of clustering solutions to evaluate\n');
fprintf('             (i.e., kmin to kmax, with step kstep). Default  value is 1.\n');
fprintf(' cvfold: number of folds for cross validation. Default value is 10.\n');
fprintf(' vo : verbose output (i.e., also saves input data to verify that all were \n');
fprintf('      read correctly. Default value is 0\n');
fprintf('\n');
fprintf(' OUTPUT: \n');
fprintf(' CIDX: sub-clustering assignments of the disease population (positive\n');
fprintf('       class).\n');
fprintf(' ARI: adjusted rand index measuring the overlap/reproducibility of\n');
fprintf('      clustering solutions across folds\n');
fprintf('\n');
fprintf(' NOTE: to compile this function do\n');
fprintf(' mcc -m  run_hydra_experiment_csv \n');
fprintf(' \n');
fprintf(' \n');
fprintf(' EXAMPLE USE: \n');
fprintf(' run_hydra_experiment_csv(''test.csv'',''.'',''covCSV'',''test_covar.csv'',''C'',0.25,''kmax'',3,''init'',3,''cvfold'',5) \n');

end
