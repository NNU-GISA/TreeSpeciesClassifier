function beta=SimpleMKL (train_kernel,labels,C,verbose)
% addpath('D:\Saeid Folder\toolbox\SVM_KM')
% addpath('D:\Saeid Folder\toolbox\simplemkl')


% verbose=1;
options.algo='svmclass'; % Choice of algorithm in mklsvm can be either
                         % 'svmclass' or 'svmreg'

                         options.stopvariation=0; % use variation of weights for stopping criterion 
options.stopKKT=0;       % set to 1 if you use KKTcondition for stopping criterion    
options.stopdualitygap=1; % set to 1 for using duality gap for stopping criterion
options.seuildiffsigma=1e-2;        % stopping criterion for weight variation 
options.seuildiffconstraint=0.1;    % stopping criterion for KKT
options.seuildualitygap=0.01;       % stopping criterion for duality gap
options.goldensearch_deltmax=1e-1; % initial precision of golden section search
options.numericalprecision=1e-8;   % numerical precision weights below this value
                                   % are set to zero 
options.lambdareg = 1e-8;          % ridge added to kernel matrix 

options.firstbasevariable='first'; % tie breaking method for choosing the base 
                                   % variable in the reduced gradient method 
options.nbitermax=500;             % maximal number of iteration  
options.seuil=0;                   % forcing to zero weights lower than this 
options.seuilitermax=10;           % value, for iterations lower than this one 

options.miniter=0;                 % minimal number of iterations 
options.verbosesvm=0;              % verbosity of inner svm algorithm 
options.efficientkernel=0;         % use efficient storage of kernels 



% [beta,w,w0,pos,nbsv,SigmaH,obj] = mklmulticlass(K,tr_labels,C,max(labels),options,verbose);
[beta,Alpsup,w0,pos,history,obj,status] = mklsvm(train_kernel,labels,C,options,verbose);


end