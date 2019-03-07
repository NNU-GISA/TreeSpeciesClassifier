function [cmd,gamma]=Cross_Validation(t,X,parallel)
%Cross_Validation SVM cross-validation
%   [cmd,gamma]=Cross_Validation(t,X,parallel) perform a cross-validation 
%   based on LibSVM support vector. 
%       Input parameters:
%            - t        : train samples
%            - X        : train lables
% (optional) - parallel : # of worker to use (if 0 revert to the old sequetial version)
% 
%       output parameters:
%            - cmd      : command string to be used with svmpredict
%            - gamma    : value of gamma
% 
% from the readme of LibSVM train:
%   If the '-v' option is specified, cross validation is
%   conducted and the returned model is just a scalar: cross-validation
%   accuracy for classification and mean-squared error for regression.
%  
% Author: Massimo Santoni 10-2014

if nargin < 2
  parallel = 0;
end

%% PARAMETERS
% start and stop the parallel pool
 startPool = true;
 poolSize = parallel; 
 profile = 'local';
 oldPool = exist('parpool','file')==0; %true if parpool doesn't exist -> use matlabpool
 
% plot graph
 plotGraph = false;

% number of values  and range of c to be evaluated
 cNum = 20;
 cMin = 0.1;
 cMax = 3e3;
 
% number of values and range of g to be evaluated
 gNum = 20;
 gMin = 0.1;
 gMax = 0.3e3;

cSpace=linspace(cMin,cMax,cNum);
gSpace=linspace(gMin,gMax,gNum);
accuracy=zeros(cNum*gNum,1);

%% parralel version
if(parallel>0)
    
    if(startPool)  %------------------------
        if(oldPool)
            if( matlabpool('size') > 0)
               matlabpool close;
            end
            matlabpool(profile, poolSize);
        else
            if( ~isempty(gcp('nocreate')))
                delete(gcp);
            end
            parpool(profile,poolSize);
        end
    end           %------------------------
    
    [p,q] = meshgrid(cSpace, gSpace);
    param = [p(:) q(:)];
    disp(['Parallel crossvalidation (' num2str(size(param,1)) ' iterations)']);
    
    parfor i=1:size(param,1)
        cmd = ['-v 5 -t 2 -c ', num2str(param(i,1)), ' -g ', num2str(param(i,2)) '-q'];
        accuracy(i) = svmtrain(double(X),double(t),cmd);    
    end
    
    if(startPool) %-----------------------
        if(oldPool)
            matlabpool close;
        else       
            delete(gcp);
        end
    end         %-------------------------
    
    [bestaccuracy,bestindex]=max(accuracy);
    bestc=param(bestindex,1);
    bestg=param(bestindex,2);
    disp(['Best parameters: -c ' num2str(bestc) ' -g ' num2str(bestg) ' | accuracy = '  num2str(bestaccuracy)]);

else
%% old sequential version
    bestcv = 1;
    for cIndex =1:cNum;
        c=cSpace(cIndex);
        for gIndex =1:gNum;
            g=gSpace(gIndex);
            cmd = ['-v 5 -t 2 -c ', num2str(c), ' -g ', num2str(g)];
            cv = svmtrain(X,t,cmd);
            if (cv >= bestcv)
                bestcv = cv; 
                bestc =c; 
                bestg=g;
                disp(['-  +  - -c ', num2str(c), ' -g ', num2str(g)]);
            end
        end
    end
end

cmd = ['-s 0 -t 2 -c ', num2str(bestc), ' -g ', num2str(bestg)];
gamma = bestg;

%% PLOT RESULT
if(plotGraph)
    x = reshape(param(:,1),gNum, cNum)';
    y = reshape(param(:,2),gNum, cNum)';
    z = reshape(accuracy,gNum, cNum)';
    surf(x,y,z)
    ylabel('g')
    xlabel('c')
end


