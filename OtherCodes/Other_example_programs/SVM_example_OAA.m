%# laod dataset
function main()
    S = load('fisheriris');
    data = zscore(S.meas);
    labels = grp2idx(S.species);

    %# cross-validate using one-vs-all approach
    opts = '-s 0 -t 2 -c 1 -g 0.25';    %# libsvm training options
    nfold = 10;
    acc = libsvmcrossval_ova(labels, data, opts, nfold);
    fprintf('Cross Validation Accuracy = %.4f%%\n', 100*mean(acc));

    %# compute final model over the entire dataset
    mdl = libsvmtrain_ova(labels, data, opts);
end

function mdl = libsvmtrain_ova(y, X, opts)
    if nargin < 3, opts = ''; end

    %# classes
    labels = unique(y);
    numLabels = numel(labels);

    %# train one-against-all models
    models = cell(numLabels,1);
    for k=1:numLabels
        
        models{k} = svmtrain(double(y==labels(k)), X, strcat(opts,' -b 1 -q'));
    end
    mdl = struct('models',{models}, 'labels',labels);
end

function [pred,acc,prob] = libsvmpredict_ova(y, X, mdl)
    %# classes
    labels = mdl.labels;
    numLabels = numel(labels);

    %# get probability estimates of test instances using each 1-vs-all model
    prob = zeros(size(X,1), numLabels);
    for k=1:numLabels
        [~,~,p] = svmpredict(double(y==labels(k)), X, mdl.models{k}, '-b 1 -q');
        prob(:,k) = p(:, mdl.models{k}.Label==1);
    end

    %# predict the class with the highest probability
    [~,pred] = max(prob, [], 2);
    %# compute classification accuracy
    acc = mean(pred == y);
end
function acc = libsvmcrossval_ova(y, X, opts, nfold, indices)
    if nargin < 3, opts = ''; end
    if nargin < 4, nfold = 10; end
    if nargin < 5, indices = crossvalidation(y, nfold); end

    %# N-fold cross-validation testing
    acc = zeros(nfold,1);
    for i=1:nfold
        testIdx = (indices == i); trainIdx = ~testIdx;
        mdl = libsvmtrain_ova(y(trainIdx), X(trainIdx,:), opts);
        [~,acc(i)] = libsvmpredict_ova(y(testIdx), X(testIdx,:), mdl);
    end
    acc = mean(acc);    %# average accuracy
end

function indices = crossvalidation(y, nfold)
    %# stratified n-fold cros-validation
    %#indices = crossvalind('Kfold', y, nfold);  %# Bioinformatics toolbox
    cv = cvpartition(y, 'kfold',nfold);          %# Statistics toolbox
    indices = zeros(size(y));
    for i=1:nfold
        indices(cv.test(i)) = i;
    end
end