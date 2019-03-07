function main()
    clear all;
    csvFileName = 'NF_IGFKM_EGF6_FULL_TREE_original.csv';
    csvSheetName = 'NF_IGFKM_EGF6_FULL_TREE_origina';
    facDivision = [repmat([6 6],1,1)]; % the start column index of different fetures sets
    dataPerClass = 50;
    noOfClasses = 4;
    trainDataPercentage = 60;
    numberOfIteration = 20;
    labelColumnIndex = 13;
    
    rbfKernel = @(X,Y,sigma) exp((-sigma)*(pdist2(X,Y,'euclidean').^2));
    %% loading data extacting train and test sampels
    Xorg = xlsread(csvFileName,csvSheetName);
    totFeatures = sum(facDivision);
    arrDim = 1:totFeatures;
    meas = Xorg(:,arrDim); species = Xorg(:,labelColumnIndex);

    data=[species meas];
    labels=species;%data(:,5); % the labels
    totFeatures = sum(facDivision);
    maxColumnIndice = totFeatures; 
    data(:,5)=[]; % the whole data

    AccArray =[];
    
    for cnt=1:1:numberOfIteration
    
        [trainDataIndiceSet, testDataIndiceSet] = getStratifiedSamples(dataPerClass*(trainDataPercentage/100), dataPerClass, noOfClasses);

        train_data = Xorg(trainDataIndiceSet,arrDim);
        train_data_label = Xorg(trainDataIndiceSet,labelColumnIndex);

        test_data = Xorg(testDataIndiceSet,arrDim);
        test_data_label = Xorg(testDataIndiceSet,labelColumnIndex);

        %% Data normalization (this part strongly affect the performance !! )  

        [~,scale_factor] = mapstd(train_data');
        test_data=mapstd('apply',test_data',scale_factor);
        train_data=mapstd('apply',train_data',scale_factor);
        train_data=train_data';
        test_data=test_data';
        
        %% Construncting 4 basis kernels with changing the paramters of RBF kernel 
        % in here the Dim. of train kernel is 15*15*4 and Dim. test kernel is
        % 135*15*4, where 4 is the number of kernels 

        kp=[0.002,0.004,0.006,0.008,0.01,0.0120,0.0140,0.0160,0.0180]; % kernel paramters 
        nok=numel(kp);  % number of kernels

          for jj=1:nok
              train_kernel(:,:,jj)=(rbfKernel(train_data,train_data,kp(jj)));     
              test_kernel(:,:,jj)=(rbfKernel(test_data,train_data,kp(jj)));
          end

        %% classification wth One-against-all strategy
        CC= 2^4; %1e1; % this is the trade-off paramter of MKL
        Coefficients=[];

        for C=1: max(test_data_label)

            OAA_train_lbl=-1*ones(size(train_data_label)); 
            OAA_train_lbl(train_data_label==C)=1;
            OAA_test_lbl=-1*ones(size(test_data_label));
            OAA_test_lbl(test_data_label==C)=1;    


            % prelocating the composite kernels in the memory, composite train kernel=
            % KTR and composite test kernel is KTS
                KTR=zeros(size(train_kernel(:,:,1)));
                KTS=zeros(size(test_kernel(:,:,1)));
            % estimating the kernel weights for each class, 
            tic
                alfa=(SimpleMKL (train_kernel,OAA_train_lbl,CC,0)); % estimating weights using SimpleMKL
                mkl_time(C)=toc; % the optimization time of MKL for class C 
                Coefficients=[Coefficients;alfa]; % the matrix of all kernel weights
             % constructing basis kernels 
                for i=1:nok
                    KTR=KTR+(alfa(i)*train_kernel(:,:,i));
                    KTS=KTS+(alfa(i)*test_kernel(:,:,i));
                end
                   KTR=[(1:size(KTR,1))', KTR]; %these two lines are for LibSVM toolbox
                   KTS=[(1:size(KTS,1))', KTS];
               % train the SVM with 5-fold CV to estimate the parametrs using
               % LibSVM toolbox
               cspace=0.01:100:2000;      
                bestcv = 0;
                        for c_values = 1:numel(cspace)
                            cmd = ['-s 0 -t 4  -v 5 -c ', num2str(cspace(c_values))];  
                                cv = svmtrain(OAA_train_lbl,KTR, cmd);
                                if (cv >= bestcv)
                                  bestcv = cv; 
                                  bestc = cspace(c_values);
                                end 
                        end
                         cmd = ['-s 0 -t 4 -c ', num2str(bestc)];
                         model=svmtrain(OAA_train_lbl,KTR,cmd);

                   % classification       
                    [~, ~,dec]=svmpredict(OAA_test_lbl,KTS,model);
                   if model.Label(1)~=1
                        dec_mat(:,C)=dec*(-1);
                     else
                        dec_mat(:,C)=dec;
                   end
        end
        [~, final_labels]=max(dec_mat,[],2);

        OA=(sum(final_labels==test_data_label))./numel(test_data_label);
        opt_time=sum(mkl_time);

        accuracyPerc = sum(final_labels == test_data_label) ./ numel(test_data_label)   %# accuracy
        C = confusionmat(test_data_label, final_labels)
        % C./repmat(sum(C,1),4,1)   C./repmat(sum(C,2),1,4)
        AccArray = [AccArray OA];
    end
    
    %disp(strcat('Overal Accuracy of  SimpleMKL = ',num2str(OA)));
    %disp(strcat('Optimiazation time of  SimpleMKL = ',num2str(opt_time)));
    disp(strcat('Overal Accuracy of  SimpleMKL = ',num2str(mean(AccArray))));
    %CoefficientsArr = [CoefficientsArr Coefficients];
end

function [trnDataIndxArr, testDataIndxArr] = getStratifiedSamples(trainDataperClass, dataPerClass, noOfClasses)
    trnDataIndxArr=[];
    testDataIndxArr=[];
    for i = 1:noOfClasses
        startIndx = (dataPerClass*(i-1))+1;
        endIndex = dataPerClass*i;
        trainDataForClass = randperm(dataPerClass,trainDataperClass) + (dataPerClass*(i-1)) ;
        testDataForClass = setdiff(startIndx:endIndex,trainDataForClass);
        trnDataIndxArr = [trnDataIndxArr trainDataForClass];
        testDataIndxArr = [testDataIndxArr testDataForClass];
    end      
end
