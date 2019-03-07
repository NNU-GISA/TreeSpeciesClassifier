function SparseSVMMain()
    csvFileName = 'NF_IGFPCA6_EGF6_FULL_TREE.csv';
    csvSheetName = 'NF_IGFPCA6_EGF6_FULL_TREE';
    fetArrColumnDivision = [repmat([6 6],1,1)]; % the start column index of different fetures sets
    calculateOptSVM_C_Param = false;
    dataPerClass = 50;
    noOfClasses = 4;
    trainDataPercentage = 60;
    totalNumberOfTrials = 10;
    
    SparseSVM(csvFileName, csvSheetName, fetArrColumnDivision, calculateOptSVM_C_Param, dataPerClass ,noOfClasses, trainDataPercentage, totalNumberOfTrials);
end

function SparseSVM(csvFileName, csvSheetName, facDivision, calculateOptSVMParam, dataPerClass ,noOfClasses, trainDataPercentage, totalNumberOfTrials)
    %clf;
    % Set random states to avoid different results each time
    %randn('state',23432); rand('state',3454);
    
    sheet = csvSheetName; % current excel file data sheet name    
    Xorg = xlsread(csvFileName,sheet);
    %Xorg = [Xorg;Xorg;Xorg;Xorg];
    totFeatures = sum(facDivision);
    maxColumnIndice = totFeatures; labelColumnIndex = totFeatures+1; 
    arrDim = 1:totFeatures;
    
    %totNoOfSamples =  dataPerClass*noOfClasses;
    %trainingSetSize = totNoOfSamples*(trainDataPercentage/100);
    
    Finaccuracuy =[];
    for i=1:1:totalNumberOfTrials
                  
        [trainDataIndiceSet, testDataIndiceSet] = getStratifiedSamples(dataPerClass*(trainDataPercentage/100), dataPerClass, noOfClasses);
    
        trainDataIndice =  trainDataIndiceSet; %randperm(totNoOfSamples,trainingSetSize); % 80% data for training
        trainDataFinal = Xorg(trainDataIndice,arrDim);
        classFinal = Xorg(trainDataIndice,labelColumnIndex);

        testDataIndice = testDataIndiceSet;  %setdiff(1:1:size(Xorg,1),trainDataIndice);
        testDatafinal = Xorg(testDataIndice,arrDim);
        testClassfinal = Xorg(testDataIndice,labelColumnIndex);
        
        optCt = getOptimalC(Xorg, trainDataIndice, 4, labelColumnIndex, arrDim, maxColumnIndice, sheet, facDivision);        
        
        Finaccuracuy = [Finaccuracuy getAccu(trainDataFinal,classFinal,testDatafinal, testClassfinal,maxColumnIndice,arrDim, sheet, facDivision, true, optCt)]; 
    end
    FinalAccuracy = mean(Finaccuracuy);
    
end

function optimalC = getOptimalC(Xorg, trainDataIndiceOriginal, crossFold, labelColumnIndex, arrDim, maxColumnIndice, sheet, facDivision)
    trainDataIndiceOriginal =  trainDataIndiceOriginal(randperm(length(trainDataIndiceOriginal)));
    %dataPerClass = length(trainDataIndiceOriginal)/4; noOfClasses =4; trainDataPercentage =100;%
    %[trnDataIndxArr, testDataIndxArr] = getStratifiedSamples(dataPerClass*(trainDataPercentage/100), dataPerClass, noOfClasses);
    
    accuracuy =[];
    cntw = 1; mAccArr =[]; Arrct1=[];
    for cc = -5:5:35
        Ct = 2^cc;
        %crossFold = 4; %size of training data should be perfectly divisible by this crossFold value
        batchSize = size(trainDataIndiceOriginal,2)/crossFold;
        %batchTrainDataIndice = zeros(batchSize,1);

        trainDataIndiceSet = {};
        for trnDataGenloopCnt = 1:1:crossFold
            trainDataIndiceSet{trnDataGenloopCnt} = trainDataIndiceOriginal(batchSize*(trnDataGenloopCnt-1)+1:batchSize*(trnDataGenloopCnt)); 
        end

        for bthCnt = 1 : 1 : crossFold

            trainDataIndice =[]; testData = []; testClass =[];
            for k =1:1:crossFold
                if(k==bthCnt)
                    d =trainDataIndiceSet{k};
                    testData = Xorg(d,:);
                    testClass = Xorg(d,labelColumnIndex);
                   %continue;
                else
                   trainDataIndice = [trainDataIndice  trainDataIndiceSet{k};];   
                end
            end

            trainDataFinal = Xorg(trainDataIndice,:);
            classFinal = Xorg(trainDataIndice,labelColumnIndex);

            %ddd = sum(classFinal==3);
            
            testDataIndice =  setdiff(1:1:size(Xorg,1),trainDataIndice);
            testDatafinal = Xorg(testDataIndice,arrDim);
            testClassfinal = Xorg(testDataIndice,labelColumnIndex);

            %testDatafinal = trainData(:,1:maxColumnIndice)*1;
            %testClassfinal =  testData(:,1:maxColumnIndice)*1;

            accuracuy = [accuracuy getAccu(trainDataFinal,classFinal,testDatafinal, testClassfinal,maxColumnIndice,arrDim, sheet, facDivision, false, Ct)];
        end 

        Arrct1(cntw) = cc;
        mAccArr(cntw) = mean(accuracuy);
        cntw = cntw+1;
    end
    
    optimalC = max(cc);
    plot(Arrct1,mAccArr,'-*');
end

function accuracuy = getAccu(trainData,class,testData, testClass, maxColumnIndice,arrDim,sheet,facDivision, isplotOn,Ct)
    %countComb =1;
    %opcell = cell(4095,2);

    for nComb = 12:1:12
        Comb = combnk(arrDim,maxColumnIndice);

        for j = 1:1:size(Comb,1)

             % One against all approach implementation
             X=trainData(:,Comb(j,:));

             class1= class;
             class1(class1~=1) = -1;
             class1(class1==1) = 1;

             [wtrain1,btrain1] = SpaseSVMClassifer(X,class1,Ct,false);

             class2= class;
             class2(class2~=2) = -1;
             class2(class2==2) = 1;

             [wtrain2,btrain2]  = SpaseSVMClassifer(X,class2,Ct,false);

             class3= class;
             class3(class3~=3) = -1;
             class3(class3==3) = 1;

             [wtrain3,btrain3] = SpaseSVMClassifer(X,class3,Ct,false);

             class4= class;
             class4(class4~=4) = -1;
             class4(class4==4) = 1;
             [wtrain4,btrain4]  = SpaseSVMClassifer(X,class4,Ct,false);

             w = mean([abs(wtrain1) abs(wtrain2) abs(wtrain3) abs(wtrain4)]',1);         
             %if(strcmp(sheet,'Sheet2') && maxColumnIndice == 12)
              % w(9)=1.2; w(10)=2.25; %w(11)=0.015;
              %w(5) = 4;
            % w(9) = 9;  w(10) =10;
             %end
               % w = abs(wtrain1);  
            w= w/max(w); %w(7) = 0.1; w(8)= 0; w(9) = 0.45;  w(10) =0.5; w(11) =0.05; w(12) =0.1;
            
            if(isplotOn)           
                figure(2);
                for cnt = 0:1:size(facDivision,2)-1;
                    startIndx = sum(facDivision(1:cnt))+1;
                    endIndx = sum(facDivision(1:cnt+1));
                    stems(w(startIndx:endIndx),sum(facDivision(1:cnt))+1);
                    hold on; 
                end
                
                fetLbls = [];
                for i = 1:1:size(facDivision,2);
                    for jj=1:1:facDivision(i)
                        
                        if(i<10)
                            fetLbls = [fetLbls; strcat('F0',num2str(i),num2str(jj))];
                        else
                            fetLbls = [fetLbls; strcat('F',num2str(i),num2str(jj))];
                        end
                    end
                end                                
                set(gca,'XTick', 1:1:maxColumnIndice+2);
                set(gca,'fontsize',24)
                set(gca,'XTickLabel',{fetLbls});
                xlabel('Features', 'FontSize', 30);
                ylabel('w', 'FontSize', 30);
            end
            
             y1=@(x) ((x*wtrain1)-btrain1);
             y2=@(x) ((x*wtrain2)-btrain2);
             y3=@(x) ((x*wtrain3)-btrain3);
             y4=@(x) ((x*wtrain4)-btrain4);

             op =[y1(testData)';y2(testData)';y3(testData)';y4(testData)'];
             [~, indx] = max(op);

             refMat = testClass;

             % Accuracy assessment
             [C,~] = confusionmat(refMat,indx)
             UserAcc = max(C')./sum(C');
             ProdAcc= max(C)./sum(C);
             accuracuy= sum(diag(C))/sum(sum(C))
             %sd = strcat(char(Comb(j,:)+64));

             %opcell{countComb,1} =  sd;
             %opcell{countComb,2} = num2str(accuracuy);           

             %countComb=countComb+1;
        end
    end

end
function stems(x,offset)
        % STEM Display a Matlab sequence, x, using a stem plot.
        %// First define sequence from [0,N-1]
        vals = 0:numel(x)-1;
        %// Now use the above and manually shift the x coordinate
        stem(vals+offset,x, 'filled','MarkerSize',15);
end

function [wtrain,btrain] = SpaseSVMClassifer(X,class,Ct, isPlotOn)
    % This is the cost of slack for SVM.
    n=size(X,1);
    d=size(X,2);
    Yt=class;%Y(1:n,1); %pick the X and Y for the particular SVM inside the loop
    Xt=X;%X(1:n,:); %X is the same and Y is different
    %Ct= 2^-2;%norm(mean(abs(Xt))); %better estimate for C
    
    % Optimization problem solver
    cvx_begin %classical svm
        variables wtrain(d) e(n) btrain t
        %dual variable alphatrain
        minimize( 0.5*(t) + Ct*sum(e)) %norm(w) almost works except it takes an extra sqrt
        subject to
            diag(Yt)*(Xt*wtrain-btrain)-1 +e > 1;   %:alphatrain;
            norm(wtrain,1) < 0.01*t;
            norm(wtrain,2) < t;
            e>0; %slack
            
    cvx_end

   % y=@(x) sign((x*wtrain)-btrain);
    %classification = y(X)';
    retFetaureSparsity = [wtrain;btrain];%'./max(wtrain);
    
    if(isPlotOn)
        % Plot the data:
        figure(1); clf;
        hh = {};
        hh{1}=plot(X(class==1,1), X(class==1,2), 'o' );
        hold all
        hh{2}=plot(X(class==-1,1), X(class==-1,2), '*' );
        legend('Class 1','Class 2');
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        legend([hh{1:2}],'Class 1','Class2');

        % Plot
        m       = wtrain; % slope
        b       = btrain;   % intercept
        grid    = linspace(-2.5,2,100);
        hh{4}   = plot( grid, (b-m(1)*grid)/m(2) );
        xlim(xl); ylim(yl);
        legend([hh{:}],'Class 1','Class2','SVM','Sparse SVM');

    end
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
