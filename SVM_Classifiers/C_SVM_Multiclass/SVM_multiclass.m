% CLEAR MATLAB ENVIRONMENT VARIABLES
clear all; close all;
%randn('state',23432); rand('state',3454);
% LOAD DATA FOR 200 TREES. 1:12 FEATURES FOR TREES, AND COLUMN 13 IS LABELS
csvFileName = 'NF_IGFKM_EGF6_FULL_TREE_original.csv'; csvSheetName = 'NF_IGFKM_EGF6_FULL_TREE_origina';
Xorg = xlsread(csvFileName,csvSheetName);
dataSamples = Xorg(:,1:6); dataLabels = Xorg(:,13);

% TRAINING/TESTING CONFIGURATION
dataPerClass = 50; trainDataPercentage = 60; noOfClasses = 4;

% NUMBER OF TRAILS
NumTrails = 20;
maxAccurcay = 0;
CmatMax =[];
accuracyArray = [];

for cnt= 1:1:NumTrails

    % STRATIFIED DATA SAMPLING
    [trainDataIndice, testDataIndice] = getStratifiedSamples(dataPerClass*(trainDataPercentage/100), dataPerClass, noOfClasses);
    % Get training set
    trainingData = dataSamples(trainDataIndice,:);
    trainingDataLabels = dataLabels(trainDataIndice);
    % Get test set
    testingData = dataSamples(testDataIndice,:);
    testingDataLabels = dataLabels(testDataIndice);

    %if(cnt == 1) % FOR THE FIRST ITERATION ONLY
        % CROSS-VALIDATION FOR OPTIMAL PARAMETER ESTIMATION 
        [cmd,~]=Cross_Validation(trainingData,trainingDataLabels,0); % 0 = no. of processor cores
    %end
    % TRAIN C-SVM MODEL
    model = svmtrain(trainingDataLabels,trainingData, cmd);

%    w = (model.sv_coef' * full(model.SVs));
  %  ww = abs(sum(w,1)/3);
  %  stem(ww);

   % w= w/max(w); %w(7) = 0.1; w(8)= 0; w(9) = 0.45;  w(10) =0.5; w(11) =0.05; w(12) =0.1;
%     facDivision = [repmat([6 6],1,1)];
%     for cnt = 0:1:size(facDivision,2)-1;
%         startIndx = sum(facDivision(1:cnt))+1;
%         endIndx = sum(facDivision(1:cnt+1));
%         stems(w(startIndx:endIndx),sum(facDivision(1:cnt))+1);
%         hold on; 
%     end
%     fetLbls = {'B_\alpha','B_l','B_k','B_w','B_s','B_n','T_v','T_d','T_c','T_l','T_\sigma','T_h'};
% 
% 
%     % fetLbls = [];
%     for i = 1:1:size(facDivision,2);
%         for jj=1:1:facDivision(i)
% 
%             if(i<10)
%               %  fetLbls = [fetLbls; strcat('F0',num2str(i),num2str(jj))];
%             else
%                % fetLbls = [fetLbls; strcat('F',num2str(i),num2str(jj))];
%             end
%         end
%     end                                
%     set(gca,'XTick', 1:1:maxColumnIndice);
%     set(gca,'fontsize',24)
%     set(gca,'XTickLabel',fetLbls);
%     set(gca,'YTick', [0 0.2 0.4 0.6 0.8 1]);
%     set(gca,'YTickLabel',{'0.0', '0.2', '0.4' '0.6', '0.8', '1.0'});
%     % xlabel('Features', 'FontSize', 30);
%     ylabel('w', 'FontSize', 30);
%     axis([0.5 12.5 0 1.0])
    
    
    % TESTING C-SVM MODEL
    [predict_label,~, ~]=svmpredict(double(testingDataLabels),testingData, model);

    % ACCURACY PERCENTAGE
    accuracyPerc = sum(predict_label == testingDataLabels) ./ numel(testingDataLabels);    %# accuracy
    C = confusionmat(testingDataLabels, predict_label); %# confusion matrix
    accuracyArray = [accuracyArray accuracyPerc];
    
    % Store the best case accuracy and C matrix
    if(accuracyPerc>maxAccurcay)
        maxAccurcay = accuracyPerc;
        CmatMax = C;
    end
    
end

disp(CmatMax);
disp(maxAccurcay);

disp('average value');
disp(mean(accuracyArray))




