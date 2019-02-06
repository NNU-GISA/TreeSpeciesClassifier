function [igffeatureArr, slicedfeatureArrTotal, numBranches] = getConiferIGFs_Method1(inFileFullName, speciesCode, allFiles, plotOn)

    igffeatureArr = zeros(size(allFiles,1),9);
    slicedfeature4Div = zeros(size(allFiles,1),32);
    slicedfeature7Div = zeros(size(allFiles,1),56);
    slicedfeature10Div = zeros(size(allFiles,1),80);
    slicedfeature13Div = zeros(size(allFiles,1),104);
    numBranches = zeros(size(allFiles,1),1);
    count = 1;
    for file = allFiles'
              
        fullFileName = strcat(inFileFullName,file.name); % get full file name

        %-----------------------------------------------------------------%
        % -- The below part of the code does LiDAR data Preprocessing --- %
        % ----------------------------------------------------------------%
        
        %Read LiDAR data
        data = LoadData(fullFileName);
        
        % Data Pre-processing and basic tree information retrieval 
        % e.g., max/min height, crown height etc.
        stData = dataPreProcessing(data);
        
        % Plot the conifer tree stem & also plots the complete LiDAR data
        showRawLiDARdata = true;  % 1 = show LiDAR data; 0 = hide LiDAR data
        plotLiDARdataWithStem(stData, showRawLiDARdata, plotOn);       
        
        % ----------------------------------------------------------------%
        % -  The below part of the code locates and plots the conifer 
        % -                        branch tips                          - %
        % ----------------------------------------------------------------%        
        plotBranchTips = false; plotConvexHull = false;
        extremeLiDARDataArray = getBoundaryPoints(stData, speciesCode, plotBranchTips, plotConvexHull, plotOn);
        
        % ----------------------------------------------------------------%
        % -- The below part of the code does inital region growing starting 
        % ----------------------- from branche tips ----------------------%
        % ----------------------------------------------------------------%
        
        % Get connected components in the LiDAR point cloud by using 3D        
        [connectComponentArray,connectComponentindexArr, extremeLiDARDataArray] = ...
            getConnectedComponetsOld(stData.lidarDataArray, extremeLiDARDataArray); 
        
        % ----------------------------------------------------------------%
        % ---- The below part of the completes the identification of -----%
        % ---- branch points not identified from region growing step -----%
        % ----------------------------------------------------------------%  
        
        % Get stem point in the direction of PC1 of branch clusters.
        [stemPointArr, D1] = getStemPointFromPCA(extremeLiDARDataArray,...
            connectComponentArray, stData.retMaxXYZ);
        
        % ----------------------------------------------------------------%
        % -------  The below part of the code calculates the IGFs ------- %
        % ----------------------------------------------------------------%
        
        % Get features from individual branches and calculates the Internal Geometric Features
        [retIGFsVec, BranchArr, nb] = getIGFFeatureValue(stData, extremeLiDARDataArray,...
            stemPointArr ,D1, connectComponentArray, connectComponentindexArr, speciesCode, plotOn);  
        numBranches(count) = nb;
        
        % Feature slicing (for slice-wise analysis)
        SlicedFeatures = FeatureSlicing(BranchArr,stData);
        slicedfeature4Div(count,:) = SlicedFeatures{1}; % 8 * 4 = 32 features
        slicedfeature7Div(count,:) = SlicedFeatures{2}; % 8 * 7 = 56 features
        slicedfeature10Div(count,:) = SlicedFeatures{3}; % 8 * 10 = 80 features
        slicedfeature13Div(count,:) = SlicedFeatures{4}; % 8 * 13 = 104 features       
        
        % Store IGFs for the tree 
        igffeatureArr(count,:) = retIGFsVec;
    
        % Print featurevalue in console and save results
        printInConsole = true;
        printResultInConsole(count, retIGFsVec, printInConsole);

        % Print plot if plotOn is true
        if(plotOn)
            hold off;
            pause(0.25); % To update plot
        end
        
        % Store the IGFs and IFs by slicing the by 4, 7, 10, and 13.
        slicedfeatureArrTotal{1} = slicedfeature4Div;
        slicedfeatureArrTotal{2} = slicedfeature7Div;
        slicedfeatureArrTotal{3} = slicedfeature10Div;
        slicedfeatureArrTotal{4} = slicedfeature13Div;
        
        count = count + 1;
    end
    
    % Save IGF values in the local drive as a .csv file.
   % csvwrite('csvlist_igf_pca.csv', igffeatureArr);
end

function printResultInConsole(count, retIGFarr, printInConsole)
    % Print results in the matlab console    
    if(printInConsole)
        if(count == 1)
            fprintf('  avgRetSlope   avgRetLen  avgPntDistToline  avgMaxEval  avgEigRatio  avgTotalPoints   avgInt  avgIntSd   avgMedian\n')
        end
        disp(strcat(num2str(retIGFarr(1)), {'    '}, num2str(retIGFarr(2)), {'        '} , num2str(retIGFarr(3)), {'        '}, num2str(retIGFarr(4)), {'    '},...
            num2str(retIGFarr(5)), {'         '}, num2str(retIGFarr(6)), {'    '}, num2str(retIGFarr(7)), {'   '}, num2str(retIGFarr(8)), {'    '},...
            num2str(retIGFarr(9))));
    end
end

function SlicedFeatures = FeatureSlicing(BranchArr,stData)
    SlicedFeatures = [];
    SliceNrArr = [4 7 10 13];
    CrownBase = stData.minCrownHeight;
    CrownTop = stData.maxtreeHeight;
    CrownHeight = CrownTop - CrownBase;
    TreeHeight = stData.treeHeight;
    TotalPointNr = size(stData.lidarDataArray,1);
    for SliceNrCounter = 1:size(SliceNrArr,2)
        SliceNr = SliceNrArr(SliceNrCounter);
        thresholds = [];
        for i = 1:(SliceNr)
            % calculate the height thresholds that will be used to
            % determine which slice a branch belongs to
            thresholds = [thresholds, i/SliceNr];
        end
        thresholds = CrownBase + thresholds*CrownHeight;
        treshTable = []; labels = [];
        % generate a label column which contains the slice the branch
        % belongs to. It will be added to the branch feature matrix
        for BranchCtr = 1:size(BranchArr,1)
            treshTable = [treshTable;BranchArr(BranchCtr,9) < thresholds];
            labels = [labels;min(find(treshTable(BranchCtr,:) == 1))];
        end
        % save the labeled branches for each different slice number
        % configuration in a cell
        LabeledBranches{SliceNrCounter} = [BranchArr,labels];
    end
    
    %calculate the average of each feature among branches in the same slice
    for OuterLoop = 1:4
        CurrentTable = LabeledBranches{OuterLoop};
        AvgIgfs = zeros(SliceNrArr(OuterLoop),6);
        AvgIntFeatures = zeros(SliceNrArr(OuterLoop),2);
        BranchCount = [];
        for SliceCtr = 1:SliceNrArr(OuterLoop)
            BranchIndexPerSlice = find(CurrentTable(:,10) == SliceCtr);
            BranchesPerSlice = CurrentTable(BranchIndexPerSlice,:);
            
            % counts the number of branches per slice
            % currently unused, but may be useful as an additional feature
            BranchCount = [BranchCount,size(BranchesPerSlice,1)];
            
            % table of IGFs and Intensity features of each slice
            IGFs = BranchesPerSlice(:,1:6);
            IntFeatures = BranchesPerSlice(:,7:8);
            
            % averaging features for the branches in each slice
            avgRetSlope = mean(IGFs(:,1));
            avgRetLen = mean(IGFs(:,2))/TreeHeight;
            avgPtToLineDst = mean(IGFs(:,3));
            avgMaxEval = mean(IGFs(:,4))/TreeHeight;
            avgEigRatio = mean(IGFs(:,5));
            if(avgEigRatio>100) % to avoid large values due to division by small number.
                avgEigRatio = 100;
            end
            avgTotalPoints = mean(IGFs(:,6))/TotalPointNr;
            AvgIgfs(SliceCtr,:) = [avgRetSlope,avgRetLen,avgPtToLineDst,avgMaxEval,avgEigRatio,avgTotalPoints];
            
            avgInt = mean(IntFeatures(:,1));
            avgIntSd = mean(IntFeatures(:,2));
            AvgIntFeatures(SliceCtr,:) = [avgInt,avgIntSd];
        end
        % save the features in a structure
        % STRUCTURE :
        % | IGFs Slice 1 | Intensity Features Slice 1| ... | IGFs Slice n | Intensity Features Slice n |
        StructuredFeatures = [];
        for j = 1:SliceNrArr(OuterLoop)
        StructuredFeatures = [StructuredFeatures,AvgIgfs(j,:),AvgIntFeatures(j,:)];
        StructuredFeatures(isnan(StructuredFeatures)) = 0; % substitute eventual NaN values with zero
        end
        SlicedFeatures{OuterLoop} = StructuredFeatures;   
    end
    
end

function extremeLiDARDataArray = getBoundaryPoints(stData, speciesCode, plotBranchTips, plotConvexHull, plotOn)
    K = boundary(stData.lidarDataArray(:,1),stData.lidarDataArray(:,2),stData.lidarDataArray(:,3),1); 
    K = unique(K);    
    tempLiDARarray = stData.lidarDataArray(K,:);     

    % Cuttoff distances for mutiple branch tip removal
    if(strcmp(speciesCode,'ar'))
        cutOffDistance = 2; 
    elseif(strcmp(speciesCode,'la'))
        cutOffDistance = 4;
    elseif(strcmp(speciesCode,'pc'))
        cutOffDistance = 3;
    else
        cutOffDistance = 2;
    end
    
    tempLiDARarray = removeRepetitions(stData, tempLiDARarray, cutOffDistance); %2:ar  0.2:la   %0.5:pc    0.5:ab
    
    % for further smoothening (remove any mutiple branch-tip candidates from the same branch)
    K1 = boundary(tempLiDARarray(:,1),tempLiDARarray(:,2),tempLiDARarray(:,3),0.5); 
    uniqueK1 = unique(K1);
    extremeLiDARDataArray = tempLiDARarray(uniqueK1,:);
    
    % Plot branch points if both plotBranchTips = true and plotOn = true
    if(and(plotBranchTips,plotOn))
        plot3(tempLiDARarray(uniqueK1,1), tempLiDARarray(uniqueK1,2), ...
        tempLiDARarray(uniqueK1,3), '.', 'Color', [1 0 0],'MarkerSize',20);
        hold on;
    end    
    
    % Plot convex hull if both plotConvexHull = true and plotOn = true
    if(and(plotConvexHull,plotOn))
        trisurf(K1,tempLiDARarray(:,1),tempLiDARarray(:,2),tempLiDARarray(:,3));
        alpha(0.2);
    end
end

function plotLiDARdataWithStem(stData, showRawLiDARParameter, plotOn)
    if(plotOn)
        clf; % clear figures            
        % plot the tree stem
        plot3([stData.retMaxXYZ(1) stData.retMaxXYZ(1)],[stData.retMaxXYZ(2)...
            stData.retMaxXYZ(2)],[0 stData.maxtreeHeight], '-o', 'Color', [1 0 0]);
        hold on;
        %To show cubic/sector grid or not ; both should not be true simulatniously 
        plotLiDARData(stData.lidarDataArrayComplete, false, false,...
        stData.htDeduction,showRawLiDARParameter, 15,stData.retMaxXYZ)

        camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; view(-45, 15); grid on;
        % maximum height
        mht = stData.maxtreeHeight + (2 - mod(stData.maxtreeHeight,2)); % round-pff to next mutiple of 2;
        % set maximum axis dimentions to be shown
        axis([-stData.maxTreeWidth stData.maxTreeWidth -stData.maxTreeWidth ...
           stData.maxTreeWidth 0 mht]);

        % Set perspective and label fonts and view angle for the 3D plot
        xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 14); 
        ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 14); 
        zlabel('Tree height','Fontname', 'Times New Roman' ,'FontSize', 14);
        set(gca,'XLim',[-6 6]);set(gca,'YLim',[-6 6]);set(gca,'ZLim',[0 mht]);
        set(gca,'XTick',-6:2:6); set(gca,'YTick',-6:2:6); set(gca,'ZTick',0:2:mht)
        % set(gca,'XTickLabel',['0';'';'1';' ';'2';' ';'3';' ';'4'])
    end
end

function [unallotedlidarDataArr, unallotedlidarIndexArr] = getUnallocatedPointIndices(stData, connectComponentArray, connectComponentArrayIndx, stemPointArr, branchTipArray)
    tmpAra =[];
    for hh =1:1:size(connectComponentArray,2)        
        tmpAra = [tmpAra; connectComponentArrayIndx{hh}];
    end
    unallotedlidarDataArr = stData.lidarDataArrayComplete(find(~ismember(stData.lidarDataArrayComplete(:,7),unique(tmpAra))), :);  
    % allocate points to cluster of the nearest branch-line
    unallotedlidarIndexArr = zeros(size(branchTipArray,1),3);
    for ii = 1: 1: size(unallotedlidarDataArr,1)        
        distArrddss = getDist2Point22(branchTipArray(:,(1:3)), stemPointArr, repmat(unallotedlidarDataArr(ii,(1:3)),[size(branchTipArray,1) 1]) );        
        [mval, mdx] = min(distArrddss);  
        %if(mval <0.5)
        unallotedlidarIndexArr(ii,:) = [ii mdx mval];        % Point       LineNumber      Distance of Point to the line
        %end
    end    
    unallotedlidarIndexArr = sortrows(unallotedlidarIndexArr,2);

end


function [retIGFsVec, BranchArr, numBranches] = getIGFFeatureValue(stData, branchTipArray, stemPointArr, D1, connectComponentArray, connectComponentArrayIndx, speciesCode, plotOn)
    
    % initilaizing variables 
    numBranches = size(branchTipArray,1);
    retIGFsArray = zeros(numBranches,9);
    
    
    % Random color generation for clusters
    retRand = distinguishable_colors(size(branchTipArray,1)); %rand(size(branchTipArray,1),3);

    % Indices of previosly unallocated liDAR point
    [unallotedlidarDataArr, unallotIndexArr] = ....
        getUnallocatedPointIndices(stData, connectComponentArray, connectComponentArrayIndx, stemPointArr, branchTipArray);
    
    % Large slopes which are very likely to be wrong.(This is a known limitation of the method)        
    if(strcmp(speciesCode,'ar'))
        slopeThreshold = 30; 
        includePointDistance = 0.5; 
    elseif(strcmp(speciesCode,'la'))
        slopeThreshold = 30; 
        includePointDistance = 0.8; 
    elseif(strcmp(speciesCode,'pc'))
        slopeThreshold = 100; 
        includePointDistance = 0.2;
    else
        slopeThreshold = 30; 
        includePointDistance = 0.8; 
    end
    
    for bCount = 1:1:numBranches
        
        % Add the LiDAR point obtained from region growing to the tmpGrp points 
        RegionGrownLiDARPoints = stData.lidarDataArray(ismember(stData.lidarDataArray(:,7), connectComponentArrayIndx{bCount}),1:3);
        
        % distance to line is less than "includePointDistance"
        unallotedLiDARPoints = unallotedlidarDataArr(unallotIndexArr(and(unallotIndexArr(:,2)==bCount, ...
            (abs(unallotIndexArr(:,3))) < includePointDistance),1),(1:3));
        
        branchPointCluster = [RegionGrownLiDARPoints; unallotedLiDARPoints; ]; % repmat(stemPointArr(jj,:),10,1)
        prosConnectedComp{bCount} = branchPointCluster;

        % Regression fit line through connected components.
        [retSlope, lineLength, pntDistToline, maxeval, eigRatio, totalPoints] = getBestFitLine(stData,branchPointCluster, stemPointArr(bCount,:), branchTipArray(bCount,(1:3)), slopeThreshold, speciesCode, true, plotOn, [0 0.5 0]); % isPlotON

        % Plot branches with slope < slopeThreshold
        if(abs(retSlope) < slopeThreshold) %  && stemPointArr(bCount,3) > 5
            if(plotOn)
                plot3(branchPointCluster(:,1),branchPointCluster(:,2),branchPointCluster(:,3),'.','MarkerSize',5,'Color',[retRand(bCount,1) retRand(bCount,2) retRand(bCount,3)]);
            end
        end
        
         % get the intensity fratures (mean and standard deviation)
        [avgInt,sdInt,medianZ] = GetIntensityFeatures(branchPointCluster,stData);        
        retIGFsArray(bCount,:) = [retSlope, lineLength, pntDistToline, maxeval,eigRatio,totalPoints,avgInt,sdInt,medianZ]; 
       
    end
    
    %Get the number of branches
    numBranches = size(retIGFsArray,1);
    
    % Remove all the rows which has 0's in at least one column
    retIGFsArray = retIGFsArray(retIGFsArray(:,5)~=-1,:);
    retIGFsArray = retIGFsArray(retIGFsArray(:,1)~=0,:);
    retIGFsArray = retIGFsArray(retIGFsArray(:,2)~=0,:);
    retIGFsArray = retIGFsArray(retIGFsArray(:,3)~=0,:);
    retIGFsArray = retIGFsArray(retIGFsArray(:,4)~=0,:);
    retIGFsArray = retIGFsArray(retIGFsArray(:,5)~=0,:);
    retIGFsArray = retIGFsArray(retIGFsArray(:,6)~=0,:);
    retIGFsArray = retIGFsArray(retIGFsArray(:,7)~=0,:);
    retIGFsArray = retIGFsArray(retIGFsArray(:,8)~=0,:);
    retIGFsArray = retIGFsArray(retIGFsArray(:,9)~=0,:);
    
    % Select branches which meets the above slopeThreshold criteria
    BranchArr  = retIGFsArray(abs(retIGFsArray(:,1)) < slopeThreshold,:);    
    retIGFsVec = sum(BranchArr,1)/size(BranchArr,1);

    % Perform height/desnity normalization for three of the features
    retIGFsVec(2) = retIGFsVec(2)/stData.treeHeight; % Branch Length
    retIGFsVec(4) = retIGFsVec(4)/stData.treeHeight; % Max  eigenvalue
    retIGFsVec(6) = retIGFsVec(6)/size(stData.lidarDataArray,1); % Total Points in a branch
    if(retIGFsVec(5)>100) % to avoid large values due to division by small number. 
        retIGFsVec(5) = 100; %(eigratio)
    end
    
end

function retEigenData = getBestEigenDirection(inputData)
    [a,~,c]=pca(inputData);
    pc = zeros(3,2);
    
    for i=1:size(a,2)
        pc(:,i) = a(:,i);%'*-c(i);
    end
    retEigPoints{1} = pc;
    retEigPoints{2} = c;    
    retEigenData = retEigPoints;    
end

function [avgInt,sdInt,medianZ] = GetIntensityFeatures(tmpGrp,stData)
tmpGrpInt = [];    
branchPoints = [];
    % height = (stData.maxtreeHeight-stData.minCrownHeight); % could also be implemented with mintreeheight instead of mincrownheight
    
    % use tmpGrp xyz information to find corresponding points in lasdata
    % and extrapolate intensity from there (since tmpGrp does not include the intensity value)
    for i=1:size(tmpGrp,1)
        target=tmpGrp(i,:);
        indicesX = find ( stData.lidarDataArray(:,1) == target(1) ); % check x
        indicesY = find ( stData.lidarDataArray(:,2) == target(2) ); % check y
        indicesZ = find ( stData.lidarDataArray(:,3) == target(3) ); % check z
        
        indices=[indicesX;indicesY;indicesZ];
        uniqueIndices=unique(indices);
        count=histc(indices,uniqueIndices); %counts the occurrence of each index. the one that appears 3 times (xyz) is the right one
        targetIndex=uniqueIndices(find(count==3),1);
        branchPoints = [branchPoints;stData.lidarDataArray(targetIndex,3),stData.lidarDataArray(targetIndex,6)]; % save z and intensity information
    end
    % calculate the mean-range of the points' z coordinate (middle point of the tree height)
    BranchMax = max(branchPoints(:,1));
    BranchMin = min(branchPoints(:,1));
    medianZ = (BranchMax+BranchMin)/2; % will be used to identify to which height zone the branch belongs to
    
    avgInt = mean(branchPoints(:,2)); % calculates the average intensity value for the branch
    sdInt = std(branchPoints(:,2)); % calculates the intensity standard deviation value for the branch
end

function getDist = getDist2Point22(point1, point2, extPoint)
  %  s= [point2(1) - point1(1);point2(2) - point1(2);point2(3) - point1(3)];
  %  M0= extPoint;
  %  M1 = point1;
  %  m2m1 = [M1(1)-M0(1) M1(2)-M0(2) M1(3)-M0(3)];   
  %  getDistu = norm(cross(m2m1,s),2)/norm(s,2);
    % distArr = getDist2Point(branchTipArray((1:100),(1:3)), stemPointArr(1:100,:), lidarDataArraytmp(1:100,(1:3)));
    getDist = zeros(size(point1,1),1);
  
    x = extPoint; %some point
    a = point1; %segment points a,b
    b = point2;
    
    v = b-a;
    w = x-a;
    
    c1 = diag(w*v');
    c2 = diag(v*v');
    
    cless = find(c1<0);
    cmore= find(c2<c1);    
    remCind = setdiff(1:1:size(point1,1), unique([cmore; cless]));
    
    if(~isempty(cless))
       getDist(cless,1) =  sqrt(sum((x(cless,:)' - a(cless,:)').^ 2))'; 
       %return;
    end
    if(~isempty(cmore))
       getDist(cmore,1) =  sqrt(sum((x(cmore,:)' - b(cmore,:)').^ 2))'; 
       %return;
    end
    if(~isempty(remCind))
        b1= c1(remCind)./c2(remCind);
        pb= a(remCind,:) + repmat(b1,[1 3]).*v(remCind,:);
        getDist(remCind,1) =  sqrt(sum((x(remCind,:)' - pb').^ 2));   
    end
    
  %  d_ab = norm(a-b);
  %  d_ax = norm(a-x);
  %  d_bx = norm(b-x);

  %  if dot(a-b,x-b)*dot(b-a,x-a)>=0
  %      A = [a,1;b,1;x,1];
  %      getDist = abs(det(A))/d_ab;        
  %  else
  %      getDist = min(d_ax, d_bx);
  %  end    
end

function isPointWithin = isWinthinArc(v1, v2, inputPoints, radius)
    isPointWithin = [];
    for i =1:1:size(inputPoints,1)
        x = inputPoints(i,1);
        y = inputPoints(i,2);
        
        isClkwisetoV1 = -v1(1)*y + v1(2)*x > 0;    
        isClkwisetoV2 = -v2(1)*y + v2(2)*x > 0;  
        
        withinRad = x*x + y*y <= radius.^2;

        if(not(isClkwisetoV1) && isClkwisetoV2 && withinRad)
            isPointWithin = [isPointWithin i];    
        end
    end
end

function [connectComponentArray,connectComponentindexArr, seedPointsArray] = getConnectedComponetsOld(lidarDataArray, seedPointsArray)

    outerLoopCnt= 1;
    noOfIterations = 10; NumNearestPoints = 5;
    %connectComponentArray = zeros(NumNearestPoints*noOfIterations,3,(compEndIndex-compStartIndex));
    I =[];
    for clstStartPointIndex = 1:1:size(seedPointsArray,1)  %length(extremeLiDARDataArray)-50:1:length(extremeLiDARDataArray)-51+noPoints
            seletedIndices = [];
        cntd =1; clusteredLiDARDataArray = zeros(NumNearestPoints*noOfIterations,3); 
        startPoint = [seedPointsArray(clstStartPointIndex,(1:3)) seedPointsArray(clstStartPointIndex,(8))];
        startPoint1 = startPoint;
         % [double(startPoint(1)) double(startPoint(2)) double(startPoint(3))]
        otherPoints = [lidarDataArray(:,(1:3)) lidarDataArray(:,(8)) ];
        otherPoints = removerows(otherPoints,'ind',I);
        for loopCnt = 1:1:noOfIterations
            %(find(and(lidarDataArray(:,3) > startPoint(3)-10, lidarDataArray(:,3) < startPoint(3)+10)),(1:3));
            [~,I] = pdist2(otherPoints,startPoint,'euclidean','Smallest',NumNearestPoints);
            temotherPoints = [otherPoints(I,:)];
            sortLiDArray = sortrows(temotherPoints,[3]); %order by 3rd column in decresing order.            
            
            a1a = otherPoints(I(2),:); b2b = startPoint1(1,(1:4));
            dst = pdist([a1a;b2b],'euclidean');
            if(dst<2)
                clusteredLiDARDataArray(cntd:cntd+(NumNearestPoints-1),:) = otherPoints(I,1:3);
               % seletedIndices = [seletedIndices; I];
                cntd = cntd + NumNearestPoints;
            end
            startPoint = sortLiDArray(NumNearestPoints,(1:4));
            otherPoints(I,:) = 0;
        end
        
        clusteredLiDARDataArray = clusteredLiDARDataArray(clusteredLiDARDataArray(:,3)~=0,:);
        connectComponentArray{outerLoopCnt} = clusteredLiDARDataArray;%zeros(NumNearestPoints*noOfIterations,3,(compEndIndex-compStartIndex));
        
        for ccCont =1:1:size(clusteredLiDARDataArray,1)
            
            aa = clusteredLiDARDataArray(ccCont,:);
            
            %idxfd =  find(and(and(aa(:,1)==lidarDataArray(:,1) ,aa(:,2)==lidarDataArray(:,2)), aa(:,3)==lidarDataArray(:,3) ));
            idx =  lidarDataArray(find(and(and(aa(:,1)==lidarDataArray(:,1) ,aa(:,2)==lidarDataArray(:,2)), aa(:,3)==lidarDataArray(:,3) )),7);
            seletedIndices = [seletedIndices; idx];
        end
        
        connectComponentIndicesArray{outerLoopCnt} = seletedIndices;
        outerLoopCnt = outerLoopCnt + 1;

    end
    
    retConnectedComp = connectComponentArray;    
    [seedPointsArray, connectComponentArray, connectComponentindexArr] = removeSparseClusters(seedPointsArray, retConnectedComp, connectComponentIndicesArray);
end


function [extremeLiDARDataArray, connectComponentArray, connectComponentindexArr] = removeSparseClusters(extremeLiDARDataArray, connectComponentArray, connectComponentindexArr)

        % Remove branch clusteres with less than 10 points (as they can
        % give misleading branch skeleton)
        idrm = [];
        for j = 1:1:size(connectComponentArray,2)
             ss = connectComponentArray{j};
             if(size(ss,1) > 10)
             %   plot3(ss(:,1), ss(:,2), ss(:,3), '*', 'Color', [rand() rand() rand()]);  
             else
                idrm = [idrm j];
             end
        end        
        extremeLiDARDataArray = removerows(extremeLiDARDataArray,'ind',idrm);
        connectComponentArray = removerows(connectComponentArray','ind',idrm)';
        connectComponentindexArr = removerows(connectComponentindexArr','ind',idrm)';

end


function [extremeLiDARDataArray, extremeLiDARDataIndex, connectComponentArray, connectComponentindexArr] = ...
    getConnectedComponets(stData, seedPointsArray, NumNearestPoints, distTheshold)
    connectComponentArray = []; invalidClusterIndx = [];
    seedPointsArrayIndex = seedPointsArray(:,6);

    cnt = 1;
    for clstStartPointIndex = 1:1:size(seedPointsArray,1)  %length(extremeLiDARDataArray)-50:1:length(extremeLiDARDataArray)-51+noPoints
        
        startPoint = seedPointsArray(clstStartPointIndex,(1:3));
        startPoint = [startPoint seedPointsArray(clstStartPointIndex,(7)) seedPointsArray(clstStartPointIndex,(6))];
        %startPoint = [startPoint stData.lidarDataDensityArr(seedPointsArrayIndex(clstStartPointIndex))];

        otherPoints = stData.lidarDataArray(:,(1:3));
        otherPoints = [otherPoints stData.lidarDataArray(:,7) stData.lidarDataArray(:,6)];
        
        % to set the threshold
        lidstartpoint = startPoint(1:3);
        otherstartpoints = seedPointsArray(:,1:3);        
        [neigdistTheshold, ~]  = pdist2(otherstartpoints, lidstartpoint,'euclidean','Smallest',2);        
        neigdistTheshold = mean(neigdistTheshold(2:2));
        
        connectComponentindexArr{clstStartPointIndex} = [];        
        connectComponentindexArr{clstStartPointIndex} = unique(getProminalPoints(startPoint, otherPoints, neigdistTheshold*10, NumNearestPoints, 30, []));
    
        %hold on;
        %ff = connectComponentindexArr{clstStartPointIndex};
        %plot3(stData.lidarDataArray(ff,1), stData.lidarDataArray(ff,2), stData.lidarDataArray(ff,3),'*','Color',[0.5 0 0]);
        
        idx = connectComponentindexArr{clstStartPointIndex};

        if(length(idx) < 4) % ignore clusters with less than 4 LiDAR points (made constant)
            invalidClusterIndx = [invalidClusterIndx clstStartPointIndex];
        else            
            newnd=[];
            for nidx = 1:1:size(idx,1);
               newnd = [newnd find(stData.lidarDataArray(:,6)==idx(nidx))];
            end
            connectComponentArray{cnt} = stData.lidarDataArray(newnd,1:3);
            cnt = cnt + 1;
        end         
    
    end

    % Remove rows with invalid clusters identified from the above step.
    extremeLiDARDataArray = removerows(seedPointsArray,'ind',invalidClusterIndx);
    extremeLiDARDataIndex = removerows(seedPointsArrayIndex,'ind',invalidClusterIndx);

end

function retPointIndx = getProminalPoints(startPoint, otherPoints, distTheshold, NumNearestPoints, maxIter, retPointIndx)
    [dist,Indx] = pdist2(otherPoints,startPoint,'euclidean','Smallest',NumNearestPoints);
    retPointIndx = [retPointIndx;Indx];
    maxIter = maxIter-1;
    for i = 2:1:2
        if(exp(dist(i)) < distTheshold) && maxIter > 0
            retPointIndx = [retPointIndx;getProminalPoints(otherPoints(retPointIndx(size(retPointIndx,1)),:), otherPoints, distTheshold, NumNearestPoints, maxIter, retPointIndx)];
        else
            retPointIndx = [];
        end
    end  
end

function retPointIndx = getProminalPointsNew(startPoint, otherPoints, distTheshold, NumNearestPoints, maxIter, retPointIndx)
    [dist,Indx] = pdist2(otherPoints(:,1:4),startPoint(1,1:4),'euclidean','Smallest',NumNearestPoints);
    retPointIndx = [retPointIndx; otherPoints(Indx, 5)];
    otherPoints  = removerows(otherPoints,'ind',Indx(1:end-1));
    maxIter = maxIter-1;
    for i = 2:1:2
        if(dist(i) < distTheshold) && maxIter > 0
            aa =  find(otherPoints(:,5)==retPointIndx(end));
            retPointIndx = [retPointIndx;getProminalPointsNew(otherPoints(aa,:), otherPoints, distTheshold, NumNearestPoints, maxIter, retPointIndx)];
        else
            retPointIndx = [];
        end
    end  
end



function retDistWeights = calculateDistWeights(lidarDataArray, neighbourhoodSize)
   [nearPoint,~] = pdist2(lidarDataArray(:,(1:3)),lidarDataArray(:,(1:3)),'euclidean','Smallest',25);  
    nearPoint(nearPoint>=neighbourhoodSize) = 0; % make dist 0 for points > threshold  neighbourhoodDist  
    nearPoint = nearPoint(2:end,:); % remove first line because it is distance to point itself.
    neigDensityArr = sum(nearPoint~=0,1)'; % find non zero elenents in the array for every column
    retDistWeights = neigDensityArr/norm(neigDensityArr); % normalization
end

function plotConnectedComponents(connectComponentArray)
    for i = 1:1:size(connectComponentArray,3);
        clusteredLiDARDataArray = connectComponentArray(:,:,i);
        plot3(clusteredLiDARDataArray(:,1), clusteredLiDARDataArray(:,2), clusteredLiDARDataArray(:,3),'*','Color',[0.5 0 0]);
        hold on; camproj perspective; rotate3d on; view(3), axis vis3d; axis equal; axis on;
    end
end

function [stemPointArr,D1] = getStemPointFromPCA(extremeLiDARDataArray, connectComponentArray,retMaxXY)
    stemPointArr = zeros(size(connectComponentArray,2),3);
    slopeValArr = [];
    D1={};
    for i = 1:1:size(connectComponentArray,2);
        % Get modified connected components, after adding prospestive stem
        [stemPointArr(i,:), ~, D1{i}] = getModifiedConnectedComp(connectComponentArray{i},extremeLiDARDataArray(i,(1:3)),retMaxXY);        
    end
end

% Plot the line 
function [retSlope, lineLength, pntDistToline, maxeval, eigRatio, totalPoints] = getBestFitLine(stData, dataArray, compStartPoint, compEndPoint, slopeThreshold, speciesCode, isPlotON, plotOn, colorArr)
    
    % Fit a line throught the cluster of points (to get branch skeleton)
    [m,p,pntDistToline] = best_fit_line(dataArray(:,1),dataArray(:,2), dataArray(:,3));
    % Calcuate the slope of the line (in Degrees)
    retSlope = subspace(p',[p(1:2) 0]')*(180/pi);

    % find the nearest point on the stem where the branch lines crosses
    % as as line start point
    stemBottom =[stData.retMaxXYZ(:,1:2) 0];
    stemTop =[stData.retMaxXYZ(:,1:2) 100];        
    P3 = [m(1)+p(1)*100 m(2)+p(2)*100 m(3)+p(3)*100];
    P4 = [m(1)+p(1)*-100 m(2)+p(2)*-100 m(3)+p(3)*-100];        
    [distToStem, nearestStemPoints] = distBtwLineSegments(stemBottom, stemTop, P3, P4);       
    %pt = nearestStemPoints{2}; 

    % find the extreme point in cluster as line end point
    [~, Dindx] = max(sqrt(sum((repmat(compStartPoint,size(dataArray,1),1) - dataArray).^ 2,2)));
    endpnt = dataArray(Dindx,:);
     %endpnt = compEndPoint; % if branch tips are sure to be end points

    %hold on;
    %plot3(pt(1),pt(2),pt(3),'*','Color','r', 'MarkerSize', 18);

    % distance from cluster centre to branch start point(i.e. near trunk)
    dst1 = pdist2(m(1,:),compStartPoint,'euclidean','Smallest',1);    
    % distnace of Cluster centre from the branch exterior points.
    dst2 = pdist2(endpnt,m(1,:),'euclidean','Smallest',1);


    % Considering that points are located at different directions around the stem.
    if(endpnt(1,1)>0)
        t=-dst2:0.2:dst1;
    else
        t= -dst1:0.2:dst2;
    end

    x1=m(1)+p(1)*t; y1=m(2)+p(2)*t; z1=m(3)+p(3)*t;

    % jugaad fix (need to correct)
    DI = max(sqrt(sum((repmat(compStartPoint,size(x1,2),1) - [x1' y1' z1']).^ 2,2)));  
    DIend = max(sqrt(sum((compStartPoint - endpnt).^ 2,2)));        
    if(DI>DIend+0.3)
        t =-t; x1=m(1)+p(1)*t; y1=m(2)+p(2)*t; z1=m(3)+p(3)*t;            
    end
    
    if(~strcmp(speciesCode,'pc'))
        trsh = 5;        
    else
        trsh =2;
    end

    % Plot line
    if(abs(retSlope) < slopeThreshold && compStartPoint(3) > trsh) 
        if(and(isPlotON,plotOn))
            %hold on
           % plot3(compStartPoint(1),compStartPoint(2),compStartPoint(3),'.','Color',[0 1 0],'MarkerSize',40);
           % plot3(endpnt(1),endpnt(2),endpnt(3),'.','Color',[1 0 0],'MarkerSize',40);
            %plot3(m(1),m(2),m(3),'.','Color',[0 0 1],'MarkerSize',40);
            lnData = [compStartPoint; endpnt];
            p1 = line(lnData(:,1),lnData(:,2),lnData(:,3), 'Color',[139/256 69/256 19/256],'LineWidth', 2); alpha(0.9);
            p1.Color(4) = 0.7;
            %plot3(x1,y1,z1,'LineWidth',3, 'Color',colorArr);
            %hold off;
        end
    end

    lineLength = (dst1 + dst2);

    % Calculate Eigen features
    retEigenData = getBestEigenDirection(dataArray);
    evec = retEigenData{1};
    eval = retEigenData{2};

    % handle cases where only one/two PC is available.
    maxeval = 0; eigRatio = 1;
    if(length(eval)==3)
        maxeval = max(eval(2:end));
        eigRatio = eval(2)/eval(3);
    elseif(length(eval)==2)
        maxeval = -1; %max(eval(2:end));
        eigRatio = -1;
    else            
        maxeval = -1; eigRatio =-1;
    end

    % total no of LiDAR point in the branch cluster
    totalPoints = size(dataArray,1);
        
end

function dataAtt = dataPreProcessing(singleTreeLiDARdata)
    % Write normalized data to a table for performance improvement 
    lidarDataArr = normalizeLiDARData(write2table(singleTreeLiDARdata));
    %lidarDataArray = lidarDataArray(randperm(size(lidarDataArray,1),9000),:);

    treeWidth =  max(lidarDataArr(:,1)); treeHeight = max(lidarDataArr(:,2)); % to calculate max wisth/breadth to set plot width
    dataAtt.maxTreeWidth = max(treeWidth, treeHeight) + 0.5; % Keep it as 5 if issues arise
    dataAtt.mintreeHeight = min(lidarDataArr(:,3));
    dataAtt.maxtreeHeight = max(lidarDataArr(:,3));        
    dataAtt.treeHeight = dataAtt.maxtreeHeight - dataAtt.mintreeHeight;        
    dataAtt.htDeduction = dataAtt.treeHeight*0.05; % to get rid of ground noise points

    % Crown height
    lidarDataArr = lidarDataArr(and(lidarDataArr(:,3) > min(lidarDataArr(:,3))+ dataAtt.htDeduction, lidarDataArr(:,3) < max(lidarDataArr(:,3))),:);
    dataAtt.minCrownHeight = getMinCrownHeight(lidarDataArr);
    dataAtt.maxCrownHeight = max(lidarDataArr(:,3));        
    dataAtt.crownHeight = dataAtt.maxCrownHeight - dataAtt.minCrownHeight;       

    %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 > 0.1),:);                
    %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 < 12),:);

    % Identify the center point of the tree (in top view).
    dataAtt.retMaxXYZ = findMaxHeightXY(lidarDataArr);

    % Original LiDAR data Array
    index = 1:1:size(lidarDataArr,1);
    lidarDataDensityArr =  calculateDistWeights(lidarDataArr,0.25);
    lidarDataArr = [lidarDataArr index' lidarDataDensityArr];
    dataAtt.lidarDataArrayComplete =  lidarDataArr(:,1:8);
    
    % Sparsified LiDAR data Array
    lidarDataArrSampled = lidarDataArr(lidarDataDensityArr(:,1) >= 0.0010,:);    
    dataAtt.lidarDataArray =  lidarDataArrSampled(:,1:8);                          
end


function eigDetails = getEigVetStemProximities(dataArr, extremeLiDARDataPoint, retMaxXY, plotPCAxisOn)
   retEigenData = getBestEigenDirection(dataArr);
   retEigPoints = retEigenData{1}; 
   retEigvalues = retEigenData{2};
          
   stemBottom = [retMaxXY(1);retMaxXY(2);0];
   stemTop = [retMaxXY(1);retMaxXY(2);retMaxXY(3)];
    
   eigDetails.distToStematCrossingArr = []; eigDetails.nearestStemPointsArr = []; eigDetails.crossingDistanceArr = [];
   eigDetails.pcVector =[];  eigDetails.pcValue =[];
   
   for ivCnt =1:1:size(retEigPoints,2)
        pcVector = retEigPoints(:,ivCnt); eigDetails.pcVector = [eigDetails.pcVector pcVector];
        pcValue = retEigvalues(ivCnt); eigDetails.pcValue = [eigDetails.pcValue pcValue];
        
        P3 = pcVector *10; % 10 is used as a large value as there exits no trees with width > 10 m in dataset
        P3(1)=P3(1)+extremeLiDARDataPoint(1); P3(2)=P3(2)+extremeLiDARDataPoint(2); P3(3) = P3(3)+extremeLiDARDataPoint(3);
        P4 = pcVector *-10; % 10 is used as a large value as there exits no trees with width > 10 m
        P4(1)=P4(1)+extremeLiDARDataPoint(1); P4(2)=P4(2)+extremeLiDARDataPoint(2); P4(3) = P4(3)+extremeLiDARDataPoint(3);
        
        if(plotPCAxisOn && ivCnt ==1)
            vect2= [P3';P4'];  
            hold on; line(vect2(:,1),vect2(:,2),vect2(:,3),'LineWidth' ,3);
        end
        
        % The function provides what is the mimimum distance between the stem line and the line formed using the PC vector direction.
        [distToStem,nearestStemPoints] = distBtwLineSegments(stemBottom, stemTop, P3, P4);
        crossingDistance = norm(nearestStemPoints{2}' - extremeLiDARDataPoint);
        
        % Store details of each of the PC vectors in array
        eigDetails.distToStematCrossingArr = [eigDetails.distToStematCrossingArr distToStem]; 
        eigDetails.nearestStemPointsArr = [eigDetails.nearestStemPointsArr; nearestStemPoints{2}'];
        eigDetails.crossingDistanceArr = [eigDetails.crossingDistanceArr crossingDistance];
   end
     
end


function [point, PricComp1Slope, D1] = getModifiedConnectedComp(dataArr,extremeLiDARDataPoint,retMaxXY)
        
      % find the closest eigen vector pointing towards the stem
        eigDetailsA = getEigVetStemProximities(dataArr, extremeLiDARDataPoint, retMaxXY, false);
      
        dx = eigDetailsA.distToStematCrossingArr;
        nsa = eigDetailsA.nearestStemPointsArr;
      %  retEigPoints = eigDetailsA.pcVector;
       % evalj = eigDetailsA.pcValue;
        D1 = eigDetailsA.crossingDistanceArr;
       %point =nsa(2,:);
          
        DL = [dx(1)*D1(1) dx(2)*D1(2) dx(3)*D1(3)];
        [maxD,maxDi] = min(DL);
        point = 0;
        
        slopeA =[];
        if(dx(1)*D1(1)==maxD)
            point = nsa(1,:);
        elseif(dx(2)*D1(2)==maxD)
            point = nsa(2,:);
        else
            point = nsa(3,:);
        end
     
        point1 = point/norm(point);
        PricComp1Slope = subspace(point1',[point1(1:2) mean(dataArr(:,3))]')*(180/pi);
 
         %point =  nsa(1,:);

%         c=[]; e=[];
%         for v =1:1:size(retEigPoints,2)
%             a = retEigPoints(:,v);
%             b = [0;0;1];
%             d = [0;0;-1];
%             c = [c acos((a'*b)/(norm(a)*norm(b)))*(180/pi)];
%             e = [e acos((a'*d)/(norm(a)*norm(d)))*(180/pi)];                
%             f= c-e;
%         end

%         [valc, indc] =  min(abs(f));
% 
%          if(abs(f(indc))>2000000)
%              DL = [D1(1) D1(2) D1(3)];
%              [maxD] = min(DL);
%              point = 0;

end

function lidarDataArray = normalizeLiDARData(lidarDataArray)
    midxPoint = (max(lidarDataArray(:,1))- min(lidarDataArray(:,1)))/2;
    midyPoint = (max(lidarDataArray(:,2))- min(lidarDataArray(:,2)))/2;
    
    lidarDataArray(:,1) = lidarDataArray(:,1)- min(lidarDataArray(:,1));
    lidarDataArray(:,2) = lidarDataArray(:,2)- min(lidarDataArray(:,2));
    lidarDataArray(:,3) = lidarDataArray(:,3)- min(lidarDataArray(:,3));

    lidarDataArray(:,1) = lidarDataArray(:,1)- midxPoint;
    lidarDataArray(:,2) = lidarDataArray(:,2)- midyPoint;
end

function lidarBranchTips = removeRepetitions(stData, extremeLiDARDaraArray, CutoffDistBetwVect)
dstarr =  pdist2(extremeLiDARDaraArray(:,1:2), stData.retMaxXYZ(1:2));
extremeLiDARDaraArray = [extremeLiDARDaraArray  dstarr];
extremeLiDARDaraArray = sortrows(extremeLiDARDaraArray,[3 8]);
    
    lidarBranchTips =[];
    
    for j=1:1:max(extremeLiDARDaraArray(:,3))       
        ss = CutoffDistBetwVect*((j*1)/max(extremeLiDARDaraArray(:,3)));
        CutoffDistBetwVect1 = CutoffDistBetwVect - ss;  %ar =1.2; la = 1.8; % CutoffDistBetwVect - 
        extremeLiDARDdataBysection = extremeLiDARDaraArray(find(and(extremeLiDARDaraArray(:,3) >= j, extremeLiDARDaraArray(:,3) <  j + 1)),:);
        
        if(size(extremeLiDARDdataBysection,1) > 0)
            DistBwExtrmPoints = tril(squareform(pdist(extremeLiDARDdataBysection(:,1:3),'euclidean')));        
            [rowIdx,~] = ind2sub(size(DistBwExtrmPoints),find(DistBwExtrmPoints < CutoffDistBetwVect1 & DistBwExtrmPoints > 0));
            extremeLiDARDdataBysection = removerows(extremeLiDARDdataBysection,'ind',unique(rowIdx));     
            lidarBranchTips = [lidarBranchTips; extremeLiDARDdataBysection];
        end
    end
  %  lidarBranchTips = removerows(extremeLiDARDaraArray,'ind',indx);
end


function [m,p,s] = best_fit_line(x,y,z)
    % x,y,z are n x 1 column vectors of the three coordinates
    % of a set of n points in three dimensions. The best line,
    % in the minimum mean square orthogonal distance sense,
    % will pass through m and have direction cosines in p, so
    % it can be expressed parametrically as x = m(1) + p(1)*t,
    % y = m(2) + p(2)*t, and z = m(3)+p(3)*t, where t is the
    % distance along the line from the mean point at m.
    % s returns with the minimum mean square orthogonal
    % distance to the line.
    % RAS - March 14, 2005

    [n,mx] = size(x); [ny,my] = size(y); [nz,mz] = size(z);
    if (mx~=1)|(my~=1)|(mz~=1)|(ny~=n)|(nz~=n)
     error('The arguments must be column vectors of the same length.')
    end
    m = [mean(x),mean(y),mean(z)];
    w = [x-m(1),y-m(2),z-m(3)]; % Use "mean" point as base
    a = (1/n)*w'*w; % 'a' is a positive definite matrix
    [u,d,v] = svd(a); % 'eig' & 'svd' get same eigenvalues for this matrix
    p = u(:,1)'; % Get eigenvector for largest eigenvalue
    s = d(2,2)+d(3,3); % Sum the other two eigenvalues
end

function retVector = rotateVector(inputVector, teeta)
    rotationMatrix = [cos(teeta) -sin(teeta) 0; sin(teeta) cos(teeta) 0; 0 0 1];
    retVector = rotationMatrix*inputVector';
end

function plotLiDARData(lidarDataArray, gridOnOffParameterCylin,  gridOnOffParameter,htDeduction,showRawLiDARParameter,zDiv,retMaxXYZ)
    if(showRawLiDARParameter==true)
        plot3(lidarDataArray(:,1), lidarDataArray(:,2), lidarDataArray(:,3),'.','MarkerSize',5,'Color',[0 0.5 0]);
    end
    
    if(gridOnOffParameterCylin)
       
        r=max(max(lidarDataArray(:,1:2)))+1;
        maxZ = max(lidarDataArray(:,3)); minZ = min(lidarDataArray(:,3));
        xamples =[];  yamples =[]; %zDiv = 15; 
        for zStep = minZ :(maxZ - minZ)/zDiv : maxZ;
            teta=-pi:0.314/4:pi;
            x=r*cos(teta) + retMaxXYZ(1);
            y=r*sin(teta) + retMaxXYZ(2);
            plot3(x,y,zeros(1,numel(x))+zStep,'Color','k', 'LineWidth',1);
            plot3(x/2,y/2,zeros(1,numel(x))+zStep,'Color','k', 'LineWidth',1);
            plot3(x/4,y/4,zeros(1,numel(x))+zStep,'Color','k', 'LineWidth',1);
            
            xamples = x(1:4:size(x,2));
            yamples = y(1:4:size(y,2));
            plot3(xamples,yamples,zeros(1,size(xamples,2))+zStep,'.','Color','k');
            
            centMat = [repmat(retMaxXYZ(:,1:2),size(xamples,2)*2,1) zeros(1,size(xamples,2)*2)'+zStep ];
            centMat([1:2:size(centMat,1)],:) = [xamples' yamples' zeros(1,size(xamples,2))'+zStep];
            
            line(centMat(:,1), centMat(:,2), centMat(:,3),'Color','k','LineWidth',1);
        end
        
        for i = 1:1:size(xamples,2)      
            plot3([xamples(i) xamples(i)],[yamples(i) yamples(i)],[6.2 maxZ], '-', 'Color', 'k', 'LineWidth',1);        
        end
        
    end
               
    if(gridOnOffParameter)
        
        maxX = max(lidarDataArray(:,1)); minX = min(lidarDataArray(:,1));
        maxY = max(lidarDataArray(:,2)); minY = min(lidarDataArray(:,2));
        maxZ = max(lidarDataArray(:,3)); minZ = min(lidarDataArray(:,3))-htDeduction;
        %xDiv = 10; yDiv = 10; zDiv = 25; 

        x = minX:(maxX - minX)/zDiv: maxX;
        y = minY:(maxY - minY)/zDiv: maxY;
        z = minZ:(maxZ - minZ)/zDiv: maxZ;

        [X1 Y1 Z1] = meshgrid(x([1 end]),y,z);
        X1 = permute(X1,[2 1 3]); Y1 = permute(Y1,[2 1 3]); Z1 = permute(Z1,[2 1 3]);
        X1(end+1,:,:) = NaN; Y1(end+1,:,:) = NaN; Z1(end+1,:,:) = NaN;

        [X2 Y2 Z2] = meshgrid(x,y([1 end]),z);
        X2(end+1,:,:) = NaN; Y2(end+1,:,:) = NaN; Z2(end+1,:,:) = NaN;

        [X3 Y3 Z3] = meshgrid(x,y,z([1 end]));   
        X3 = permute(X3,[3 1 2]); Y3 = permute(Y3,[3 1 2]); Z3 = permute(Z3,[3 1 2]);
        X3(end+1,:,:) = NaN; Y3(end+1,:,:) = NaN; Z3(end+1,:,:) = NaN;

      %#figure('Renderer','opengl')
        h = line([X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)], [Z1(:);Z2(:);Z3(:)]);
        set(h, 'Color',[10/256,10/256,10/256], 'LineWidth',0.5, 'LineStyle','-')
        
    end
   % hold off;
end

function returnData =  LoadData(fullFileName)
    returnData = lasdata(fullFileName);
end

function retTable = write2table(lasFile)
    retTable = zeros(size(lasFile.x,1),6);
    retTable(:,1) = lasFile.x;
    retTable(:,2) = lasFile.y;
    retTable(:,3) = lasFile.z;
    retTable(:,4) = get_classification(lasFile);    
    retTable(:,5) = lasFile.get_return_number;
    retTable(:,6) = get_intensity(lasFile);    
end

function [distance, dd] = distBtwLineSegments(p1, p2, p3, p4)

    u = p1 - p2; 
    v = p3 - p4;
    w = p2 - p4;
    
    a = dot(u,u);
    b = dot(u,v);
    c = dot(v,v);
    d = dot(u,w);
    e = dot(v,w);
    D = a*c - b*b;
    sD = D;
    tD = D;
    
    SMALL_NUM = 0.00000001;
    
    % compute the line parameters of the two closest points
    if (D < SMALL_NUM)  % the lines are almost parallel
        sN = 0.0;       % force using point P0 on segment S1
        sD = 1.0;       % to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    else                % get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0)   % sc < 0 => the s=0 edge is visible       
            sN = 0.0;
            tN = e;
            tD = c;
        elseif (sN > sD)% sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        end
    end
    
    if (tN < 0.0)            % tc < 0 => the t=0 edge is visible
        tN = 0.0;
        % recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        elseif (-d > a)
            sN = sD;
        else
            sN = -d;
            sD = a;
        end
    elseif (tN > tD)       % tc > 1 => the t=1 edge is visible
        tN = tD;
        % recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        elseif ((-d + b) > a)
            sN = sD;
        else 
            sN = (-d + b);
            sD = a;
        end
    end
    
    % finally do the division to get sc and tc
    if(abs(sN) < SMALL_NUM)
        sc = 0.0;
    else
        sc = sN / sD;
    end
    
    if(abs(tN) < SMALL_NUM)
        tc = 0.0;
    else
        tc = tN / tD;
    end
    
    % get the difference of the two closest points
    dP = w + (sc * u) - (tc * v);  % = S1(sc) - S2(tc)

    distance = norm(dP);
    outV = dP;
    
    varargout(1) = {outV};      % vector connecting the closest points
    varargout(2) = {p2+sc*u};   % Closest point on object 1 
    varargout(3) = {p4+tc*v};   % Closest point on object 2
    
    dd = varargout;
    
end

function crownHt = getMinCrownHeight(lidarDataArray)
   maxTreeHeight = max(lidarDataArray(:,3));
   crownHt = min(lidarDataArray(:,3));
   for i = maxTreeHeight:-1:1
       templidararray = lidarDataArray(find(and(lidarDataArray(:,3)<i,(lidarDataArray(:,3)>i-1))),1:3);
        if(size(templidararray,1)<10)
            crownHt = i;
            break;
        end
   end
end

function retMaxXYZ = findMaxHeightXY(lidarDataArray)
    [~,indx] = max(lidarDataArray(:,3));
     retMaxXYZ = lidarDataArray(indx,1:3);
end

function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)

    % We can solve for T with T = A\B.

    % Input checks.
    error(nargchk(2,5,nargin))

    % Adjustments when fewer than five arguments are supplied.
    switch nargin
        case 2
            robust = true;
            x2 = x1;
            y2 = y1;
            self_intersect = true;
        case 3
            robust = x2;
            x2 = x1;
            y2 = y1;
            self_intersect = true;
        case 4
            robust = true;
            self_intersect = false;
        case 5
            self_intersect = false;
    end

    % x1 and y1 must be vectors with same number of points (at least 2).
    if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
            length(x1) ~= length(y1)
        error('X1 and Y1 must be equal-length vectors of at least 2 points.')
    end
    % x2 and y2 must be vectors with same number of points (at least 2).
    if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
            length(x2) ~= length(y2)
        error('X2 and Y2 must be equal-length vectors of at least 2 points.')
    end


    % Force all inputs to be column vectors.
    x1 = x1(:);
    y1 = y1(:);
    x2 = x2(:);
    y2 = y2(:);

    % Compute number of line segments in each curve and some differences we'll
    % need later.
    n1 = length(x1) - 1;
    n2 = length(x2) - 1;
    xy1 = [x1 y1];
    xy2 = [x2 y2];
    dxy1 = diff(xy1);
    dxy2 = diff(xy2);

    % Determine the combinations of i and j where the rectangle enclosing the
    % i'th line segment of curve 1 overlaps with the rectangle enclosing the
    % j'th line segment of curve 2.
    [i,j] = find(repmat(min(x1(1:end-1),x1(2:end)),1,n2) <= ...
        repmat(max(x2(1:end-1),x2(2:end)).',n1,1) & ...
        repmat(max(x1(1:end-1),x1(2:end)),1,n2) >= ...
        repmat(min(x2(1:end-1),x2(2:end)).',n1,1) & ...
        repmat(min(y1(1:end-1),y1(2:end)),1,n2) <= ...
        repmat(max(y2(1:end-1),y2(2:end)).',n1,1) & ...
        repmat(max(y1(1:end-1),y1(2:end)),1,n2) >= ...
        repmat(min(y2(1:end-1),y2(2:end)).',n1,1));

    % Force i and j to be column vectors, even when their length is zero, i.e.,
    % we want them to be 0-by-1 instead of 0-by-0.
    i = reshape(i,[],1);
    j = reshape(j,[],1);

    % Find segments pairs which have at least one vertex = NaN and remove them.
    % This line is a fast way of finding such segment pairs.  We take
    % advantage of the fact that NaNs propagate through calculations, in
    % particular subtraction (in the calculation of dxy1 and dxy2, which we
    % need anyway) and addition.
    % At the same time we can remove redundant combinations of i and j in the
    % case of finding intersections of a line with itself.
    if self_intersect
        remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
    else
        remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
    end
    i(remove) = [];
    j(remove) = [];

    % Initialize matrices.  We'll put the T's and B's in matrices and use them
    % one column at a time.  AA is a 3-D extension of A where we'll use one
    % plane at a time.
    n = length(i);
    T = zeros(4,n);
    AA = zeros(4,4,n);
    AA([1 2],3,:) = -1;
    AA([3 4],4,:) = -1;
    AA([1 3],1,:) = dxy1(i,:).';
    AA([2 4],2,:) = dxy2(j,:).';
    B = -[x1(i) x2(j) y1(i) y2(j)].';

    % Loop through possibilities.  Trap singularity warning and then use
    % lastwarn to see if that plane of AA is near singular.  Process any such
    % segment pairs to determine if they are colinear (overlap) or merely
    % parallel.  That test consists of checking to see if one of the endpoints
    % of the curve 2 segment lies on the curve 1 segment.  This is done by
    % checking the cross product
    %
    %   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
    %
    % If this is close to zero then the segments overlap.

    % If the robust option is false then we assume no two segment pairs are
    % parallel and just go ahead and do the computation.  If A is ever singular
    % a warning will appear.  This is faster and obviously you should use it
    % only when you know you will never have overlapping or parallel segment
    % pairs.

    if robust
        overlap = false(n,1);
        warning_state = warning('off','MATLAB:singularMatrix');
        % Use try-catch to guarantee original warning state is restored.
        try
            lastwarn('')
            for k = 1:n
                T(:,k) = AA(:,:,k)\B(:,k);
                [unused,last_warn] = lastwarn;
                lastwarn('')
                if strcmp(last_warn,'MATLAB:singularMatrix')
                    % Force in_range(k) to be false.
                    T(1,k) = NaN;
                    % Determine if these segments overlap or are just parallel.
                    overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
                end
            end
            warning(warning_state)
        catch err
            warning(warning_state)
            rethrow(err)
        end
        % Find where t1 and t2 are between 0 and 1 and return the corresponding
        % x0 and y0 values.
        in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
        % For overlapping segment pairs the algorithm will return an
        % intersection point that is at the center of the overlapping region.
        if any(overlap)
            ia = i(overlap);
            ja = j(overlap);
            % set x0 and y0 to middle of overlapping region.
            T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
                min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
            T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
                min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
            selected = in_range | overlap;
        else
            selected = in_range;
        end
        xy0 = T(3:4,selected).';

        % Remove duplicate intersection points.
        [xy0,index] = unique(xy0,'rows');
        x0 = xy0(:,1);
        y0 = xy0(:,2);

        % Compute how far along each line segment the intersections are.
        if nargout > 2
            sel_index = find(selected);
            sel = sel_index(index);
            iout = i(sel) + T(1,sel).';
            jout = j(sel) + T(2,sel).';
        end
    else % non-robust option
        for k = 1:n
            [L,U] = lu(AA(:,:,k));
            T(:,k) = U\(L\B(:,k));
        end

        % Find where t1 and t2 are between 0 and 1 and return the corresponding
        % x0 and y0 values.
        in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
        x0 = T(3,in_range).';
        y0 = T(4,in_range).';

        % Compute how far along each line segment the intersections are.
        if nargout > 2
            iout = i(in_range) + T(1,in_range).';
            jout = j(in_range) + T(2,in_range).';
        end
    end
end

% function [extremeLiDARDataArray, extremeLiDARDataIndices] = locateExteriorPointsEfficient(stData, zDiv, maxAngle, angleIncrement, removeSimilarRows, neighbourCutoff)
%     %extremeLiDARDataArray = zeros(zDiv*maxAngle/angleIncrement,size(lidarDataArray,2)); %cell2table(cell(0,5));
%     extremeLiDARDataArray =[]; maxZ = max(stData.lidarDataArray(:,3)); minZ = min(stData.lidarDataArray(:,3));
%     extremeLiDARDataIndices = [];
%     for zStep = minZ :(maxZ - minZ)/zDiv : maxZ-(maxZ - minZ)/zDiv;
%             clippedLiDARpointsIndices = find(and(stData.lidarDataArray(:,3) >= zStep , stData.lidarDataArray(:,3) < zStep+(maxZ - minZ)/zDiv));
%             data = [stData.lidarDataArray(clippedLiDARpointsIndices,1:3) stData.index(clippedLiDARpointsIndices)'];    
%             
%             for angle = 0: angleIncrement: maxAngle-angleIncrement;
%                 ang1 = angle; ang2 = angle+angleIncrement;
%                 r = stData.maxTreeWidth;
%                 v1 = [r*cos(ang1*pi/180) r*sin(ang1*pi/180)]; v2 = [r*cos(ang2*pi/180) r*sin(ang2*pi/180)];
%                 isPointWithin = isWinthinArc(v1, v2, data(:,1:2), r);             
%                 
%                 P2 = data(isPointWithin,1:4); 
%                 P1 = [ones(size(P2,1),1)*stData.retMaxXYZ(1) ones(size(P2,1),1)*stData.retMaxXYZ(2) data(isPointWithin,3)];
%                
%                 if(size(P1,1)>0)
%                     [~, Dindx] = max(sqrt(sum((P1 - P2(:,1:3)).^ 2,2)));
%                     indxx = P2(Dindx,4);
%                     extremeLiDARDataIndices = [extremeLiDARDataIndices indxx];
%                 end
%                
%             end
%     end
%     
%     extremeLiDARDataIndices = unique(extremeLiDARDataIndices); 
%     edr = stData.lidarDataArray(extremeLiDARDataIndices,1:3);
%     
%     extremeLiDARDaraArray = edr;    
%     if(removeSimilarRows)
%         extremeLiDARDaraArray = removeRepetitions(edr, neighbourCutoff);
%     end
%     
%     k = boundary(extremeLiDARDaraArray(:,1),extremeLiDARDaraArray(:,2),extremeLiDARDaraArray(:,3), 0.8);
%     k =unique(k(:));
%     extremeLiDARDataArray = extremeLiDARDaraArray(k,:);
%     extremeLiDARDataIndices = extremeLiDARDataIndices(k);
%     
%     %plot3(extremeLiDARDaraArray(:,1),extremeLiDARDaraArray(:,2),extremeLiDARDaraArray(:,3),'.','Color','r','MarkerSize',20);    
%     %plot3(stData.lidarDataArray(:,1),stData.lidarDataArray(:,2),stData.lidarDataArray(:,3),'.','Color','k');    
% end


% function lidarBranchTips = locateExteriorPoints1(lidarDataArray, maxXYZ, xDiv, yDiv, zDiv)
%    
%    maxXYZs =  [max(lidarDataArray(:,1)) max(lidarDataArray(:,2)) max(lidarDataArray(:,3))];
%    minXYZs =  [min(lidarDataArray(:,1)) min(lidarDataArray(:,2)) min(lidarDataArray(:,3))];
% 
%    zDiv = zDivi; xDiv = xDivi; yDiv = yDivi; xIncr = (maxXYZs(1)-minXYZs(1))/xDiv; yIncr = (maxXYZs(2)-minXYZs(2))/yDiv;
%    
%    maxXinSection = []; maxYinSection = []; zcnt =1; zcntArr=[];
%    tipPointIndxArray = []; sectionXLimitValArr = []; sectionYLimitValArr = []; sectionXLimitValIncrArr = []; sectionYLimitValIncrArr = [];
%    for zStep = maxXYZs(3) : -(maxXYZs(3) - minXYZs(3))/zDiv : minXYZs(3)
%        
%         clippedPointIndices = find(and(lidarDataArray(:,3) >= zStep, lidarDataArray(:,3) < zStep+(maxXYZs(3) - minXYZs(3))/zDiv));        
%         data = [lidarDataArray(clippedPointIndices,1:3) clippedPointIndices];
%         
%          maxXinSection = [maxXinSection max(data(:,1))];
%          %minXinSection = min(data(:,1)); divXSect = (maxXinSection-minXinSection)/zDiv;        
%          maxYinSection = [maxYinSection max(data(:,2))];
%          %minYinSection = min(data(:,2)); divYSect = (maxYinSection-minYinSection)/zDiv;         
% %         for xSect = minXinSection: divXSect : maxXinSection
% %             for ySect = minYinSection : divYSect : maxYinSection
% %                 
% %             end            
% %         end
%        for xSect = minXYZs(1): xIncr: maxXYZs(1)
%            for ySect = minXYZs(2): yIncr : maxXYZs(2)
%                %cc = ySect
%                %ySect = ySect;
%                boxDataIndices = find( and(  and(  and(data(:,1) >= xSect, data(:,1) < (xSect+xIncr) ), data(:,2) >= ySect ), data(:,2) < ySect+yIncr));
%                dataBox  = data(boxDataIndices,:);
%                
%                P1 = dataBox(:,1:3);
%                P2 = [zeros(size(dataBox,1),1) zeros(size(dataBox,1),1) dataBox(:,3)];
%                
%                if(size(P1,1)>0)
%                     [dmax, Dindx] = max(sqrt(sum((P1 - P2).^ 2,2)));               
%                     tipPointIndxArray = [tipPointIndxArray dataBox(Dindx,4)];
%                     sectionXLimitValArr = [sectionXLimitValArr maxXinSection];
%                     sectionYLimitValArr= [sectionYLimitValArr maxYinSection];
%                     sectionXLimitValIncrArr = [sectionXLimitValIncrArr xIncr];
%                     sectionYLimitValIncrArr = [sectionYLimitValIncrArr yIncr];
%                     zcntArr = [zcntArr zcnt];
%                end
%            end            
%        end
%        zcnt = zcnt+1;
%    end
%    
%    [tipPointIndxArray, idx, icx] = unique(tipPointIndxArray);
%    tipPointArray = [lidarDataArray(tipPointIndxArray,(1:3)) tipPointIndxArray' sectionXLimitValArr(idx)' sectionYLimitValArr(idx)' sectionXLimitValIncrArr(idx)' sectionYLimitValIncrArr(idx)' zcntArr'];
%    
%    pnts = [];
%    for im = 1: 1: size(unique(tipPointArray(:,9)),1)
%        tmPPnts = tipPointArray(tipPointArray(:,9)==im,1:3);
%        k = boundary(tmPPnts(:,1),tmPPnts(:,2),0.1);
%        pnts = [pnts;tmPPnts(k,1:3)];
%    %    if(and( and( and(abs(tipPointArray(im,1))<=abs(tipPointArray(im,5)), abs(tipPointArray(im,1))>(abs(tipPointArray(im,5))-abs(tipPointArray(im,7)))), abs(tipPointArray(im,2)) <= abs(tipPointArray(im,6)) ), abs(tipPointArray(im,2)) > (abs(tipPointArray(im,6))-abs(tipPointArray(im,8))) ))
%    %         id = [id im];
%    %    end
%    end
%    
%    %tipPointArray = tipPointArray(id',:);
%    
%    
%    lidarBranchTips = removeRepetitions1(pnts, 1);   %1:ar  1:la  0.5:pc   0.5:ab
%   % lidarBranchTips = removeRepetitions1(pnts, 1);   
%    
%    
% end

% function [retConnectedComp, seletedIndices] = getConnectedComponetsold(lidarDataArray, seedPointsArray, compStartIndex, compEndIndex)
% 
%     outerLoopCnt= 1;
%     noOfIterations = 12; NumNearestPoints = 2;
%     connectComponentArray = {};%zeros(NumNearestPoints*noOfIterations,3,(compEndIndex-compStartIndex));
%     
%     for clstStartPointIndex = compStartIndex:1:compEndIndex  %length(extremeLiDARDataArray)-50:1:length(extremeLiDARDataArray)-51+noPoints
%         cntd =1; clusteredLiDARDataArray = zeros(NumNearestPoints*noOfIterations,3); 
%         startPoint = seedPointsArray(clstStartPointIndex,(1:3));
%         seletedIndices = [];
%          % [double(startPoint(1)) double(startPoint(2)) double(startPoint(3))]
%         otherPoints = lidarDataArray(:,(1:3));
%         for loopCnt = 1:1:noOfIterations
%             %(find(and(lidarDataArray(:,3) > startPoint(3)-10, lidarDataArray(:,3) < startPoint(3)+10)),(1:3));
%             [nearPoint,I] = pdist2(otherPoints,startPoint,'euclidean','Smallest',NumNearestPoints);
%             seletedIndices = [seletedIndices I];
%             
%             temotherPoints = [otherPoints(I,:) lidarDataArray(I,6)];
%             sortLiDArray = sortrows(temotherPoints,[4]); %order by 3rd column in decresing order.            
%             clusteredLiDARDataArray(cntd:cntd+(NumNearestPoints-1),:) = otherPoints(I,:);
% 
%             otherPoints(I,:) = 0;%zeros(length(I),3);
%             startPoint = sortLiDArray(NumNearestPoints,(1:3));
%             cntd = cntd + NumNearestPoints;
%         end
%         
%         connectComponentArray{outerLoopCnt} = clusteredLiDARDataArray;
%         connectComponentIndicesArray{clstStartPointIndex} = seletedIndices;
%         outerLoopCnt = outerLoopCnt + 1;
%     end
%     
%     %plotConnectedComponents(connectComponentArray);
%     retConnectedComp = connectComponentArray;
% end

% function voxelSize = getVoxelSize(stData, maxAngle, angleIncrement, removeSimilarRows, neighbourCutoff, showGraph)
% 
%         cntVal =[]; varvalMean =[]; zDiviArr =[];
%         mVari = []; %mMeani=[];
%         for zDiv = 5:5:30 % There is no point starting from low values such as 0
%             
%             % locate the exterior points (new method)
%             [extremeLiDARDataArraySet, extremeLiDARDataIndicesSet] =...
%             locateExteriorPointsEfficient(stData, zDiv, maxAngle, ...
%             angleIncrement, removeSimilarRows, neighbourCutoff);
%             
%             % locate the exterior points (old method)
%             % [extremeLiDARDataArraySet, extremeLiDARDataIndicesSet] = ...
%             % locateExteriorPoints(stData.lidarDataArray, zDiv, ...
%             % maxAngle, angleIncrement);
%         
%             varval =[]; %meanVal =[];
%             
%             for k = 1:1: size(extremeLiDARDataArraySet,1)
%                 pnt =  extremeLiDARDataArraySet(k,1:3);
% 
%                 [D,I] = pdist2(extremeLiDARDataArraySet(:,1:3),pnt,'euclidean','Smallest',10);
%                 Dm = D(2:10);
%                 varval = [varval var(Dm)];
%                 %meanVal = [meanVal mean(Dm)];
%             end
%             mVari = [mVari mean(varval)];
%             %mMeani = [mMeani mean(meanVal)];
%             varvalMean = [varvalMean mean(varval)]; varvalMeani = varvalMean/norm(varvalMean);
%             cntVal = [cntVal size(extremeLiDARDataArraySet,1)]; cntVali = cntVal/norm(cntVal);
%             zDiviArr = [zDiviArr zDiv];
%         end
%         
%         x = linspace(5,30,6);
%         p = polyfit(x,cntVali,4);
%         f1 = polyval(p,x);
%         q = polyfit(x,varvalMeani,2);
%         f2 = polyval(q,x); 
% 
%         if(showGraph)
%             figure(2);
%             hold on;
%             %plot(zDiviArr', cntVali','*', 'Color', 'r'); 
%             plot(x,f1,'r--')
%             
%             %plot(zDiviArr', varvalMeani','*', 'Color', 'b');
%             plot(x,f2,'b--')
%             xlabel('Divisions'); ylabel('Normalized values');            
%             legend('Number of branch points','Varinace of dist between branch point')
%             figure(1);
%         end
%         [X0,~] = intersections(x,f1,x,f2,'ROBUST');
%         
%         voxelSize = X0; % 4 is a buffer increment
% end


% function isValid = isValidExtremePoint(lidarDataArray,startPoint,thresholdvalue)
%     [nearPoint,~] = pdist2(lidarDataArray,startPoint,'euclidean','Smallest',2);          
%     isValid = 1;
%     if(nearPoint(2)>thresholdvalue)
%         isValid = 0;
%     end    
% end

% function [extremeLiDARDataArray, extremeLiDARDataIndex] = locateExteriorPoints(lidarDataArray, zDiv, maxAngle, angleIncrement)
%     count = 1;
%     extremeLiDARDataArray = zeros(zDiv*maxAngle/angleIncrement,size(lidarDataArray,2)); %cell2table(cell(0,5));
%     extremeLiDARDataIndex = zeros(zDiv*maxAngle/angleIncrement,1);
%     %extremeLIDARDataAngle = zeros(zDiv*maxAngle/angleIncrement,1); 
%     maxZ = max(lidarDataArray(:,3)); minZ = min(lidarDataArray(:,3));
%     
%     for zStep = minZ :(maxZ - minZ)/zDiv : maxZ-(maxZ - minZ)/zDiv;
%         for angle = 0 : angleIncrement : maxAngle;
%             
%          lidarDataArray1 = rotateLiDARdata(lidarDataArray(:,(1:3)), angle*(pi/180));
%             clippedLiDARpointsIndices = find(and(lidarDataArray(:,3) >= zStep , lidarDataArray(:,3) < zStep+(maxZ - minZ)/zDiv));
%             
%             %To record right extreme point in a particualr ZDiv
%             [~,tmpIndex] = max(lidarDataArray1(clippedLiDARpointsIndices,1));
%             temparr = lidarDataArray(clippedLiDARpointsIndices,:);
%             extremeLiDARDataArray(count,:) = temparr(tmpIndex,:); % lidarDataArray1(tmpIndex,:);
%             extremeLiDARDataIndex(count,:) =  zStep;
%             %extremeLIDARDataAngle(count,:) = angle; %??
%             count = count + 1;
%             
%             % To record left extreme point in a particualr ZDiv
%             [~,tmpIndex] = min(lidarDataArray1(clippedLiDARpointsIndices,1));
%             temparr = lidarDataArray(clippedLiDARpointsIndices,:);            
%             extremeLiDARDataArray(count,:) = temparr(tmpIndex,:); 
%             extremeLiDARDataIndex(count,:) =  zStep;
%             %extremeLIDARDataAngle(count,:) = 180 + angle; %??
%             
%             count = count + 1;  
%         end
%     end
%     
%     [~, ind] = unique(extremeLiDARDataArray,'rows');    
%     extremeLiDARDataArray = extremeLiDARDataArray(ind,:);
%     extremeLiDARDataIndex = extremeLiDARDataIndex(ind,:);
%     %extremeLIDARDataAngle = extremeLIDARDataAngle(ind,:);
%     
%     isValidArr = zeros(size(extremeLiDARDataArray,1),1);
%     for i = 1:1:size(extremeLiDARDataArray,1)
%         isValidArr(i,1) = isValidExtremePoint(lidarDataArray(:,(1:3)),extremeLiDARDataArray(i,(1:3)), 2); %2 is threshold value
%     end
% 
%     extremeLiDARDataArray = extremeLiDARDataArray(isValidArr==1,1:3);
%     extremeLiDARDataIndex = extremeLiDARDataIndex(isValidArr==1,:)';
%     %extremeLIDARDataAngle = extremeLIDARDataAngle(isValidArr==1,:);
%     
% end

% function circle(x,y,r)
%     %x and y are the coordinates of the center of the circle
%     %r is the radius of the circle
%     %0.01 is the angle step, bigger values will draw the circle faster but
%     %you might notice imperfections (not very smooth)
%     ang=0:0.01:2*pi; 
%     xp=r*cos(ang);
%     yp=r*sin(ang);
%     plot(x+xp,y+yp);
% end


% function [extremeLiDARDataArray, extremeLiDARDataIndex, extremeLIDARDataAngle] = removeRepetitions_oldd(extremeLiDARDaraArray,extremeLiDARDataIndex, extremeLIDARDataAngle,  CutoffDistBetwVect)
%     indx = [];
%     DistBwExtrmPoints = tril(squareform(pdist(extremeLiDARDaraArray(:,1:3),'euclidean')));
%     dims = [size(DistBwExtrmPoints,1),size(DistBwExtrmPoints,2)];   
%     SIndx = find(DistBwExtrmPoints < CutoffDistBetwVect & DistBwExtrmPoints > 0);
% 
%     for i = 1:size(DistBwExtrmPoints,1)*size(DistBwExtrmPoints,2)
%        if ismember(i,SIndx) == 1
%            [row,col] = ind2sub(dims, i);  
%            indx = [indx row];
%        end
%     end
%     indx = unique(indx);
%     extremeLiDARDataArray =removerows(extremeLiDARDaraArray,'ind',indx);
%     extremeLiDARDataIndex= extremeLiDARDataIndex(setdiff((1:1: length(extremeLiDARDaraArray)),indx));
%     extremeLIDARDataAngle = extremeLIDARDataAngle(setdiff((1:1: length(extremeLiDARDaraArray)),indx));
% 
% end
% 
% function rotateLiDARdata = rotateLiDARdata(lidarDataArray, teeta)
%     for i = 1: length(lidarDataArray)
%         lidarPosVector = [lidarDataArray(i,1) lidarDataArray(i,2) lidarDataArray(i,3)];
%         lidarNewPosVector = rotateVector(lidarPosVector,teeta);
%         lidarDataArray(i,1)  = lidarNewPosVector(1);
%         lidarDataArray(i,2)  = lidarNewPosVector(2);
%         lidarDataArray(i,3)  = lidarNewPosVector(3);
%     end
%     rotateLiDARdata = lidarDataArray;
% end

% function getDist = getDist2Point(point1, point2, extPoint)
%   %  s= [point2(1) - point1(1);point2(2) - point1(2);point2(3) - point1(3)];
%   %  M0= extPoint;
%   %  M1 = point1;
%   %  m2m1 = [M1(1)-M0(1) M1(2)-M0(2) M1(3)-M0(3)];   
%   %  getDistu = norm(cross(m2m1,s),2)/norm(s,2);
%     
%     x = extPoint; %some point
%     a = point1; %segment points a,b
%     b = point2;
%     
%     v = b-a;
%     w = x-a;
%     c1 = w*v';
%     c2 = v*v';
%     
%     if(c1<=0)
%        getDist =  sqrt(sum((x - a).^ 2)); 
%        return;
%     elseif (c2<=c1)
%        getDist =  sqrt(sum((x - b).^ 2)); 
%        return;
%     else
%         b1= c1/c2;
%         pb= a + b1*v;
%         getDist =  sqrt(sum((x - pb).^ 2));   
%     end
%     
%   %  d_ab = norm(a-b);
%   %  d_ax = norm(a-x);
%   %  d_bx = norm(b-x);
% 
%   %  if dot(a-b,x-b)*dot(b-a,x-a)>=0
%   %      A = [a,1;b,1;x,1];
%   %      getDist = abs(det(A))/d_ab;        
%   %  else
%   %      getDist = min(d_ax, d_bx);
%   %  end    
% end


% 
% function lidarBranchTips = removeRepetitionsOld(extremeLiDARDaraArray, CutoffDistBetwVect)
%    
%     DistBwExtrmPoints = tril(squareform(pdist(extremeLiDARDaraArray(:,1:3),'euclidean')));
%     %dims = [size(DistBwExtrmPoints,1),size(DistBwExtrmPoints,2)];   
%     %SIndx = find(DistBwExtrmPoints < CutoffDistBetwVect & DistBwExtrmPoints > 0);
%     
%     [I,~] = ind2sub(size(DistBwExtrmPoints),find(DistBwExtrmPoints < CutoffDistBetwVect & DistBwExtrmPoints > 0));
%     indx = unique(I);
%     lidarBranchTips = removerows(extremeLiDARDaraArray,'ind',indx);
% end

%     indicestobeRemoved = [];
%     for i =1:1:size(connectComponentArray,2)        
%         tmpAr = connectComponentArrayIndx{i};
%         ccArr = [];
%         for k = 1:1:size(tmpAr,1)
%            indx = stData.index(tmpAr(k));
%            ccArr(k) = stData.indexComplete(indx);
%         end
%         connectComponentArrayIndxNew{i} = ccArr';
%         
%         %stData == stData.indexComplete
%         indicestobeRemoved = [indicestobeRemoved; connectComponentArrayIndxNew{i}];
%     end
% 
%     indicestobeRemoved = unique(indicestobeRemoved);  
% 
%     
%     if(cIntorempved == indicestobeRemoved)
%        df =0; 
%     end
  
% function featureArr = printResultInConsole(count, retIGFarr, printInConsole)
%     % Print results in the matlab console    
%     if(printInConsole)
%         if(count == 1)
%             disp('  avgRetSlope   avgRetLen   avgPntDistToline   avgMaxEval   avgEigRatio  avgTotalPoints');
%         end
%         disp([retIGFarr.avgRetSlope,retIGFarr.avgRetLen, retIGFarr.avgPntDistToline,...
%             retIGFarr.avgMaxEval,retIGFarr.avgEigRatio,retIGFarr.avgTotalPoints]);
%     end
%     featureArr(count,:) = [retIGFarr.avgRetSlope, retIGFarr.avgRetLen, ...
%         retIGFarr.avgPntDistToline, retIGFarr. avgMaxEval,retIGFarr. ...
%         avgEigRatio,retIGFarr. avgTotalPoints];
% end

%copyfile(fullfile(matlabroot,'extern','examples','mx','mxcreatecharmatrixfromstr.c'),'.','f')

%     subplot(2,3,1);
%     [N1,X1] = hist(retSlopeArr,100);
%     subplot(2,3,2);
%     [N2,X2] = hist(retLenArr/TWidth,100);
%     subplot(2,3,3);
%     [N3,X3] = hist(pntDistTolineArr,100);
%     subplot(2,3,4);
%     [N4,X4] = hist(maxevalArr/TreeHeight,100);
%     subplot(2,3,5);
%     [N5,X5] = hist(eigRatioArr,100);
%     subplot(2,3,6);
%     [N6,X6] = hist(totalPointsArr/TotalNumberOfPoints,100);
%     
%     [~,ai] = max(retSlopeArr); 
%     retIGFarr.avgRetSlope1 = mean(X1(ai-1: ai+1));
%     [~,bi] = max(retLenArr/TWidth); retIGFarr.avgRetLen1 =  mean(X2(bi-1: bi+1));
%     [~,ci] = max(pntDistTolineArr); retIGFarr.avgPntDistToline1 =  mean(X3(ci-1: ci+1));
%     [~,di] = max(maxevalArr/TreeHeight); retIGFarr.avgMaxEval1 =  mean(X4(di-1: di+1));
%     [~,ei] = max(eigRatioArr); retIGFarr.avgMaxEval1 =  mean(X5(ei-1: ei+1));
%     [~,fi] = max(totalPointsArr/TotalNumberOfPoints); retIGFarr.avgTotalPoints1 = mean(X6(fi-1: fi+1));
