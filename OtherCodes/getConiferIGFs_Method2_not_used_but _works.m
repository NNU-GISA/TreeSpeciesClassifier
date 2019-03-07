function igffeatureArr = getConiferIGFs_Method2(inFileFullName, speciesCode, allFiles, isPlotOnAll)
    % Set paths
    % inFilepath = 'C:\My_Files\1_PhD_Research\1_PhD_Research_Topic\1_LiDAR_Forest_Applications\1_Conifer_Species_Classification\1_Matlab_Programs\CropPruneLiDARData\2_PruneData\PruneOutput\ar\';
    %inFilepath = 'D:\CSC\PruneOutput\ar\';
    exePath='C:\My_Files\3_Software_Packages\Lidar_Tools\LiDAR_CPP_Tools\LasTools\bin\';
    %exePath = 'D:\Studies\LiDAR_Research\LASTools\bin\';
    outFilepath = 'C:\My_Files\1_PhD_Research\1_PhD_Research_Topic\1_LiDAR_Forest_Applications\3_Research_Files\1_Conifer_Species_Classification\1_Matlab_Files\1_ConiferSpeciesDetection\';
        
    %randn('state',23432); rand('state',3454); 
       
    %files = dir(strcat(inFilepath,'*.las'));
    igffeatureArr = zeros(size(allFiles,1),6);
    count = 1;
    %isPlotOnAll = false;
    
    for file = allFiles'
        tic
        if(isPlotOnAll)
            try
                clf(Figure1,'RESET'); 
                clf(Figure2,'RESET');
            catch
                clf;
            end
        end
        sdf =file.name;
        % Read LiDAR data
        fullFileName = strcat(inFileFullName,sdf); % get full file name
        singleTreeLiDARdata = LoadData(fullFileName);

        % LiDAR data Pre-processing and basic tree information retrieval
        stData = dataPreProcessing(singleTreeLiDARdata);
        
        % plot LiDAR data (original space)      
        pltLiDARTree = true;
        plotLiDARData(stData, pltLiDARTree, isPlotOnAll);

        % Plot projected/slabbed LiDAR data
        plotProjectedData = true; 
        slabbedCoordinates = getProjectedLiDARdata(stData, plotProjectedData, isPlotOnAll);
        %slabbedCoordinatesold = getProjectedLiDARdataold(stData, plotProjectedData);

                
        if(strcmp(speciesCode,'ar'))
            proximThreshold = 1;
        elseif(strcmp(speciesCode,'la'))
            proximThreshold = 1.6;
        elseif(strcmp(speciesCode,'ab'))
            proximThreshold = 1.2;
        else
            proximThreshold = 0.3;
        end
        
        % Detect and Plot branch tips from DEM (in projected space)         
        removeProximalTips = true; plotBranchTips = false; % proximThreshold = 1; 
        branchTips = getBranchTips(slabbedCoordinates(:,1:3), removeProximalTips, proximThreshold, outFilepath, exePath, plotBranchTips, isPlotOnAll);
        retRand = distinguishable_colors(size(branchTips.branchTipPoints,1));
         
        % IDENTIFY LIDAR-BRANCH POINT CLUSTERS 
        pltKMeanResult = false; pltParaboloid = false; pltBranchCluster = true;
        [cIndx, clusteredPointArray, peakPointArr] = getBranchClusters(stData, slabbedCoordinates, branchTips, pltKMeanResult, pltParaboloid, pltBranchCluster, retRand, isPlotOnAll);

        
        % Get stem point in the direction of PC1 of branch clusters.
        [stemPointArr, ~] = getStemPointFromPCA(peakPointArr, clusteredPointArray, stData.retMaxXYZ);

        % GET INTERNAL GEOMETRIC FEATURES
        slopeThreshold = 60; plotBranchLines=true;
        igffeatureArr(count,:) = getBranchIGFs(stData, cIndx, stemPointArr, peakPointArr, slopeThreshold, plotBranchLines, isPlotOnAll);
        count = count + 1;
        toc
    end
    csvwrite('csvlist_igfs_all_trees.csv',igffeatureArr);
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

function retIGFsVec = getBranchIGFs(stData, cIndx, compStartPoint, branchTips, slopeThreshold, plotBranchLines, isPlotOnAll)
    numBranches = size(cIndx,2);
    retIGFsArray = zeros(numBranches,6);
    
    if(and(plotBranchLines,isPlotOnAll))
        figure(2);
        subplot(1,2,2);
        hold on;
    end

    for i = 1:1:numBranches
        newIndxList = cIndx{i};
        if(size(newIndxList,1)>3)

            branchCluster = stData.lidarDataArray(newIndxList,1:3);
            compStartP = compStartPoint(i,:);
            compEndP = branchTips(i,1:3);

            % Get best fit line for the modified connected component        stData, dataArray, compStartPoint, compEndPoint, slopeThreshold, isPlotON, plotOn, colorArr
            [retSlope, lineLength, pntDistToline, maxeval, eigRatio, totalPoints] = getBestFitLine(stData, branchCluster, compStartP, compEndP, slopeThreshold, plotBranchLines, isPlotOnAll); % isPlotON
            %[lineSlopeArr(i,1), lineLengthArr(i,1), avgDistToLineArr(i,1)] = getBestFitLine(branchCluster, ignoreSlope, plotBranchLines, isPlotOnAll, [0.5 0 0]); % isPlotON

            if(eigRatio>100) % to avoid large values due to division by small number. 
                eigRatio = -1; %(eigratio)
            end
            
            retIGFsArray(i,:) = [retSlope, lineLength, pntDistToline, maxeval,eigRatio,totalPoints]; 
        end
    end   

    newEigRatio = retIGFsArray(and(retIGFsArray(:,5)>0,retIGFsArray(:,5)<10),5); 
    newEigRatio = sum(newEigRatio)/size(newEigRatio,1); 
    
    % Select branches which meets the above slopeThreshold criteria
    BranchArr  = retIGFsArray(abs(retIGFsArray(:,1)) < slopeThreshold,:);    
    retIGFsVec = sum(BranchArr,1)/size(BranchArr,1);

    % Perform height/desnity normalization for three of the features
    retIGFsVec(2) = retIGFsVec(2)/stData.treeHeight; % Branch Length
    retIGFsVec(4) = retIGFsVec(4)/stData.treeHeight; % Max  eigenvalue
    retIGFsVec(5) = newEigRatio;
    retIGFsVec(6) = retIGFsVec(6)/size(stData.lidarDataArray,1); % Total Points in a branch

    %branchIGFs = [AvgSlope AvglineLength distAvg maxeigArr eifArr avgNumberPoints]; % same order as PCA based one
    
    if(and(plotBranchLines, isPlotOnAll))
        figure(2);
        subplot(1,2,2);
        hold off;
    end
end

function [cIndx ,clusteredPointArray, peakPointArr] = getBranchClusters(stData, slabbedCoordinates, branchTips, plotKMeanResult, pltParaboloid, pltBranchCluster, retRand, isPlotOnAll)
    cIndx={};
    
    if(and(or(plotKMeanResult, or(pltParaboloid, pltBranchCluster)),isPlotOnAll))
        figure(2);
        subplot(1,2,1);
        hold on;
        xlabel('Tree Height','Fontname', 'Times New Roman' ,'FontSize', 16); 
        ylabel('Distance to Reference Point','Fontname', 'Times New Roman' ,'FontSize', 16); 
        zlabel('Crown Radius','Fontname', 'Times New Roman' ,'FontSize', 16);
        set(gca,'XLim',[0 stData.treeHeight]);set(gca,'YLim',[-20 20]);set(gca,'ZLim',[0 stData.maxTreeWidth]);
        set(gca,'XTick',0:2:stData.treeHeight); set(gca,'YTick',-20:2:20); set(gca,'ZTick',0:1:stData.maxTreeWidth)
    end
    
    %retClustIndxArray = getConnectedComponets(slabbedCoordinates, branchTips.branchTipPoints);
    [clusteredPointArray, retClustIndxArray, peakPointIndexArr, peakPointArr] = do_Kmeans(stData, slabbedCoordinates, size(branchTips.branchTipPoints,1), branchTips.branchTipPoints, plotKMeanResult);
    
    for j1 = 1:1:size(retClustIndxArray,2)
        branchPointIndex = retClustIndxArray{j1};            
        branchTipsIndx = peakPointIndexArr(j1);
        
        if(size(branchPointIndex,1)>1)

                inputBClustArr = slabbedCoordinates(branchPointIndex,:);
                   
                xShift =  mean(inputBClustArr(:,1)); 
                yShift =  mean(inputBClustArr(:,2)); 
                maxZVal = max(inputBClustArr(:,3));

                minVal = max(inputBClustArr(:,3));
                inputBClustArr(:,3) = (inputBClustArr(:,3)-minVal);

                % GET PARABOLOID PARAMETERS
                isParabloid = true; % change if another geometric shape is used
                coeffArr(:,j1) = fitParabloid(inputBClustArr);
                plotShape(coeffArr(:,j1), xShift, yShift, maxZVal, isParabloid, pltParaboloid);
                
                inputBClustArr(:,1) = inputBClustArr(:,1)-mean(inputBClustArr(:,1));
                inputBClustArr(:,2) = inputBClustArr(:,2)-mean(inputBClustArr(:,2));

                slabbedCoordinates1 = slabbedCoordinates;
                slabbedCoordinates1(:,3) = (slabbedCoordinates(:,3)-minVal);

                samePointinOrigPntCld = slabbedCoordinates(find(slabbedCoordinates(:,4)==branchTipsIndx),1:3);
                samePointintmpPntCld =  inputBClustArr(find(inputBClustArr(:,4)==branchTipsIndx),1:3);               
                shiftinSamePoint = samePointinOrigPntCld - samePointintmpPntCld;
                slabbedCoordinates1(:,1) = slabbedCoordinates(:,1) - shiftinSamePoint(1);
                slabbedCoordinates1(:,2) = slabbedCoordinates(:,2) - shiftinSamePoint(2);

                % FIND INDICES OF POINTS INCLUDED WITHIN THE ELLIPSOID AND STORE IT IN THE "cIndx" STRUCTURE
                aa = [];
                slabbedCoordinates1 = slabbedCoordinates1(and(abs(slabbedCoordinates1(:,1))<2.5, abs(slabbedCoordinates1(:,2))<2.5),:);
                for pt =1:1:size(slabbedCoordinates1,1);                    
                    isIncluded = isPointIncluded(slabbedCoordinates1(pt,1:3),coeffArr(1,j1),coeffArr(2,j1),coeffArr(3,j1),1,'elparaboloid');
                    if(~isIncluded)
                       % plot3(slabbedCoordinates1(pt,1), slabbedCoordinates1(pt,2), slabbedCoordinates1(pt,3),'.', 'MarkerSize' ,40, 'Color',[0 0 0.5]);
                       aa = [aa slabbedCoordinates1(pt,4)];
                    end
                end                    
                cIndx{j1} = find(ismember(slabbedCoordinates(:,4),aa));
        end
    end

    % FIND CLUSTER CENTERS AND INDICES OF POINTS THAT WERE ALREADY
    % ASSIGNED TO ONE OF THE CLUSTERS.
    AllConvexPointIndexSet =[];
    newClusCenter =[];
     for kk= 1:1:size(cIndx,2)           
        pointSetIndexCurrent = cIndx{kk};
        AllConvexPointIndexSet = [AllConvexPointIndexSet; pointSetIndexCurrent];            
        XX = slabbedCoordinates(pointSetIndexCurrent,1);
        YY = slabbedCoordinates(pointSetIndexCurrent,2);            
        newClusCenter = [newClusCenter; mean(XX) mean(YY)];
     end        
    %FIND INDICES OF THE REMAINING UNASSIGNED POINTS
    remainingPoints = setdiff([1:1:size(slabbedCoordinates,1)],AllConvexPointIndexSet);

    %ASSIGN PREVIOUSLY UNALLOCATED POINTS TO NEAREST CLUSTER (BASED ON THE 
    %DISTANCE TO ITS CLUSTER CENTER)        
    for kk= 1:1:size(remainingPoints,2)           
        pointCurrent = slabbedCoordinates(kk,1:2);               
        minDist = pdist2(pointCurrent,newClusCenter);
        [minminDist, minminIndx] = min(minDist);
        if(minminDist <0.5)
            newIndxList = cIndx{minminIndx};  
            newIndxList = [newIndxList; slabbedCoordinates(kk,4)];
            cIndx{minminIndx} = newIndxList;
        end         
    end
    
    if(and(pltBranchCluster,isPlotOnAll))
        clf;
        for bcIdx = 1:1:size(cIndx,2) 
            indxList = cIndx{bcIdx};
            subplot(1,2,1);
            hold on;
            plot3(slabbedCoordinates(indxList,1), slabbedCoordinates(indxList,2), slabbedCoordinates(indxList,3),'.','MarkerSize',15,'Color',retRand(bcIdx,:));
            subplot(1,2,2);
            hold on;
            plot3(stData.lidarDataArray(indxList,1), stData.lidarDataArray(indxList,2), stData.lidarDataArray(indxList,3),'.','Color',retRand(bcIdx,:));
            %[K,V] = boundary(stData.lidarDataArray(newIndxList,1),stData.lidarDataArray(newIndxList,2),stData.lidarDataArray(newIndxList,3),0.5);
            %trisurf(K,stData.lidarDataArray(newIndxList,1),stData.lidarDataArray(newIndxList,2),stData.lidarDataArray(newIndxList,3));
            
            plot3(stData.lidarDataArray(peakPointIndexArr(bcIdx),1), stData.lidarDataArray(peakPointIndexArr(bcIdx),2), stData.lidarDataArray(peakPointIndexArr(bcIdx),3),'.','MarkerSize',20,'Color',[1 0 0]);   
            
        end
    end
    
    if(and(or(plotKMeanResult, or(pltParaboloid, pltBranchCluster)),isPlotOnAll))
        subplot(1,2,1);
        camproj perspective; rotate3d on; axis vis3d; axis square; axis on; grid on; view(60, 45);
        xlabel('Tree Height','Fontname', 'Times New Roman' ,'FontSize', 16); 
        ylabel('Distance to Reference Point','Fontname', 'Times New Roman' ,'FontSize', 16); 
        zlabel('Crown Radius','Fontname', 'Times New Roman' ,'FontSize', 16);
        set(gca,'XLim',[0 stData.treeHeight]);set(gca,'YLim',[-20 20]);set(gca,'ZLim',[0 stData.maxTreeWidth]);
        set(gca,'XTick',0:2:stData.treeHeight); set(gca,'YTick',-20:2:20); set(gca,'ZTick',0:1:stData.maxTreeWidth)
        axis([0 stData.treeHeight -20 20 0 stData.maxTreeWidth]);
        subplot(1,2,2);
        axis([-stData.maxTreeWidth stData.maxTreeWidth -stData.maxTreeWidth stData.maxTreeWidth stData.mintreeHeight stData.maxtreeHeight]);
        camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-45, 15)
        xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 16); 
        ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 16); 
        zlabel('Tree Height','Fontname', 'Times New Roman' ,'FontSize', 16);
        set(gca,'XLim',[-stData.maxTreeWidth stData.maxTreeWidth]);set(gca,'YLim',[-stData.maxTreeWidth stData.maxTreeWidth]);set(gca,'ZLim',[0 stData.treeHeight]);
        set(gca,'XTick',-stData.maxTreeWidth:2:stData.maxTreeWidth); set(gca,'YTick',-stData.maxTreeWidth:2:stData.maxTreeWidth); set(gca,'ZTick',0:2:stData.treeHeight)
    end
end
% stData, dataArray, compStartPoint, compEndPoint, slopeThreshold, isPlotON, plotOn, colorArr
function [retSlope, lineLength, pntDistToline, maxeval, eigRatio, totalPoints] = getBestFitLine(stData, dataArray, compStartPoint, compEndPoint, slopeThreshold, isPlotON, plotAllOn)
    
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
    % endpnt = compEndPoint; % if branch tips are sure to be end points

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

    % Plot line
    if(abs(retSlope) < slopeThreshold) % && compStartPoint(3) > 5
        if(and(isPlotON,plotAllOn))
            %hold on
           % plot3(compStartPoint(1),compStartPoint(2),compStartPoint(3),'.','Color',[0 1 0],'MarkerSize',40);
           % plot3(endpnt(1),endpnt(2),endpnt(3),'.','Color',[1 0 0],'MarkerSize',40);
            %plot3(m(1),m(2),m(3),'.','Color',[0 0 1],'MarkerSize',40);
            lnData = [compStartPoint; endpnt];
            p1 = line(lnData(:,1),lnData(:,2),lnData(:,3), 'Color',[139/256 69/256 19/256],'LineWidth', 3); alpha(0.9);
            p1.Color(4) = 0.7;
            %plot3(x1,y1,z1,'LineWidth',3, 'Color',colorArr);
            %hold off;
        end
    end

    lineLength = (dst1 + dst2);

    % Calculate Eigen features
    retEigenData = getBestEigenDirection(dataArray);
    %evec = retEigenData{1};
    eval = retEigenData{2};

    % handle cases where only one/two PC is available.
    maxeval = 0; eigRatio = 1;
    if(length(eval)==3)
        maxeval = max(eval(2:end));
        eigRatio = eval(2)/eval(3);
    elseif(length(eval)==2)
        maxeval = max(eval(2:end));
        eigRatio = 1;
    else            
        maxeval = 0; eigRatio =0.5;
    end

    % total no of LiDAR point in the branch cluster
    totalPoints = size(dataArray,1);
        
end

% function [retSlopeArr, retLengthArr, avgDistToLine] = getBestFitLine(dataArray, slopeThreshold, isPlotON, isPlotOnAll, colorArr)    
% 
%         % Fit a line throught the cluster of points (to get branch skeleton)
%         [m,p,avgDistToLine] = best_fit_line(dataArray(:,1),dataArray(:,2), dataArray(:,3));
% 
%         dist_matrix  = squareform(pdist(dataArray, 'euclidean'));
%         [maxNum, maxIndex] = max(dist_matrix(:));
% 
%         % For considering that points are located at different directions around the stem.        
%         t=-maxNum/2:0.2:maxNum/2;
%         x1=m(1)+p(1)*t;
%         y1=m(2)+p(2)*t;
%         z1=m(3)+p(3)*t;
% 
%         % Get slope of the line with respect to the X-Y plane
%         %retSlopeArr = getLineSlope(x1(1), x1(2), y1(1), y1(2), z1(1), z1(2));
%         retSlopeArr = subspace(p',[p(1:2) 0]')*(180/pi);
%         
%         retLengthArr =  maxNum;
%         
%         % Plot line
%         if(and(isPlotON,isPlotOnAll))
%             if(abs(retSlopeArr) < slopeThreshold)
%                 plot3(x1,y1,z1,'LineWidth',3, 'Color',colorArr);
%             end
%            % hold off;
%         end
% end

function plotShape(coeffArr, xmid, ymid, maxZVal, isParabloid, pltParaboloid)
    a = coeffArr(1); b = coeffArr(2); c = coeffArr(3);
    
    X1= [-2.5:0.1:2.5];
    Y1= [-2.5:0.1:2.5];
    [X,Y] = meshgrid(X1,Y1);
    
    Z = [];
    if(isParabloid)
        % Plot Parabloid
        Z = myparabloid(X,Y,a,b,c);  
    else
        % Plot myCone
        Z = myCone(X,Y,a,b,c);  
    end 
    if(pltParaboloid)
        subplot(1,2,1);
        surf(X+xmid,Y+ymid,Z+maxZVal);
    end
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


% function [retSlope, lineLength, pntDistToline, maxeval, eigRatio, totalPoints] = getBestFitLine(stData, dataArray, compStartPoint, compEndPoint, slopeThreshold, isPlotON, colorArr)
%     
%         % Fit a line throught the cluster of points (to get branch skeleton)
%         [m,p,pntDistToline] = best_fit_line(dataArray(:,1),dataArray(:,2), dataArray(:,3));
%         % Calcuate the slope of the line (in Degrees)
%         retSlope = subspace(p',[p(1:2) 0]')*(180/pi);
%         
%         % find the nearest point on the stem where the branch lines crosses
%         % as as line start point
%         stemBottom =[stData.retMaxXYZ(:,1:2) 0];
%         stemTop =[stData.retMaxXYZ(:,1:2) 100];        
%         P3 = [m(1)+p(1)*100 m(2)+p(2)*100 m(3)+p(3)*100];
%         P4 = [m(1)+p(1)*-100 m(2)+p(2)*-100 m(3)+p(3)*-100];        
%         [distToStem, nearestStemPoints] = distBtwLineSegments(stemBottom, stemTop, P3, P4);        
%         %pt = nearestStemPoints{2}; 
%         
%         % find the extreme point in cluster as line end point
%         [~, Dindx] = max(sqrt(sum((repmat(compStartPoint,size(dataArray,1),1) - dataArray).^ 2,2)));  
%         endpnt = dataArray(Dindx,:); 
%          %endpnt = compEndPoint; % if branch tips are sure to be end points
%         
%         %hold on;
%         %plot3(pt(1),pt(2),pt(3),'*','Color','r', 'MarkerSize', 18);
%         
%         % distance from cluster centre to branch start point(i.e. near trunk)
%         dst1 = pdist2(m(1,:),compStartPoint,'euclidean','Smallest',1);    
%         % distnace of Cluster centre from the branch exterior points.
%         dst2 = pdist2(endpnt,m(1,:),'euclidean','Smallest',1);
%         
%         
%         % Considering that points are located at different directions around the stem.
%         if(endpnt(1,1)>0)
%             t=-dst2:0.2:dst1;
%         else
%             t= -dst1:0.2:dst2;
%         end
%         
%         x1=m(1)+p(1)*t; y1=m(2)+p(2)*t; z1=m(3)+p(3)*t;
%         
%         % jugaad fix (need to correct)
%         DI = max(sqrt(sum((repmat(compStartPoint,size(x1,2),1) - [x1' y1' z1']).^ 2,2)));  
%         DIend = max(sqrt(sum((compStartPoint - endpnt).^ 2,2)));        
%         if(DI>DIend+0.3)
%             t =-t; x1=m(1)+p(1)*t; y1=m(2)+p(2)*t; z1=m(3)+p(3)*t;            
%         end
%         
%         % Plot line
%         if(abs(retSlope) < 45 && compStartPoint(3) > 5) 
%             if(strcmp(isPlotON,'true'))
%                 %hold on
%                % plot3(compStartPoint(1),compStartPoint(2),compStartPoint(3),'.','Color',[0 1 0],'MarkerSize',40);
%                % plot3(endpnt(1),endpnt(2),endpnt(3),'.','Color',[1 0 0],'MarkerSize',40);
%                 %plot3(m(1),m(2),m(3),'.','Color',[0 0 1],'MarkerSize',40);
%                 lnData = [compStartPoint; endpnt];
%                 p1 = line(lnData(:,1),lnData(:,2),lnData(:,3), 'Color',[139/256 69/256 19/256],'LineWidth', 3); alpha(0.9);
%                 p1.Color(4) = 0.7;
%                 %plot3(x1,y1,z1,'LineWidth',3, 'Color',colorArr);
%                 %hold off;
%             end
%         end
%         
%         lineLength = (dst1 + dst2);
%         
%         % Calculate Eigen features
%         retEigenData = getBestEigenDirection(dataArray);
%         evec = retEigenData{1};
%         eval = retEigenData{2};
%         
%         % handle cases where only one/two PC is available.
%         maxeval = 0; eigRatio = 1;
%         if(length(eval)==3)
%             maxeval = max(eval(2:end));
%             eigRatio = eval(2)/eval(3);
%         elseif(length(eval)==2)
%             maxeval = max(eval(2:end));
%             eigRatio = 1;
%         else            
%             maxeval = 0; eigRatio =1;
%         end
%         
%         % total no of LiDAR point in the branch cluster
%         totalPoints = size(dataArray,1);
%         
% end

function isIncluded = isPointIncluded(point,a,b,c,r,shapeType)   
    isIncluded = 0;%a=1.4;b=1.4;c=12.4;
    x = point(1); y = point(2); z = point(3);
    Z =[];
    if (strcmp(shapeType, 'elcone'));
        Z = myCone(x,y,a,b,c);
    elseif (strcmp(shapeType, 'elparaboloid'));
        Z = myparabloid(x,y,a,b,c);
    elseif(strcmp(shapeType,'elCylinder'));        
        dis = sqrt(sum((x - y).^ 2));
        if(dis > r)
            Z = z+1;
        else
            Z = z-1;
        end
    end
    if(z>Z)
        isIncluded = 1;
    elseif(z==Z)
        isIncluded = 2;
    end
end

function coeffVal  = fitParabloid(inputPoints)
    x =inputPoints(:,1); y = inputPoints(:,2); z = inputPoints(:,3);

    x = x- mean(x);
    y = y- mean(y);

    A = [x.*x y.*y];
    
    isParabloid = true;    
    if(isParabloid)
        coeff =  A\z; % for hyperboloid
        a = coeff(1); b = coeff(2);
        c = mean(z(:)./ ( (power(x(:),2)*a^2) + (power(y(:),2)*b^2) ) ); %old 
        
%         if(a>b)
%             temp=b;
%             a=b;
%             a=temp;
%         end
%         if(c<-3)
%             c = -3;
%         end
%         %a = 21; b = 10; c = 0.05;
%         %a = ev(1); b = ev(2); c = 60;
%         %c=0.2;
% 
%         if(abs(a/b) > 2)
%         a= 25; b=10; c =-0.05;
%         elseif(abs(b/a) > 2)
%         a= 25; b=10; c =-0.05;
%         end
         a= 25; b=10; c =-0.04;
        
    else 
        coeff =  A\(z.*z);  % for cone
        a = coeff(1); b = coeff(2);
        c = mean(sqrt( power(z(:),2)./ ( (power(x(:),2)*a^2) + (power(y(:),2)*b^2) ) )); % cone
        %a = 25; b = 10; c = 0.2;
    end

    coeffVal = [a b c];
end

function z = myparabloid(x,y,a,b,c)
   z = ( ( (power(x,2)*a^2) + (power(y,2)*b^2) ) * c);
  %Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
end

function z = myCone(x,y,a,b,c)
  z = sqrt( ( (power(x,2)*a^2) + (power(y,2)*b^2) ) * c^2);
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
    lidarDataDensityArr = ones(size(lidarDataArr,1),1);  %calculateDistWeights(lidarDataArr,0.25);
    lidarDataArr = [lidarDataArr index' lidarDataDensityArr];
    dataAtt.lidarDataArrayComplete =  lidarDataArr(:,1:8);
    
    % Sparsified LiDAR data Array
    lidarDataArrSampled = lidarDataArr(lidarDataDensityArr(:,1) >= 0.0010,:);    
    dataAtt.lidarDataArray =  lidarDataArrSampled(:,1:8);                          
end

% function dataAtt = dataPreProcessing(singleTreeLiDARdata)
% % Write normalized data to a table for performance improvement 
%         lidarDataArr = normalizeLiDARData(write2table(singleTreeLiDARdata));
%         %lidarDataArray = lidarDataArray(randperm(size(lidarDataArray,1),9000),:);
%                
%         treeWidth =  max(lidarDataArr(:,1)); treeHeight = max(lidarDataArr(:,2)); % to calculate max wisth/breadth to set plot width
%         dataAtt.maxTreeWidth = max(treeWidth, treeHeight) + 0.5; % Keep it as 5 if issues arise
%         dataAtt.mintreeHeight = min(lidarDataArr(:,3));
%         dataAtt.maxtreeHeight = max(lidarDataArr(:,3));        
%         dataAtt.treeHeight = dataAtt.maxtreeHeight - dataAtt.mintreeHeight;        
%         dataAtt.htDeduction = dataAtt.treeHeight*0.05; % to get rid of ground noise points
%         
%         % Crown height
%         lidarDataArr = lidarDataArr(find(and(lidarDataArr(:,3) > min(lidarDataArr(:,3))+ dataAtt.htDeduction, lidarDataArr(:,3) < max(lidarDataArr(:,3)))),:);
%         dataAtt.minCrownHeight = getMinCrownHeight(lidarDataArr);
%         dataAtt.maxCrownHeight = max(lidarDataArr(:,3));        
%         dataAtt.crownHeight = dataAtt.maxCrownHeight - dataAtt.minCrownHeight;       
%               
%         %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 > 0.1),:);                
%         %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 < 12),:);
%         
%         % Identify the center point of the tree (in top view).
%         dataAtt.retMaxXYZ = findMaxHeightXY(lidarDataArr);
%         
%         lidarDataDensityArr =  calculateDistWeights(lidarDataArr,0.25); 
%         lidarDataArr = lidarDataArr(find(lidarDataDensityArr(:,1) > 0.0025),:); 
%         
%         lidarDataDensityArr =  calculateDistWeights(lidarDataArr,0.25); %get point in very close neightbourhood
%                      
%         index = 1:1:size(lidarDataArr,1); 
%         
%         
%         % % Density priuned LiDAR data
%         selIndx = find(lidarDataDensityArr(:,1) > 0.0025);        
%         dataAtt.lidarDataArray =  lidarDataArr(selIndx,:);                
%         dataAtt.index = index(selIndx);  % 1:1:size(dataAtt.subSampledlidarDataArray,1); 
%         dataAtt.lidarDataDensityArr = lidarDataDensityArr(selIndx,:);        
%         dataAtt.lidarDataArray = [dataAtt.lidarDataArray dataAtt.index' dataAtt.lidarDataDensityArr];
%         
%         % Complete LiDAR data
%         dataAtt.lidarDataArrayComplete =  lidarDataArr;              
%         dataAtt.indexComplete = index;
%         dataAtt.lidarDataDensityArrComplete = lidarDataDensityArr;       
%         dataAtt.lidarDataArrayComplete = [dataAtt.lidarDataArrayComplete dataAtt.indexComplete' dataAtt.lidarDataDensityArrComplete];
%         
%         
%         %dataAtt.lidarDataArray = lidarDataArray(randperm(size(lidarDataArray,1),9000),:); % for reducing the no. of samples for testing only.
%         
%       %  dataAtt.lidarDataDensityArr =  calculateDistWeights(lidarDataArr,0.25); % get point in very close neightbourhood
%   
% end

function branchTips = getBranchTips(slabbedCoordinates, removeProximalTips, proximThreshold, outFilepath, exePath, plotBranchTips, isPlotOnAll)
        % Make DEM from LAS
        makeTXTFromLASMAIN(slabbedCoordinates,outFilepath);      
        isSuccess =  makeLASFromTXT(exePath,'txt2las',outFilepath,'TempFile.txt',outFilepath,'Out.las');
        if(isSuccess == 0)
            %disp(strcat('Sucessfully created .las file'));
            fclose('all');
            delete(strcat(outFilepath,'TempFile.txt'));
           % h = msgbox('Operation Completed Successfully','Success');
        end      
        isSuccess = makeDEMFromLAS(exePath,'las2dem',outFilepath,'Out.las',outFilepath,'ABC.tif');
        
        % Peak detection in the projected space
        peakPoints = Peak_detection(slabbedCoordinates, strcat(outFilepath,'ABC.tif'));
        %znew = ((peakPoints(:,3)/max(peakPoints(:,3)))*(4.00-0.31)) + 0.31;
        pPoints = [peakPoints(:,1) peakPoints(:,2) peakPoints(:,3)];
        indxt = [];
        for lnr = 1:1:size(pPoints,1)        
            [~,I] = pdist2(slabbedCoordinates,pPoints(lnr,:),'euclidean','Smallest',1);        
            indxt = [indxt I];
        end        
        % Branch Tips Points Indices
        branchTips.indxt = indxt';
        branchTips.branchTipPoints = slabbedCoordinates(branchTips.indxt,:);
        if(removeProximalTips)
            % Remove any identified exterior branch end points which are very close to each other.
            [branchTips.branchTipPoints, branchTips.indxt] = removeRepetitions([branchTips.branchTipPoints branchTips.indxt],  proximThreshold);
        end
        
        if(and(plotBranchTips,isPlotOnAll))
           subplot(1,2,2);
           hold on;
           plot3(slabbedCoordinates(branchTips.indxt,1), slabbedCoordinates(branchTips.indxt,2), slabbedCoordinates(branchTips.indxt,3),'.','MarkerSize',30,'Color',[1 0 0]);   
           hold off;
        end
end

function projLiDARdata = getProjectedLiDARdataold(stData, plotProjectedData)
    projLiDARdata = zeros(size(stData.lidarDataArray,1),7);
    % Find extreme branch point
    retVal = extremeBranchPoint(stData.lidarDataArray, stData.retMaxXYZ);
    for i = 1:1:size(stData.lidarDataArray,1)            
        inputPoint = stData.lidarDataArray(i,1:3);            
        Ax = retVal(1); Ay = retVal(2); Bx = -retVal(1); By = -retVal(2); X = inputPoint(1); Y = inputPoint(2);        
        position = sign((Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax)); % +1 on one side, -1 on the other side.
        [sectorLength, disttocentre] = getSectorLengthold(retVal, inputPoint, [0 0]); %% disttoCentre = Y, sectorlength = X,  inputPointZ = Z  
        projLiDARdata(i,1:7) = [inputPoint(3) position*sectorLength  disttocentre 0 inputPoint];    
    end
    if(plotProjectedData)
        %camproj perspective; rotate3d on; axis vis3d; axis square; axis on; view(-45, 0); grid on;% axis equal; 
        axis([-stData.maxTreeWidth stData.maxTreeWidth stData.mintreeHeight stData.maxtreeHeight]);
        plot3(projLiDARdata(:,1), projLiDARdata(:,2), projLiDARdata(:,3),'.','Color',[0 0.5 0]);
        xlabel('Y Axis'); ylabel('Y Axis'); zlabel('Tree Height');
    end
    projLiDARdata(:,4) = [1:1:size(stData.lidarDataArray,1)]';
end

function projLiDARdata = getProjectedLiDARdata(stData, plotProjectedData, isPlotOnAll)
    % Find extreme branch point
    retVal = extremeBranchPoint(stData.lidarDataArray, stData.retMaxXYZ);
         
    inputPoint = stData.lidarDataArray(:,1:3);            
    Ax = retVal(1); Ay = retVal(2); Bx = -retVal(1); By = -retVal(2); X = inputPoint(:,1); Y = inputPoint(:,2);        
    position = sign((Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax)); % +1 on one side, -1 on the other side.
    [sectorLength, disttocentre] = getSectorLength(retVal, inputPoint, [0 0]); %% disttoCentre = Y, sectorlength = X,  inputPointZ = Z  
    projLiDARdata = [inputPoint(:,3) position.*sectorLength  disttocentre zeros(size(position,1),1) inputPoint];    

    if(and(plotProjectedData,isPlotOnAll))
        subplot(1,2,2); 
        hold on;
        camproj perspective; rotate3d on; axis vis3d; axis square; axis on; view(60, 45); grid on;% axis equal; 
        axis([-stData.maxTreeWidth stData.maxTreeWidth stData.mintreeHeight stData.maxtreeHeight]);
        plot3(projLiDARdata(:,1), projLiDARdata(:,2), projLiDARdata(:,3),'.','Color',[0 0.5 0]);
        xlabel('Tree Height','Fontname', 'Times New Roman' ,'FontSize', 16); 
        ylabel('Distance to Reference Point','Fontname', 'Times New Roman' ,'FontSize', 16); 
        zlabel('Crown Radius','Fontname', 'Times New Roman' ,'FontSize', 16);
        set(gca,'XLim',[0 stData.treeHeight]);set(gca,'YLim',[-20 20]);set(gca,'ZLim',[0 stData.maxTreeWidth]);
        set(gca,'XTick',0:2:stData.treeHeight); set(gca,'YTick',-20:2:20); set(gca,'ZTick',0:1:stData.maxTreeWidth)
        axis([0 stData.treeHeight -20 20 0 stData.maxTreeWidth]);
        hold off;
    end
    projLiDARdata(:,4) = [1:1:size(stData.lidarDataArray,1)]';
end


% function plotBranchClusters(stData, retClustIndxArray)  
%     for i=1:1:size(retClustIndxArray,2) 
%        hold on;
%        tempIndx = retClustIndxArray{i};
%        plot3(stData.lidarDataArray(tempIndx,1),stData.lidarDataArray(tempIndx,2),stData.lidarDataArray(tempIndx,3),'*','Color',[rand() rand() rand()]);
%     end
%     camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; view(-45, 15); 
%     axis([-stData.maxTreeWidth stData.maxTreeWidth -stData.maxTreeWidth stData.maxTreeWidth stData.mintreeHeight stData.maxtreeHeight]);
% end

function isSuccess = makeDEMFromLAS(exePath,exeName,inFilePath,inFileName,outFilePath,outFileName)
    command = char(strcat(exePath,exeName,{' -i '},inFilePath,inFileName,{' -o '},outFilePath,outFileName,{' -step 0.05 '})); %-keep_class 1 2 3 4 5 6 7
    [status,cmdout] = system(command); 
    isSuccess = status;
end

function makeTXTFromLASMAIN(lidarDataArray,outFilepath)        
    makeTXTFromLAS(lidarDataArray, strcat(outFilepath,'TempFile'),'.txt');
end

function isSuccess = makeLASFromTXT(exePath,exeName,inFilePath,inFileName,outFilePath,outFileName)
    command = char(strcat(exePath,exeName,{' -i '},inFilePath,inFileName,{' -o '},outFilePath,outFileName,{''}));
    [status,cmdout] = system(command); 
    isSuccess = status;
end

function makeTXTFromLAS(lidarDataArray,fullFileName,format)
    dlmwrite(strcat(fullFileName,format),lidarDataArray,'delimiter',' ','precision',10);
end

function [clusteredPointArray, retClustIndxArray, peakPointIndexArr, peakPointArr]= do_Kmeans(stData, slabbedCoordinates, numCluster,spoints, plotResult)

    clusteredPointArray = {}; peakPointArr = []; peakPointIndexArr = [];
    opts = statset('Display','final');
    [cidx, vect] = kmeans(slabbedCoordinates(:,1:2), numCluster, 'Distance','city', 'Options',opts, 'Start',spoints(:,1:2));
    cnt=1;
    dist_matrix  = squareform(pdist(vect, 'euclidean'));
    indices_to_set = true(size(dist_matrix));
    indices_to_set = triu(indices_to_set);
    dist_matrix(indices_to_set) = inf;
    [v1,v2] = find(and(dist_matrix < 1,dist_matrix ~=0));
    for i =1:1:size(v1,1)
        cidx(cidx==v1(i)) = cidx(v2(i))+1000;
    end

    uniquecidx = unique(cidx);
    index2Remove = [];
    for i = 1:1:size(uniquecidx,1)
        retClustIndxArray{i} = find(cidx==uniquecidx(i));
        if(size(retClustIndxArray{i},1) < 5)
           index2Remove = [index2Remove i];
        end
        aavsv =slabbedCoordinates(retClustIndxArray{i},1:4);
        clusteredPointArray{i} = stData.lidarDataArray(cidx==uniquecidx(i),1:3);
        [vald,idd] = max(slabbedCoordinates(retClustIndxArray{i},3));        
        peakPointIndex{i} = slabbedCoordinates(aavsv(idd,4),4);
        peakPointIndexArr = [peakPointIndexArr; peakPointIndex{i}];
        peakPointArr = [peakPointArr;  stData.lidarDataArray(ismember(stData.lidarDataArrayComplete(:,7), peakPointIndex{i}),:)];
    end

    if(plotResult)
        figure(1);
        for jindx = 1:1:size(retClustIndxArray,2)
            idx = retClustIndxArray{jindx};
            %plot3(slabbedCoordinates(idx,1), slabbedCoordinates(idx,2), slabbedCoordinates(idx,3),'.','MarkerSize',10,'Color',[rand() rand() rand()]);            
            hold on
            
            if(~isempty(peakPointIndex{jindx}))
                ppnt = slabbedCoordinates(ismember(slabbedCoordinates(:,4),peakPointIndex{jindx}),1:3);
                plot3(ppnt(1), ppnt(2), ppnt(3),'.','MarkerSize',40,'Color',[0 0 1]);
            end
            
        end
    end
    clusteredPointArray =  removerows(clusteredPointArray','ind',index2Remove)';
    retClustIndxArray =  removerows(retClustIndxArray','ind',index2Remove)';
    peakPointIndexArr = removerows(peakPointIndexArr,'ind',index2Remove)';
    peakPointArr =  removerows(peakPointArr,'ind',index2Remove);
end

% function [retClustIndxArray, peakPointIndex]= do_Kmeans(stData, slabbedCoordinates, numCluster,spoints, plotResult)
% 
%     opts = statset('Display','final');
%     [cidx, vect] = kmeans(slabbedCoordinates(:,1:2), numCluster, 'Distance','city', 'Options',opts, 'Start',spoints(:,1:2));
%     cnt=1;
%     dist_matrix  = squareform(pdist(vect, 'euclidean'));
%     indices_to_set = true(size(dist_matrix));
%     indices_to_set = triu(indices_to_set);
%     dist_matrix(indices_to_set) = inf;
%     [v1,v2] = find(and(dist_matrix < 1,dist_matrix ~=0));
%     for i =1:1:size(v1,1)
%         cidx(cidx==v1(i)) = cidx(v2(i))+1000;
%     end
% 
%     %val = unique(cidx);    
%     for i = 1:1:size(vect,1)
%         retClustIndxArray{i} = find(cidx==i);
%         aavsv =slabbedCoordinates(retClustIndxArray{i},1:4);
%         [vald,idd] = max(slabbedCoordinates(retClustIndxArray{i},3));        
%         peakPointIndex{i} = slabbedCoordinates(aavsv(idd,4),4);
%     end
%     if(plotResult)
%         for jindx = 1:1:size(retClustIndxArray,2)
%             idx = retClustIndxArray{jindx};
%             plot3(slabbedCoordinates(idx,1), slabbedCoordinates(idx,2), slabbedCoordinates(idx,3),'.','MarkerSize',10,'Color',[rand() rand() rand()]);            
%         end
%     end
% end

function F = myfun(p0,data)
    k = 1:1:size(data,1);
    F = power(data(k,1),2)./p0(1)^2 + power(data(k,2),2)./p0(2)^2 - data(k,3)./p0(3);
end

function [sectorLength, disttocentre0] = getSectorLengthold(refPoint, inputPoint, centrePoint)
    disttocentre0 = pdist2(inputPoint(:,1:2),centrePoint);
    disttocentre = pdist2(refPoint(:,1:2),centrePoint);    
    refAxisPointXY = getPointAtDistance(refPoint(:,1:2),centrePoint,disttocentre);    
    angleBwVextors = getAngleBetweenVectorsold(refAxisPointXY, inputPoint(:,1:2));
    sectorLength = 2*pi*disttocentre*(angleBwVextors/360);
end

function [sectorLength, disttocentre0] = getSectorLength(refPoint, inputPoint, centrePoint)
    disttocentre0 = pdist2(inputPoint(:,1:2),centrePoint);
    disttocentre = pdist2(refPoint(:,1:2),centrePoint);    
    refAxisPointXY = getPointAtDistance(refPoint(:,1:2),centrePoint,disttocentre);    
    angleBwVextors = real(getAngleBetweenVectors(refAxisPointXY, inputPoint(:,1:2)));
    sectorLength = 2*pi*disttocentre*(angleBwVextors/360);
end

function angleBwVextors = getAngleBetweenVectorsold(vec1, vec2)
    sdf = vec1*vec2';
    dfsd = (norm(vec1)*norm(vec2));
    angleBwVextors =  acosd((vec1*vec2')/(norm(vec1)*norm(vec2)));
end

function angleBwVextors = getAngleBetweenVectors(vec1, vec2)
    numwe = (vec1*vec2')';    
    denm = norm(vec1).*sqrt(vec2(:,1).^2 + vec2(:,2).^2);
    angleBwVextors =  acosd(numwe./denm);
end


function retPoint = getPointAtDistance(point1,point2,distVal)
    d = sqrt((point1(2)-point1(1))^2 + (point2(2) - point1(2))^2);
    r = distVal / d;
    x3 = r * point2(1) + (1 - r) * point1(1); % find point that divides the segment
    y3 = r * point2(2) + (1 - r) * point1(2); % into the ratio (1-r):r
    retPoint = [x3 y3];
end

function retVal = extremeBranchPoint(lidarDataArray, startPoint)
    otherPoints = lidarDataArray(:,(1:2));
    [~,I] = pdist2(otherPoints,startPoint(1:2),'euclidean','largest',1);
    %[maxV,maxI] = max(D);
    retVal = lidarDataArray(I,(1:3));
end

function [prosConnectedComp, retRand] = getProspConnectCompData(lidarDataArray, branchTipArray, lineSlopeAndIndexArr,includePointDistance)
    
    retRand = rand(size(branchTipArray,1),3);
    stemPointArr  = lineSlopeAndIndexArr(:,4:6);

    lidarAssignIndexArr = zeros(size(branchTipArray,1),3);
    for ii = 1: 1: size(lidarDataArray,1)
        distArr=zeros(size(branchTipArray,1),3);
        for jj = 1: 1: size(branchTipArray,1)
            distArr(jj,1) = ii; 
            distArr(jj,2) = jj; 
            distArr(jj,3) = getDist2Point(branchTipArray(jj,(1:3)), stemPointArr(jj,:), lidarDataArray(ii,(1:3)));
        end
        %distArr(distArr(:,3)>2,3) = inf;
        distArr = sortrows(distArr, 3); % Sort by distance column
        lidarAssignIndexArr(ii,:) = distArr(1,1:3);
    end
    
    lidarAssignIndexArr = sortrows(lidarAssignIndexArr,2);
    
    for jj = 1: 1: size(branchTipArray,1)        
        prosConnectedComp{jj} = lidarDataArray(lidarAssignIndexArr(and(lidarAssignIndexArr(:,2)==jj,abs(lidarAssignIndexArr(:,3))<includePointDistance),1),(1:3));
        tmpGrp = prosConnectedComp{jj};

       % hold on;        
        %plot3(tmpGrp(:,1),tmpGrp(:,2),tmpGrp(:,3),'*','Color',[retRand(jj,1) retRand(jj,2) retRand(jj,3)]);
       % hold off;
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

function getDist = getDist2Point(point1, point2, extPoint)
    x = extPoint; %some point
    a = point1; %segment points a,b
    b = point2;
    
    v = b-a;
    w = x-a;
    c1 = w*v';
    c2 = v*v';
    if(c1<=0)
       getDist =  sqrt(sum((x - a).^ 2)); 
       return;
    elseif (c2<=c1)
       getDist =  sqrt(sum((x - b).^ 2)); 
       return;
    else
        b1= c1/c2;
        pb= a + b1*v;
        getDist =  sqrt(sum((x - pb).^ 2));   
    end
  
end

function [extremeLiDARDataArray, extremeLiDARDataIndex, extremeLIDARDataAngle] = locateExteriorPoints(lidarDataArray, zDiv, maxAngle, angleIncrement)
    count = 1;
    extremeLiDARDataArray = zeros(zDiv*maxAngle/angleIncrement,size(lidarDataArray,2)); %cell2table(cell(0,5));
    extremeLiDARDataIndex = zeros(zDiv*maxAngle/angleIncrement,1);
    extremeLIDARDataAngle = zeros(zDiv*maxAngle/angleIncrement,1); 
    maxZ = max(lidarDataArray(:,3)); minZ = min(lidarDataArray(:,3));
    
    for zStep = minZ :(maxZ - minZ)/zDiv : maxZ-(maxZ - minZ)/zDiv;
        for angle = 0 : angleIncrement : maxAngle;
            
         lidarDataArray1 = rotateLiDARdata(lidarDataArray(:,(1:3)), angle*(pi/180));
            clippedLiDARpointsIndices = find(and(lidarDataArray(:,3) >= zStep , lidarDataArray(:,3) < zStep+(maxZ - minZ)/zDiv));
            
            %To record right extreme point in a particualr ZDiv
            [~,tmpIndex] = max(lidarDataArray1(clippedLiDARpointsIndices,1));
            temparr = lidarDataArray(clippedLiDARpointsIndices,:);
            extremeLiDARDataArray(count,:) = temparr(tmpIndex,:); % lidarDataArray1(tmpIndex,:);
            extremeLiDARDataIndex(count,:) =  zStep;
            extremeLIDARDataAngle(count,:) = angle; %??
            count = count + 1;
            
            % To record left extreme point in a particualr ZDiv
            [~,tmpIndex] = min(lidarDataArray1(clippedLiDARpointsIndices,1));
            temparr = lidarDataArray(clippedLiDARpointsIndices,:);            
            extremeLiDARDataArray(count,:) = temparr(tmpIndex,:); 
            extremeLiDARDataIndex(count,:) =  zStep;
            extremeLIDARDataAngle(count,:) = 180 + angle; %??
            
            count = count + 1;  
        end
    end
    
    [~, ind] = unique(extremeLiDARDataArray,'rows');    
    extremeLiDARDataArray = extremeLiDARDataArray(ind,:);
    extremeLiDARDataIndex = extremeLiDARDataIndex(ind,:);
    extremeLIDARDataAngle = extremeLIDARDataAngle(ind,:);
    
    isValidArr = zeros(size(extremeLiDARDataArray,1),1);
    for i = 1:1:size(extremeLiDARDataArray,1)
        isValidArr(i,1) = isValidExtremePoint(lidarDataArray(:,(1:3)),extremeLiDARDataArray(i,(1:3)), 2); %2 is threshold value
    end

    extremeLiDARDataArray = extremeLiDARDataArray(isValidArr==1,:);
    extremeLiDARDataIndex = extremeLiDARDataIndex(isValidArr==1,:);
    extremeLIDARDataAngle = extremeLIDARDataAngle(isValidArr==1,:);
    
end

function retPointIndx = getProminalPoints(startPoint, otherPoints, distTheshold, NumNearestPoints, maxIter, retPointIndx)
    [dist,Indx] = pdist2(otherPoints,startPoint,'euclidean','Smallest',NumNearestPoints);
    retPointIndx = [retPointIndx;Indx];
    maxIter = maxIter-1;
    for i = 2:1:2
        if(dist(i) < distTheshold) && maxIter > 0
            retPointIndx = [retPointIndx;getProminalPoints(otherPoints(retPointIndx(size(retPointIndx,1)),:), otherPoints, distTheshold, NumNearestPoints, maxIter, retPointIndx)];
        else
            retPointIndx = [];
        end
    end  
end

function retConnectedComp = getConnectedComponets(lidarDataSlab, seedPointsArray)

    %outerLoopCnt= 1;
    NumNearestPoints = 6; distTheshold = 1.5;
    %connectComponentArray = zeros(NumNearestPoints*noOfIterations,3,(compEndIndex-compStartIndex));
    %seletedIndices = [];
    
    for clstStartPointIndex = 1:1:size(seedPointsArray,1)  %length(extremeLiDARDataArray)-50:1:length(extremeLiDARDataArray)-51+noPoints
        %cntd =1; clusteredLiDARDataArray = zeros(NumNearestPoints*noOfIterations,3); 
        startPoint = seedPointsArray(clstStartPointIndex,(1:3));
        %startPoint1 = startPoint;
         % [double(startPoint(1)) double(startPoint(2)) double(startPoint(3))]
        otherPoints = lidarDataSlab(:,(1:3));
        indexList{clstStartPointIndex} = [];
        
        
        %(find(and(lidarDataArray(:,3) > startPoint(3)-10, lidarDataArray(:,3) < startPoint(3)+10)),(1:3));

        indexList{clstStartPointIndex} = unique(getProminalPoints(startPoint, otherPoints, distTheshold, NumNearestPoints, 10, []));
        ff =0;
        %[dist,I] = pdist2(otherPoints,startPoint,'euclidean','Smallest',NumNearestPoints);
        %temotherPoints = [otherPoints(I,:) lidarDataSlab(I,6)];
        %sortLiDArray = sortrows(temotherPoints,[4]); %order by 3rd column in decresing order.            

        %a1a = otherPoints(I(2),:); b2b = startPoint1(1,(1:3));
        %dst = pdist([a1a;b2b],'euclidean');
        %if(dst<2)
        %    clusteredLiDARDataArray(cntd:cntd+(NumNearestPoints-1),:) = otherPoints(I,:);                
        %    cntd = cntd + NumNearestPoints;
        %end
        %startPoint = sortLiDArray(NumNearestPoints,(1:3));
        %otherPoints(I,:) = 0;

        
        %clusteredLiDARDataArray = clusteredLiDARDataArray(clusteredLiDARDataArray(:,3)~=0,:);
        %connectComponentArray{outerLoopCnt} = clusteredLiDARDataArray;%zeros(NumNearestPoints*noOfIterations,3,(compEndIndex-compStartIndex));        
        %outerLoopCnt = outerLoopCnt + 1;

    end
    
    retConnectedComp = indexList;
end

function plotConnectedComponents(connectComponentArray)
    for i = 1:1:size(connectComponentArray,3);
        clusteredLiDARDataArray = connectComponentArray(:,:,i);
        plot3(clusteredLiDARDataArray(:,1), clusteredLiDARDataArray(:,2), clusteredLiDARDataArray(:,3),'*','Color',[0.5 0 0]);
        hold on; camproj perspective; rotate3d on; view(3), axis vis3d; axis equal; axis on;
    end
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

function rotateLiDARdata = rotateLiDARdata(lidarDataArray, teeta)
    for i = 1: length(lidarDataArray)
        lidarPosVector = [lidarDataArray(i,1) lidarDataArray(i,2) lidarDataArray(i,3)];
        lidarNewPosVector = rotateVector(lidarPosVector,teeta);
        lidarDataArray(i,1)  = lidarNewPosVector(1);
        lidarDataArray(i,2)  = lidarNewPosVector(2);
        lidarDataArray(i,3)  = lidarNewPosVector(3);
    end
    rotateLiDARdata = lidarDataArray;
end


function [lidarBranchT, remainingIndx] = removeRepetitions(extremeLiDARDaraArray, CutoffDistBetwVect)
   
    DistBwExtrmPoints = tril(squareform(pdist(extremeLiDARDaraArray(:,1:3),'euclidean')));
    %dims = [size(DistBwExtrmPoints,1),size(DistBwExtrmPoints,2)];   
    %SIndx = find(DistBwExtrmPoints < CutoffDistBetwVect & DistBwExtrmPoints > 0);
    
    [I,~] = ind2sub(size(DistBwExtrmPoints),find(DistBwExtrmPoints < CutoffDistBetwVect & DistBwExtrmPoints > 0));
    indx = unique(I);
    lidarBranchTips = removerows(extremeLiDARDaraArray,'ind',indx);
    
    lidarBranchT = lidarBranchTips(:,1:3);
    remainingIndx = lidarBranchTips(:,4); 
end

% function [extrLiDARDataArray, remainingIndx] = removeRepetitionsOld(extremeLiDARDaraArray,  CutoffDistBetwVect)
%     indx = [];   
%         
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
% 
%     extremeLiDARDaraArray = removerows(extremeLiDARDaraArray,'ind',indx);
%     
%     extrLiDARDataArray = extremeLiDARDaraArray(:,1:3);
%     remainingIndx = extremeLiDARDaraArray(:,4);    
% end

function retVector = rotateVector(inputVector, teeta)
    rotationMatrix = [cos(teeta) -sin(teeta) 0; sin(teeta) cos(teeta) 0; 0 0 1];
    retVector = rotationMatrix*inputVector';
end

function retVector = rotateVectorXZ(inputVector, teeta)
    rotationMatrix = [1 0 0; 0 cos(teeta) -sin(teeta); 0 sin(teeta) cos(teeta)];
    retVector = rotationMatrix*inputVector';
end

function plotLiDARData(stData, isPlotOn, isPlotOnAll)
    if(and(isPlotOn, isPlotOnAll))
        figure(1); 
        subplot(1,2,1); 
        hold on;
        camproj perspective; rotate3d on; axis vis3d; axis equal; grid on; axis on; view(-45, 15); axis([-stData.maxTreeWidth stData.maxTreeWidth -stData.maxTreeWidth stData.maxTreeWidth stData.mintreeHeight stData.maxtreeHeight]);
        plot3(stData.lidarDataArray(:,1), stData.lidarDataArray(:,2), stData.lidarDataArray(:,3),'.','Color',[0 0.5 0]);
        hold on;    
        %if(true)       
        plot3([0 0],[0 0],[0 stData.maxtreeHeight],'-o','Color',[1 0 0]);
        %end
        hold off;
        xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 16); 
        ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 16); 
        zlabel('Tree Height','Fontname', 'Times New Roman' ,'FontSize', 16);
        set(gca,'XLim',[-stData.maxTreeWidth stData.maxTreeWidth]);set(gca,'YLim',[-stData.maxTreeWidth stData.maxTreeWidth]);set(gca,'ZLim',[0 stData.treeHeight]);
        set(gca,'XTick',-stData.maxTreeWidth:2:stData.maxTreeWidth); set(gca,'YTick',-stData.maxTreeWidth:2:stData.maxTreeWidth); set(gca,'ZTick',0:2:stData.treeHeight)
    end
end

function returnData =  LoadData(fullFileName)
    returnData = lasdata(fullFileName);
end

function retTable = write2table(lasFile)
    retTable = zeros(size(lasFile.x,1),5);
    retTable(:,1) = lasFile.x;
    retTable(:,2) = lasFile.y;
    retTable(:,3) = lasFile.z;
    retTable(:,4) = get_classification(lasFile);    
    retTable(:,5) = lasFile.get_return_number;
    retTable(:,6) = get_intensity(lasFile);
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

function isValid = isValidExtremePoint(lidarDataArray,startPoint,thresholdvalue)
    [nearPoint,~] = pdist2(lidarDataArray,startPoint,'euclidean','Smallest',2);          
    isValid = 1;
    if(nearPoint(2)>thresholdvalue)
        isValid = 0;
    end    
end

function retDistWeights = calculateDistWeights(lidarDataArray, neighbourhoodSize)
   [nearPoint,~] = pdist2(lidarDataArray(:,(1:3)),lidarDataArray(:,(1:3)),'euclidean','Smallest',100);  
    nearPoint(nearPoint>=neighbourhoodSize) = 0; % make dist 0 for points > threshold  neighbourhoodDist  
    nearPoint = nearPoint(2:end,:); % remove first line because it is distance to point itself.
    neigDensityArr = sum(nearPoint~=0,1)'; % find non zero elenents in the array for every column
    retDistWeights = neigDensityArr/norm(neigDensityArr); % normalization
end

function retMaxXYZ = findMaxHeightXY(lidarDataArray)
    [~,indx] = max(lidarDataArray(:,3));
     retMaxXYZ = lidarDataArray(indx,1:3);
end

function Z = Peak_detection(sc, demFullPath)
    dsmImage = imread(demFullPath);
    dsmImage(isnan(dsmImage)==1)=1000;
    dsmImage(dsmImage==-9999)=1000;
    minVal=min(min(dsmImage));
    dsmImage(dsmImage==1000)=0;
    dsmImage(dsmImage~=0) = dsmImage(dsmImage~=0) - minVal;
    dsmImage = (dsmImage/max(max(dsmImage)))*255;

    %j=1;
    for i = 20:1:20
    myfilter = fspecial('gaussian',[i i], 1);
    myfiltereddsmImage = imfilter(dsmImage, myfilter, 'replicate');
    %peakIndices=FastPeakFind(myfiltereddsmImage);
    peakIndices=FastPeakFind(myfiltereddsmImage, 0, myfilter ,2, 1, 'abcc');

    cnt=1;
    peakValues=0;
    for k =1:2:length(peakIndices)-1
        rowInd(cnt) = peakIndices(k+1);
        colInd(cnt) = peakIndices(k);
        peakValues(cnt) = myfiltereddsmImage(peakIndices(k+1), peakIndices(k));
        cnt=cnt+1;
    end
   
    sz = size(dsmImage);
     
    aam = max(sc(:,1)); %29.27
    aami = min(sc(:,1)); %6.06
    bbm = max(sc(:,2)); %12.57
    bbmi = min(sc(:,2)); %-12.57
    ccm = max(sc(:,3)); %54.00
    ccmi = min(sc(:,3));% 50.31

    
    peakIndices(1:2:end) = ((peakIndices(1:2:end)/sz(2))*(aam-aami))+aami;
    peakIndices(2:2:end) = ((peakIndices(2:2:end)/sz(1))*(bbm - bbmi))+ bbmi;
    znew = ((peakValues/max(peakValues))*(ccm-ccmi)) + ccmi;
            
    peakIndicesNew = [peakIndices(1:2:end) peakIndices(2:2:end)];
    
    Z = [peakIndicesNew znew'];
%     if(true==true)
%     
%         subplot(1,2,j);
%         surf(double(myfiltereddsmImage)); %alpha(.3);
%         hold on;
%         sdf =peakIndices(1:2:end);
%         stem3(peakIndices(1:2:end),peakIndices(2:2:end),Z,'Marker','s',...
%                             'MarkerEdgeColor','m',...
%                             'MarkerFaceColor','r');
% 
%         j=j+1;
%     end
    
    end

end
% cropn = imadjust(crop,[0.3 1],[0 1]);       
%addpath(genpath('C:\My_Files\Topic_of_Research\1_LiDAR_Research\Lidar_Tools\LASData_Matlab_Tool'))
% set(0,'format', 'long');\

%if(abs(a)/abs(b)>5)
%    tmp = min(a,b);
 %   a = tmp;
 %   b = tmp;
%elseif(abs(b)/abs(a)>5)
   % tmp = min(abs(a),abs(b));
   % a = tmp;
   % b = tmp;
%end

%f = @(x)parameterfun_modello(x,stData.lidarDataArr,1,stData.crownHeight,stData.maxTreeWidth);
%lsqnonlin(f,1,5,30)

%         figure(2); pData = true;     
%         if(plotResult)
%             for jindx = 1:1:size(retClustIndxArray,2)
%                 idx = retClustIndxArray{jindx};
%                 %hold on;
%                % if(jindx ==52)
%                     plot3(slabbedCoordinates(idx,1), slabbedCoordinates(idx,2), slabbedCoordinates(idx,3),'.','MarkerSize',20,'Color',[rand() rand() rand()]);            
%               % end
%             end
%         hold off;
%         camproj perspective; rotate3d on; axis vis3d; axis normal; axis on; view(-45, 15); 
%         %axis([-stData.maxTreeWidth stData.maxTreeWidth -stData.maxTreeWidth stData.maxTreeWidth stData.mintreeHeight stData.maxtreeHeight]);
%         axis([5 30 -10 10 0 5]);
%         end