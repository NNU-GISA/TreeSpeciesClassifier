function egffeatureArr = getEGFs(inFilepath,species, allFiles,plotOn)

    count1 = 1;
    egffeatureArr = zeros(size(allFiles,1),6);
    
    for file = allFiles'
                
        fullFileName = strcat(inFilepath,file.name); % get full file name

        %-----------------------------------------------------------------%
        % -- The below part of the code does LiDAR data Preprocessing --- %
        % ----------------------------------------------------------------%
        
        %Read LiDAR data
        data = LoadData(fullFileName);
        
        % LiDAR data Pre-processing and basic tree information retrieval
        stData = dataPreProcessing(data);
        
       % try
         %   selectedIndices = randperm(length(stData.lidarDataArray),1000); 
          %  lidarDataArray = stData.lidarDataArray(selectedIndices,:);
       % catch        
            lidarDataArray = stData.lidarDataArray;            
      %  end

        lidarDataArray2 = lidarDataArray;
        treeHeight = stData.treeHeight;        
        retMaxXY = stData.retMaxXYZ;       
        crownHeight =  stData.crownHeight;
       
        heighdtMax = max(lidarDataArray(:,3));
        heightMin = min(lidarDataArray(:,3));
        ht = heighdtMax - heightMin;

        X= stData.lidarDataArrSMall(:,(1:3));
        x0 = [0 ;0; ht];
        a0 =[0;0;25];
        phi0 = 1*(atan(stData.maxTreeWidth/ht)*(180/pi));
        r0 = stData.maxTreeWidth;
        tolp = 0.001;
        tolg = 0.001;
       
        
        X(:,3)= -X(:,3) + max(X(:,3));
       % Cone around the z-axis, point at the origin   
       [x0n, an, phin, rn, d, sigmah, conv, Vx0n, Van, uphin, ... 
                urn, GNlog, a, R0, R] = lscone(X, x0, a0, phi0, r0, ... 
                tolp, tolg);

        if(conv == 0)
            r = 3; m = 9.34;
        else
            r = rn; m = phin;
        end
       
        [R,A] = meshgrid(linspace(0,r*2,110),linspace(0,2*pi,41));
        X1 = R .* cos(A) + x0n(1); %+ sigmaXc;    
        Y1 = R .* sin(A) + x0n(2);
        Z1 = -R*(m*0.5)+ treeHeight;
        
        fet0 = (1/3)*pi*((r)^2)*ht; % volume of cone % fet0 = (1/3)*pi*(maxx/2)^2*ht; %

        if(plotOn)
            clf; % clear figures   
            hold on; rotate3d on; axis equal;
            axis([-6 6 -6 6 0 treeHeight]);
            xlabel('X Axis'); ylabel('Y Axis'); zlabel('Tree Height');
            alpha(0.5); view(3);   
            rotate3d on; %axis equal;
            plot3(lidarDataArray2(:,1),lidarDataArray2(:,2),lidarDataArray2(:,3),'.','Color',[0 0.5 0])
            hold on;
            mesh(X1,Y1,Z1);
            hold on; 
            plot3([retMaxXY(1) retMaxXY(1)],[retMaxXY(2) retMaxXY(2)],[0 ht],'-o');
            camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-45, 15); axis([-6 6 -6 6  0 treeHeight]);
        end
        
        [K,V] = boundary(lidarDataArray2(:,1),lidarDataArray2(:,2),lidarDataArray2(:,3),0); %0.5:ar  0.2:la   %0.5:pc    0.5:ab
        
        if(plotOn)
            %trisurf(K,lidarDataArray2(:,1),lidarDataArray2(:,2),lidarDataArray2(:,3))
            axis([-10 10 -10 10 0 ht]);
            camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(-45, 15); axis([-10 10 -10 10 0 treeHeight]); alpha(0.5);
            hold on;
        end

        d1 = zeros(size(lidarDataArray2,1),3);        
        for i = 1:1:size(lidarDataArray2,1)
            tmp = getDistOfPoint2PlaneNew(lidarDataArray2(K(:,1),1:3), lidarDataArray2(K(:,2),1:3), lidarDataArray2(K(:,3),1:3), lidarDataArray2(i,1:3));
            tmp(tmp==0) = 999999999999999;
            %[aa,bb]=min(tmp);
            [d1(i,2), d1(i,3)]=min(tmp);
            d1(i,1) = i;
        end

        fet1 = V/size(K,1);
        fet2 = (V-fet0)/V;
        if(conv == 0)
            if(strcmp(species,'ar'))
                fet3 = 0.8747; fet4 = 1.0797;  % when the algorithm does not converge (average value)
            elseif(strcmp(species,'la'))
                fet3 = 0.9876; fet4 = 1.2357;  % when the algorithm does not converge (average value)
            elseif(strcmp(species,'pc'))
                fet3 = 0.2488; fet4 = 0.3222;  % when the algorithm does not converge (average value)
            else
                fet3 = 0.7237; fet4 = 0.9124;  % when the algorithm does not converge (average value)
            end
        else            
            fet3 = sum(abs(d))/size(d,1);
            fet4 = sigmah;
        end
        fet5 = crownHeight/treeHeight;
        fet6 = mean(d1(:,2));
        egffeatureArr(count1,:) = [fet1 fet2 fet3 fet4 fet5 fet6];
        count1 = count1 + 1;
        
        if(plotOn)
            pause(1); 
        end
    end
      
    csvwrite('csvlistconefit.csv',egffeatureArr);
    
end 

function retCoefVec = eqPlanenew(A,B,C)
    a = (B(:,2)-A(:,2)).*(C(:,3)-A(:,3)) - (C(:,2)-A(:,2)).*(B(:,3)-A(:,3));
    b = -((B(:,1)-A(:,1)).*(C(:,3)-A(:,3)) - (C(:,1)-A(:,1)).*(B(:,3)-A(:,3)));
    c = (B(:,1)-A(:,1)).*(C(:,2)-A(:,2)) - (C(:,1)-A(:,1)).*(B(:,2)-A(:,2));  
    d = -(a.*A(:,1) + b.*A(:,2) + c.*A(:,3));
    retCoefVec = [a;b;c;d];
end

function distOfPoint2Plane =  getDistOfPoint2PlaneNew(p1,p2,p3,pext)
    retCoefVec = eqPlanenew(p1,p2,p3);
    %distOfPoint2Plane = abs(retCoefVec(1)*pext(1) + retCoefVec(2)*pext(2) + retCoefVec(3)*pext(3) + retCoefVec(4))/sqrt(sumsqr(retCoefVec(1:3)));
    dfg = size(retCoefVec,1)/4;
    dd =  reshape(retCoefVec,dfg,4);
    
    num = diag(abs(repmat([pext 1], dfg,1)*dd'));
    dem = sqrt(dd(:,1).^2 + dd(:,2).^2 + dd(:,3).^2); 
    
    distOfPoint2Plane = num./dem;
end

function dataAtt = dataPreProcessing(singleTreeLiDARdata)

% Write normalized data to a table for performance improvement 
        lidarDataArr = normalizeLiDARData(write2table(singleTreeLiDARdata));
        
        % Tree height
        dataAtt.maxTreeWidth = max(max(lidarDataArr(:,1)), max(lidarDataArr(:,2))); 
        dataAtt.mintreeHeight = min(lidarDataArr(:,3));
        dataAtt.maxtreeHeight = max(lidarDataArr(:,3));        
        dataAtt.treeHeight = dataAtt.maxtreeHeight - dataAtt.mintreeHeight;        
        dataAtt.htDeduction = dataAtt.treeHeight*0.1;
        
        % Crown height
        lidarDataArr = lidarDataArr(find(and(lidarDataArr(:,3) > min(lidarDataArr(:,3))+ dataAtt.htDeduction, lidarDataArr(:,3) < max(lidarDataArr(:,3)))),:);
        dataAtt.minCrownHeight = getMinCrownHeight(lidarDataArr);
        dataAtt.maxCrownHeight = max(lidarDataArr(:,3));        
        dataAtt.crownHeight = dataAtt.maxCrownHeight - dataAtt.minCrownHeight;       
        
        if(size(lidarDataArr,1)>2000)
            dataAtt.lidarDataArrSMall = lidarDataArr(randperm(size(lidarDataArr,1),2000),:);
        else
            dataAtt.lidarDataArrSMall = lidarDataArr;
        end
        
        %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 > 0.1),:);                
        %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 < 12),:);
        
        % Identify the center point of the tree (in top view).
        dataAtt.retMaxXYZ = findMaxHeightXY(lidarDataArr);
        
        
        dataAtt.lidarDataArray = lidarDataArr;
        %dataAtt.lidarDataArray = lidarDataArray(randperm(size(lidarDataArray,1),9000),:); % for reducing the no. of samples for testing only.
        
        dataAtt.lidarDataDensityArr = zeros(size(lidarDataArr));% calculateDistWeights(lidarDataArr,30);
        dataAtt.index = 1:1:size(lidarDataArr,1);
        %dataAtt.lidarDataDensityArr = calculateDistWeights(lidarDataArr,10);

        
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
end

function WriteLidar2Txt(lidarDataArray,fullFileName,format)
    dlmwrite(strcat(fullFileName,format),lidarDataArray,'delimiter',' ','precision',10);
end

function isSuccess = makeLASFromTXT(exePath,exeName,inFilePath,inFileName,outFilePath,outFileName)
    command = char(strcat(exePath,exeName,{' -i '},inFilePath,inFileName,{' -o '},outFilePath,outFileName,{''}));
    [status,cmdout] = system(command); 
    isSuccess = status;
end

function retMaxXYZ = findMaxHeightXY(lidarDataArray)
    [~,indx] = max(lidarDataArray(:,3));
     retMaxXYZ = lidarDataArray(indx,1:3);
end