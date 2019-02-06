function mainFunctionOfMethord()
    % Define the filepath here and file name variable in this code. Give this filepath and name as input to the getConiferIGFs and
    % getConiferEGFs function and (remove filepath from them) to avoid repetition of filepath and filename data.

    %IGFs_proposed = 6 Internal Geometric Features derived from the proposed model
    %IGFs_proposed = 6 Internal Geometric Features derived from the SoA model
    %EGFs = 6 External Geometric Features

    % Provide full path of the folder containing the .las files
    % The full path should not have '-' (hyphens) in it %
    
    % Proved full name of specific .las file if needed, else for considering 
    % all the .las files the folder (on by one) provide '*.las'.
    
    FolderPath = 'C:\My_Files\1_PhD_Research\1_PhD_Research_Topic\1_LiDAR_Forest_Applications\3_Research_Files\1_Conifer_Species_Classification\1_Matlab_Files\1_ConiferSpeciesDetection\Conifer_classifiation_Files\';
    InputFilePath = strcat(FolderPath,'LiDARDataSingleTrees\temp_test_sample\');
    speciesFolder = {'ar' 'la' 'pc' 'ab'}; % LiDAR data of the four different species.    
    OutFolder = strcat(FolderPath,'generatedCSVs\');
    
    IGFsProposed=[]; IGFSoAs=[]; EGFs=[]; Labels=[]; 

    % Show 3D plots true/false
    plotOn = true;
    
    for i=1:size(speciesFolder,2)
        % Provide specific .las file if needed. Else for considering all 
        % the .las files the folder (on by one) provide '*.las'
        inFilepath = char(strcat(InputFilePath,speciesFolder(i),'\')); % cycles through all the input folders
        files = dir(strcat(inFilepath,'*.las'));
        disp(strcat('Currently working on:', {' '}, num2str(speciesFolder{i})) );
        
        % Get the IGFs from the Proposed method        
        if(plotOn)
            figure(1);
        end        
        [A, ASliced, numBranches] = getConiferIGFs_Proposed(inFilepath, speciesFolder{i}, files, plotOn);
        %A = getConiferIGFs_Method2(inFilepath, speciesFolder{i}, files, plotOn);
        
        % 6 IGFs from entire tree (proposed method)   
        IGFsProposed = [IGFsProposed; A(:,1:6)];
        B = getConiferIGFsSoA(inFilepath,files,numBranches, plotOn); 
        IGFSoAs = [IGFSoAs; B];
        
        %Get the EGFs from the respective Cone fit & Convex hull
        if(plotOn)
            figure(3);
        end
        C = getEGFs(inFilepath,speciesFolder{i},files, plotOn); 
        EGFs = [EGFs; C];

        % Label column values
        Labels=[Labels; repmat(i,size(C,1),1)];
    end
    
    % Feature value normalization
    Norm_IGFsPCA = normalize(IGFsProposed); 
    Norm_IGFKMs = normalize(IGFSoAs);
    Norm_EGFs = normalize(EGFs);    
    
    % Generate Excel Files with feature values for experiments
    Normalized_features_set0=[Norm_EGFs Labels];  % For experiments 1
    Normalized_features_set1=[Norm_IGFKMs Labels];  % For experiments 2 
    Normalized_features_set2=[Norm_IGFsPCA Labels];  % For experiments 3
    Normalized_features_set3=[Norm_IGFKMs, Norm_EGFs, Labels];  % For experiments 4
    Normalized_features_set4=[Norm_IGFsPCA, Norm_EGFs, Labels];  % For experiments 5
     
    % Print result in CSV
    csvwrite(strcat(OutFolder,'Norm_EGFs.csv'),Normalized_features_set0);
    csvwrite(strcat(OutFolder,'Norm_IGFKMs.csv'),Normalized_features_set1);
    csvwrite(strcat(OutFolder,'Norm_IGFsPCA.csv'),Normalized_features_set2);
    csvwrite(strcat(OutFolder,'Norm_IGFKM_EGFs.csv'),Normalized_features_set3);
    csvwrite(strcat(OutFolder,'Norm_IGFsPCA_EGFs.csv'),Normalized_features_set4);
    
    % The analysis of the results is done using the Matlab code in the 
    % SVM_Classifiers folder.
end

% Function for normalization column by column
function Norm=normalize(A)
    Norm_A=[];
    for i=1:size(A,2)
        Norm_col=(A(:,i)-min(A(:,i)))/(max(A(:,i))-min(A(:,i)));
        Norm_A=[Norm_A,Norm_col];
    end
    Norm=Norm_A;
end

% function mainFunctionOfMethord()
%     % Define the filepath here and file name variable in this code. Give this filepath and name as input to the getConiferIGFs and
%     % getConiferEGFs function and (remove filepath from them) to avoid repetition of filepath and filename data.
% 
%     %IGFs_proposed = Internal Geometric Features derived from the proposed model
%     %IGFs_proposed = Internal Geometric Features derived from the SoA model
%     %EGFs = Internal Geometric Features derived
% 
%     % Provide full path of the folder containing the .las files
%     % The full path should not have '-' (hyphens) in it %
%     
%     % Proved full name of specific .las file if needed, else for considering 
%     % all the .las files the folder (on by one) provide '*.las'.
%     
%     FolderPath = 'C:\My_Files\1_PhD_Research\1_PhD_Research_Topic\1_LiDAR_Forest_Applications\3_Research_Files\1_Conifer_Species_Classification\1_Matlab_Files\1_ConiferSpeciesDetection\Conifer_classifiation_Files\LiDARDataSingleTrees\';
%     InputFilePath = strcat(FolderPath,'temp_test_sample\');
%     speciesFolder = {'ar' 'la' 'pc' 'ab'}; % LiDAR data of the four different species.    
%     OutFolder = strcat(FolderPath,'generatedCSVs\');
%     
%     IGFsProposed=[]; IGFSoAs=[]; EGFs=[]; Labels=[]; 
%     %Intensity_Features=[];
%     %IGF_4Div=[];  IGF_7Div=[];  IGF_10Div=[]; IGF_13Div=[];
%     
%     % Show 3D plots true/false
%     plotOn = false;
%     
%     for i=1:size(speciesFolder,2)
%         % Provide specific .las file if needed. Else for considering all 
%         % the .las files the folder (on by one) provide '*.las'
%         inFilepath = char(strcat(InputFilePath,speciesFolder(i),'\')); % cycles through all the input folders
%         files = dir(strcat(inFilepath,'*.las'));
%         disp( strcat('Currently working on:', {' '}, num2str(speciesFolder{i})) );
%         
%         % Get the IGFs from Proposed method        
%         if(plotOn)
%             figure(1);
%         end
%         
%         [A, ASliced, numBranches] = getConiferIGFs_Proposed(inFilepath, speciesFolder{i}, files, plotOn);
%         %A = getConiferIGFs_Method2(inFilepath, speciesFolder{i}, files, plotOn);
%         
%         % 6 IGFs from entire tree (proposed method)   
%         IGFsProposed = [IGFsProposed; A(:,1:6)];
%         
%         %2 intensity feature from entire tree (proposed method)
%         % Intensity_Features = [Intensity_Features; A(:,7:8)];
%         
%         %6-IGFs + 2-IntensityFeatures from 4 tree divisions (i.e. 8 * 4 = 32 features)
%         % IGF_4Div = [IGF_4Div; ASliced{1}];
%         %6-IGFs + 2-IntensityFeatures from 7 tree divisions (i.e. 8 * 7 = 56 features)
%         % IGF_7Div = [IGF_7Div; ASliced{2}];
%         %6-IGFs + 2-IntensityFeatures from 10 tree divisions (i.e. 8 * 10 = 80 features)
%         % IGF_10Div = [IGF_10Div; ASliced{3}];
%         %6-IGFs + 2-IntensityFeatures from 13 tree divisions (i.e. 8 * 13 = 104 features)
%         % IGF_13Div = [IGF_13Div; ASliced{4}];
% 
% %         %Get the IGFs from State-of-the-art method
% %         if(plotOn)
% %             figure(2);
% %         end  
%          B = getConiferIGFsSoA(inFilepath,files,numBranches, plotOn); 
%          IGFSoAs = [IGFSoAs; B];
%         
%         %Get the EGFs from Cone fit & Convex hull
%         if(plotOn)
%             figure(3);
%         end
%         C = getEGFs(inFilepath,speciesFolder{i},files, plotOn); 
%         EGFs = [EGFs; C];
% 
%         % Label column values
%         Labels=[Labels; repmat(i,size(C,1),1)];
%     end
%     
%     % Feature value normalization
%     Norm_IGFsPCA = normalize(IGFsProposed);
%     %Norm_Int_Features = normalize(Intensity_Features);    
%     %Norm_IGF_4Div = normalize(IGF_4Div);
%     %Norm_IGF_7Div = normalize(IGF_7Div);
%     %Norm_IGF_10Div = normalize(IGF_10Div);
%     %Norm_IGF_13Div = normalize(IGF_13Div);    
%     Norm_IGFKMs = normalize(IGFSoAs);
%     Norm_EGFs = normalize(EGFs);    
%     
%     % Generate Excel Files with feature values for experiments
%     Normalized_features_set0=[Norm_EGFs Labels];  % For experiments 1
%     Normalized_features_set1=[Norm_IGFKMs Labels];  % For experiments 2 
%     Normalized_features_set2=[Norm_IGFsPCA Labels];  % For experiments 2
%     Normalized_features_set3=[Norm_IGFKMs, Norm_EGFs, Labels];  % For experiments 3
%     Normalized_features_set4=[Norm_IGFsPCA, Norm_EGFs, Labels];  % For experiments 4
%     %Normalized_features_set13=[Norm_IGF_13Div, Norm_EGFs, Labels];  % For experiments 5
%     %Normalized_features_set1=[Norm_EGFs Labels];  % For experiments 6
%     %Normalized_features_set2=[Norm_IGFsPCA Labels];  % For experiments 7
%     %Normalized_features_set3=[Norm_Int_Features Labels];  % For experiments 8
%      
%     % Print result in CSV
%     
%     
%     
%     csvwrite(strcat(OutFolder,'Norm_EGFs.csv'),Normalized_features_set0);
%     csvwrite(strcat(OutFolder,'Norm_IGFKMs.csv'),Normalized_features_set1);
%     csvwrite(strcat(OutFolder,'Norm_IGFsPCA.csv'),Normalized_features_set2);
%     csvwrite(strcat(OutFolder,'Norm_IGFKM_EGFs.csv'),Normalized_features_set3);
%     csvwrite(strcat(OutFolder,'Norm_IGFsPCA_EGFs.csv'),Normalized_features_set4);
% end
% 
% % Function for normalization column by column
% function Norm=normalize(A)
%     Norm_A=[];
%     for i=1:size(A,2)
%         Norm_col=(A(:,i)-min(A(:,i)))/(max(A(:,i))-min(A(:,i)));
%         Norm_A=[Norm_A,Norm_col];
%     end
%     Norm=Norm_A;
% end


