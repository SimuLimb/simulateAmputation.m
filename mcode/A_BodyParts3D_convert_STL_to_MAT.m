clear; close all; clc;

%% Description
% This code parses all STL files in the /stl folder and exports them as MAT
% files to the /mat folder.
%
% This code requires the GIBBON MATLAB toolbox
% <www.gibboncode.org>

%% Control parameters

% Path names
projectFolder = fileparts(fileparts(mfilename('fullpath'))); %Main code path
loadFolder=fullfile(projectFolder,'data','BodyParts3D','stl'); %The STL loading folder
saveFolder=fullfile(projectFolder,'data','BodyParts3D','mat'); %The MAT saving folder

%%

%Get list of all STL files
fileList = dir(fullfile(loadFolder,['*','stl']));
fileList={fileList(1:end).name};

Dc=readtable(fullfile(projectFolder,'data','BodyParts3D','FMA_ID_label_obj.csv'));
P=lower(Dc.FMAID);

%Convert all STL files to MAT
numFiles=numel(fileList); %Number of files to parse
for q=1:1:numFiles %Loop over all files
    %Get file name
    [~,fileNameClean,~]=fileparts(fileList{q}); %File name without path or extension
    
    %Prepare STL file name
    fileNameLoad=fullfile(loadFolder,fileList{q}); %STL file name
    
    %Prepare save file name
    fileNameSave=fullfile(saveFolder,[fileNameClean,'.mat']); %MAT file name
    
    %Get FMAID code
    FMAID=sscanf(fileNameClean,'FMA%d');
    
    %Get preferred name
    logicFMAID=P==FMAID;
    preferredName=Dc.preferredName{logicFMAID};
    
    %Parse STL import
    disp(['Parsing file ', num2str(q),' of ',num2str(numFiles),', ',sprintf('%3.0f',round(100*q/numFiles)),'% done, ',fileNameClean,', ',preferredName])    
    try %MATLAB's stlread
        TR = stlread(fileNameLoad);
        v=TR.Points;
        f=TR.ConnectivityList;
    catch %GIBBON's importer
        [stlStruct] = import_STL_bin(fileNameLoad);
        v=stlStruct.solidVertices{1};
        f=stlStruct.solidFaces{1};
    end
    
    %Merge vertices
    [f,v]=mergeVertices(f,v);
    
    %Build model structure
    model.sourceName=fileList{q};
    model.preferredName=preferredName;
    model.faces=f;
    model.vertices=v;
    
    %Export MAT file
    save(fileNameSave,'-struct','model');    
end

%% 
% _*SimuLimb footer text*_ 
% 
% License: <https://github.com/SimuLimb/simulateAmputation.m/blob/main/LICENSE>
%   
% Copyright (C) 2006-2021 Kevin Mattheus Moerman and the SimuLimb contributors
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
