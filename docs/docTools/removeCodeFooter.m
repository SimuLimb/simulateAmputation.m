clear; close all; clc;

%% Path names

toolboxPath=fileparts(fileparts(fileparts(mfilename('fullpath'))));
docPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'docs');
libPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'mcode');
docToolsPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'docs','docTools');

pathNames={docToolsPath,libPath,docPath};

%%

fileExtension='.m'; %Extension to remove boilerplate from

footerTargetStart='% _*SimuLimb footer text*_ '; %Target for removal

%N.B. Removal is up to end from the start so make sure the target is
%appropriately set!!!!!

%% Removing footer text

numPaths=numel(pathNames); 
for q_path=1:1:numPaths
    pathName=pathNames{q_path};     
    files = dir(fullfile(pathName,'*.m'));
    files={files(1:end).name};
    files=sort(files(:));
    numFiles=numel(files);    
    for q_file=1:1:numFiles
        fileName=fullfile(pathName,files{q_file});
        [T_now]=txtfile2cell(fileName);
        targetStartIndex = find(strcmp(footerTargetStart,T_now));       
        if ~isempty(targetStartIndex)            
            targetStartIndex=targetStartIndex(end); %Keep last occurance
            T_now=T_now(1:targetStartIndex-1-1); %NB -1 is used to remove %% above target 
            cell2txtfile(fileName,T_now,0,0);
        end
    end
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
