clear; close all; clc;

%% Description 
% This code parses a selection of MAT files in the /mat folder for surface
% meshes, and remeshes them at a desired point spacing. The remeshed
% surfaces are exported to the /post folder. 
%
% This code requires the GIBBON MATLAB toolbox 
% <www.gibboncode.org>

%% Plotting settings
lineWidth=0.5;

%% Control parameters

% Path names
projectFolder = fileparts(fileparts(mfilename('fullpath'))); %Main code path
loadFolder=fullfile(projectFolder,'data','BodyParts3D','mat'); %The MAT loading folder
saveFolder=fullfile(projectFolder,'data','BodyParts3D','post'); %The MAT saving folder for processed data

fileNames_FMA={'FMA16586','FMA24474','FMA24477','FMA24480','FMA24486'};

pointSpacings=[4 4 4 3 3];

optionStructRemesh.disp_on=0; % Turn off command window text display

saveOn=0;

%%
    
for q=1:1:numel(fileNames_FMA)
    
    %% Import mesh
    fileName_FMA=fileNames_FMA{q};    
    fileName_mat=fullfile(loadFolder,[fileName_FMA,'.mat']);
    
    model=load(fileName_mat);
    F=model.faces;
    V=model.vertices;
    
    %% Remesh
    optionStructRemesh.pointSpacing=pointSpacings(q); %Set desired point spacing
    [Fn,Vn]=ggremesh(F,V,optionStructRemesh);
    
    %% Visualisation
    
    cFigure;
    gtitle([fileName_FMA,' ',model.preferredName])
    subplot(1,2,1); 
    title('Raw');
    gpatch(F,V,'w','r',1,lineWidth);
    axisGeom;
    camlight headlight;
    
    subplot(1,2,2);
    title('Remeshed');
    gpatch(Fn,Vn,'w','g',1,lineWidth);
    axisGeom;
    camlight headlight;
    gdrawnow;
        
    %% Saving
    
    modelNew.source=model;
    modelNew.faces=Fn;
    modelNew.vertices=Vn;
    modelNew.pointSpacing=pointSpacings(q);
    
    if saveOn==1
        %Create save name with lowercase letters and underscores instead of spaces
        saveNameMesh=regexprep(lower(model.preferredName),' ','_');        
        saveName_mat=fullfile(saveFolder,[saveNameMesh,'.mat']);
        save(saveName_mat,'-struct','modelNew')
    end
    
end