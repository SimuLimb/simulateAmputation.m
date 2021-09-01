clear; close all; clc; 

%% Description 
% This code processes the MAT file for the skin surface (found in the /mat
% folder), since it contains both the inner and outer skin surfaces. The
% inner surface is removed and the outer surface is processed to produce a
% single "water tight" mesh. The surface is also remeshed using ggremesh.
% The remeshed surface is exported to the /post folder. 
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

fileName_FMA='FMA7163';
saveNameMesh='skin_leg_right';

pointSpacing=4; 
cutHeight=675;
saveOn=1;

%%

fileName_mat=fullfile(loadFolder,[fileName_FMA,'.mat']);

model=load(fileName_mat);
F=model.faces;
V=model.vertices;
        
%%

logicVertices=V(:,1)<0 & V(:,3)<=cutHeight+mean(patchEdgeLengths(F,V));
logicFaces=any(logicVertices(F),2);
logicFaces=triSurfLogicSharpFix(F,logicFaces); 
Fs=F(logicFaces,:); 
[Fs,Vs]=patchCleanUnused(Fs,V);

optionStruct.outputType='label';
[G,~,groupSize]=tesgroup(Fs,optionStruct);
[~,indMax]=max(groupSize);
logicKeep=G==indMax;
Fs=Fs(logicKeep,:); 
[Fs,Vs]=patchCleanUnused(Fs,Vs);

snapTolerance=mean(patchEdgeLengths(Fs,Vs))/100; 
n=vecnormalize([0 0 1]); %Normal direction to plane
P=[0 0 cutHeight]; %Point on plane
[Fc,Vc,~,logicSide]=triSurfSlice(Fs,Vs,[],P,n,snapTolerance);
[Fc,Vc]=patchCleanUnused(Fc(logicSide,:),Vc);

%% Construct alpha shape
% This deteriorates surface quality (bridges concave regions) but is a
% termporary "quick-and-dirty" work-around to obtain outer skin surface
% only. 

shp = alphaShape(Vc,max(patchEdgeLengths(Fc,Vc)),'HoleThreshold',500,'RegionThreshold',1);
[Fa,Va] = boundaryFacets(shp); %Get boundary faces of alpha shape
[Fa,Va] = patchCleanUnused(Fa,Va); %Remove unused vertices

%% Remesh alpha shape using |ggremesh|
optionStructRemesh.pointSpacing=pointSpacing; %Set desired point spacing
optionStructRemesh.disp_on=1; % Turn off command window text display
[Fn,Vn]=ggremesh(Fa,Va,optionStructRemesh);

%% Visualisation

cFigure; 
subplot(1,2,1); hold on; 
gpatch(F,V,'w','none',0.25);
gpatch(Fn,Vn,'bw','none');
axisGeom; 
camlight headlight; 

subplot(1,2,2); hold on; 
gpatch(Fn,Vn,'bw','k',1,lineWidth);
axisGeom; 
camlight headlight; 
gdrawnow; 

%% Store model in structure

modelNew.source=fileName_FMA;
modelNew.faces=Fn;
modelNew.vertices=Vn;
modelNew.pointSpacing=pointSpacing; 
modelNew.cutHeight=cutHeight;

%% Save model

if saveOn==1
    saveName_mat=fullfile(saveFolder,[saveNameMesh,'.mat']);
    save(saveName_mat,'-struct','modelNew')
end
