clear; close all; clc;

%% Description 
% This code processes skin and bone surface meshes for the lower limb to
% simulate amputation. Amputation is "simulated" in the sense that a 3D
% surface mesh is obtained with a shape that may be expected from a typical
% amputation. The current code is for transtibial amputation. 
%
% This code requires the GIBBON MATLAB toolbox 
% <www.gibboncode.org>

%% Control parameters 

projectFolder = fileparts(fileparts(mfilename('fullpath')));
loadFolder=fullfile(projectFolder,'data','BodyParts3D','post'); 
saveFolder_stl=fullfile(projectFolder,'data','BodyParts3D','post_stl'); 
saveFolder=loadFolder; 
saveNameGeom='BodyParts3D_right_leg_transtibial_amp';
select_amputation_case='tt';%'tt';'tf' % where the trans-tibial(tt); trans-femoral (tf)

switch select_amputation_case
     case 'tf'
         fileNames={'right_femur','right_leg_skin','right_leg_muscles'};
     case 'tt'
        fileNames={'right_femur','right_tibia','right_fibula','right_patella','right_leg_skin','right_leg_muscles'};
end
    
saveOn=0;

%The percentage for trans-femoral (above knee) amputation
amputationPercentageFemur=25;
%Distance parameters for trans-tibial (below knee) amputation
amputationDistanceTibia=180; 
amputationDistanceFibula=amputationDistanceTibia-30;
amputationDistanceSkin=amputationDistanceTibia-25;
distalExcess=50;
distalExcess_muscle_flap=25;
boneRoundFactor=1;
radiusEnd=distalExcess+(amputationDistanceTibia-amputationDistanceSkin);
taperHeigth=100; 
taperThreshold=15; %Regions more distance than threshold are tappered
taperFraction=0.3; 
amputationDistances=[amputationDistanceTibia amputationDistanceFibula amputationDistanceSkin amputationDistanceSkin];
tibiaDistalCutRadialFactor=0.5;
numBezierPoints=50;
nCut=vecnormalize([0 0 1]); %Normal direction to plane
topCropOffset=10;
topCropOffset_tf=30;
%% Loading surfaces into cell array
%Allocate cell arrays
FT=cell(1,numel(fileNames));
VT=cell(1,numel(fileNames));
CT=cell(1,numel(fileNames));

%Loop over all surfaces
for q=1:1:numel(fileNames)
    
    %Import mesh
    fileName=fileNames{q};    
    fileName_mat=fullfile(loadFolder,[fileName,'.mat']);
    model=load(fileName_mat);
    F=model.faces; %Faces
    V=model.vertices; %Vertices 
    C=q*ones(size(F,1),1); %Color label
    
    %Store mesh data in cell array
    FT{q}=F;
    VT{q}=V;   
    CT{q}=C;
end

%% 
% Visualize surface and landmarks
cFigure; hold on; 
gpatch(FT,VT,CT,'none',0.5);
axisGeom; camlight headlight; 
colormap gjet; icolorbar; 
gdrawnow; 


%% Compute landmarks for tt amputation
switch select_amputation_case
    case 'tf'
        %% Compute the femur length for tf amputation
        V_bone_femur=VT{1};
        [Min_femur,I_min_femur]=min(V_bone_femur(:,3));
        [Max_femur,I_max_femur]=max(V_bone_femur(:,3));
        femurLength = Max_femur - Min_femur; %femur length
        c=[1 2 3]; %Indices (=labels) for surfaces to cut

        %%
        % Visualize surface and landmarks
        cFigure; hold on;
        gpatch(FT,VT,CT,'none',0.5);
        plotV(V_bone_femur(I_min_femur,:),'k.','MarkerSize',35);
        axisGeom; camlight headlight;
        colormap gjet; icolorbar;
        gdrawnow;

    case 'tt'
        P_patella_centroid=triSurfCentroid(FT{4},VT{4}); %Centroid of the patella
        c=[2 3 5 6]; %Indices (=labels) for surfaces to cut
        %%
        % Visualize surface and landmarks
        cFigure; hold on;
        gpatch(FT,VT,CT,'none',0.5);
        plotV(P_patella_centroid,'k.','MarkerSize',35);
        axisGeom; camlight headlight;
        colormap gjet; icolorbar;
        gdrawnow;
end

%% Cut surfaces
%Allocate cell arrays
FT_amp=FT; 
VT_amp=VT; 
CT_amp=CT; 

for q=1:numel(c)  

    switch select_amputation_case 
        case 'tt'
            cutLevelNow=P_patella_centroid(:,3)-amputationDistances(q);
        case 'tf'
            cutLevelNow=Min_femur+amputationPercentageFemur/100*femurLength;
    end
  
    %Get surface
    F=FT_amp{c(q)};
    V=VT_amp{c(q)};
    
    %Use triSurfSlice to process cut
    snapTolerance=mean(patchEdgeLengths(F,V))/100;    
    P=[0 0 cutLevelNow]; %Point on plane
    [Fc,Vc,~,logicSide]=triSurfSlice(F,V,[],P,nCut,snapTolerance);
    [Fc,Vc]=patchCleanUnused(Fc(~logicSide,:),Vc); %Clean-up mesh
    Cc=c(q)*ones(size(Fc,1),1); %Adjusted color data (shorter list after cut)
    
    %Store processed mesh data in cell arrays
    FT_amp{c(q)}=Fc;
    VT_amp{c(q)}=Vc;
    CT_amp{c(q)}=Cc;
end

%%
% Visualization

cFigure; hold on; 
title('Cut features');
% gpatch(FT,VT,'w','none',0.1);
gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
axisGeom; camlight headlight; 
colormap gjet; icolorbar; 
gdrawnow; 

%% Set up skin taper parameterisation
switch select_amputation_case
    case 'tf'
        %Femur
        Ft=FT_amp{1};
        Vt=VT_amp{1};
        %Skin
        Fs=FT_amp{2};
        Vs=VT_amp{2};
        %Muscle
        Fm=FT_amp{3};
        Vm=VT_amp{3};        
    case 'tt'
        %Tibia
        Ft=FT_amp{2};
        Vt=VT_amp{2};
        %Skin
        Fs=FT_amp{5};
        Vs=VT_amp{5};
        %Muscle
        Fm=FT_amp{6};
        Vm=VT_amp{6};     
end

        
%Get cut bone end curve
Ebt=patchBoundary(Ft,Vt);
indBt=edgeListToCurve(Ebt);
indBt=indBt(1:end-1);
P_end_centroid=mean(Vt(indBt,:),1);
%Skin
D=Vs(:,3)-(min(Vs(:,3))+taperHeigth);
D(D>0)=0;
D=abs(D);
D=D./max(D);
D=D.^2;
        
[Dp,indMin]=minDist(Vs(:,[1 2]),Vt(indBt,[1 2]));

Dp(Dp<=taperThreshold)=taperThreshold;
Dp=Dp-min(Dp);
Dp=Dp./max(Dp);

R=Vt(indBt(indMin),[1 2])-Vs(:,[1 2]);
R=R.*D.*taperFraction.*Dp;
        
%% Process skin taper
Vs1=Vs; %Initialize as original
Vs(:,[1 2])=Vs(:,[1 2])+R; %Push to produce taper

%Overide skin surface with tapered surface
switch select_amputation_case
    case 'tf'
        FT_amp{2}=Fs;
        VT_amp{2}=Vs;
        
    case 'tt'
        %Tibia
        FT_amp{5}=Fs;
        VT_amp{5}=Vs;
end

cFigure; hold on;
title('taper morphing');
gpatch(FT_amp,VT_amp,'w','none',0.5);
% gpatch(Fs,Vs,'w','none',0.5);
quiverVec(Vs1,R);
colormap gjet; colorbar;
axisGeom; camlight headlight;
gdrawnow;


%% Process skin distal end closure
P_end=P_end_centroid-[0 0 distalExcess];
Ebs=patchBoundary(Fs,Vs);
[Fsb,Vsb,indBs,Nd,P1,P2,P3,XB,YB,ZB]=smoothBezierClose(Fs,Vs,Ebs,P_end,numBezierPoints);

%%
% Visualization 

cFigure; hold on; 
title('Distal end closure');
gpatch(FT_amp,VT_amp,'w','none',1);
gpatch(Fs,Vs,'w','none',0.25);
hp=gpatch(Fsb,Vsb,Vsb(:,3),'none',1); hp.FaceColor='interp';
% gpatch(Ebs,Vs,'none','b',1,3);

quiverVec(Vs(indBs,:),Nd(indBs,:),radiusEnd/4,'k');

plotV(P1,'r.-','MarkerSize',15,'LineWidth',2);
plotV(P2,'g.-','MarkerSize',15,'LineWidth',2);
plotV(P3,'b.-','MarkerSize',15,'LineWidth',2);
for q=1:size(XB,2)
    plotV([XB(:,q) YB(:,q) ZB(:,q)],'k.-','LineWidth',0.5,'MarkerSize',5);
end

plotV(P_end,'k.','MarkerSize',25);

axisGeom; camlight headlight; 
gdrawnow; 


%% Merge skin components and remesh 
pointSpacing=mean(patchEdgeLengths(Fs,Vs));
[Fs,Vs]=joinElementSets({Fs,Fsb},{Vs,Vsb});
[Fs,Vs]=mergeVertices(Fs,Vs);

optionStructRemesh.pointSpacing=pointSpacing; %Set desired point spacing
optionStructRemesh.disp_on=0; % Turn off command window text display
[Fs,Vs]=ggremesh(Fs,Vs,optionStructRemesh);

%%
switch select_amputation_case
    case 'tf'
        FT_amp{2}=Fs;
        VT_amp{2}=Vs;
        CT_amp{2}=2*ones(size(Fs,1),1);
        
    case 'tt'
        FT_amp{5}=Fs;
        VT_amp{5}=Vs;
        CT_amp{5}=5*ones(size(Fs,1),1);
end


cFigure; 
subplot(1,2,1); hold on; 
% gpatch(FT,VT,'w','none',0.25);
gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
axisGeom; camlight headlight; 
colormap gjet; icolorbar; 

subplot(1,2,2); hold on; 
gpatch(FT_amp,VT_amp,'w','k',1);
axisGeom; camlight headlight; 
gdrawnow; 


%% Process muscle taper
Dm=Vm(:,3)-(min(Vm(:,3))+taperHeigth);
Dm(Dm>0)=0;
Dm=abs(Dm);
Dm=Dm./max(Dm);
Dm=Dm.^2;
        
[Dp,indMin]=minDist(Vm(:,[1 2]),Vt(indBt,[1 2]));

Dp(Dp<=taperThreshold)=taperThreshold;
Dp=Dp-min(Dp);
Dp=Dp./max(Dp);

Rm=Vt(indBt(indMin),[1 2])-Vm(:,[1 2]);
Rm=Rm.*Dm.*taperFraction.*Dp;

Vm1=Vm; %Initialize as original
Vm(:,[1 2])=Vm(:,[1 2])+Rm; %Push to produce taper

%Overide muscle surface with tapered surface
switch select_amputation_case
    case 'tf'
        %Muscle
        FT_amp{3}=Fm;
        VT_amp{3}=Vm;
        
    case 'tt'
        %Muscle
        FT_amp{6}=Fm;
        VT_amp{6}=Vm;
end

cFigure; hold on;
title('taper morphing');
gpatch(FT_amp,VT_amp,'w','none',0.5);
quiverVec(Vm1,Rm);
colormap gjet; colorbar;
axisGeom; camlight headlight;
gdrawnow;

%% Process muscle distal end closure
P_end=P_end_centroid-[0 0 distalExcess_muscle_flap];
Ebm=patchBoundary(Fm,Vm);
[Fmb,Vmb,indBm,Nd,P1,P2,P3,XB,YB,ZB]=smoothBezierClose(Fm,Vm,Ebm,P_end,numBezierPoints);
%%
% Visualization 
cFigure; hold on; 
title('Distal end closure');
gpatch(FT_amp,VT_amp,'w','none',0.25);
gpatch(Fm,Vm,'w','none',1);
hp=gpatch(Fmb,Vmb,Vmb(:,3),'none',1); hp.FaceColor='interp';

quiverVec(Vm(indBm,:),Nd(indBm,:),radiusEnd/4,'k');

plotV(P1,'r.-','MarkerSize',15,'LineWidth',2);
plotV(P2,'g.-','MarkerSize',15,'LineWidth',2);
plotV(P3,'b.-','MarkerSize',15,'LineWidth',2);
for q=1:size(XB,2)
    plotV([XB(:,q) YB(:,q) ZB(:,q)],'k.-','LineWidth',0.5,'MarkerSize',5);
end

plotV(P_end,'k.','MarkerSize',25);

axisGeom; camlight headlight; 
gdrawnow; 

%% Merge skin components and remesh 
pointSpacing=mean(patchEdgeLengths(Fm,Vm));
[Fm,Vm]=joinElementSets({Fm,Fmb},{Vm,Vmb});
[Fm,Vm]=mergeVertices(Fm,Vm);

optionStructRemesh.pointSpacing=pointSpacing; %Set desired point spacing
optionStructRemesh.disp_on=0; % Turn off command window text display
[Fm,Vm]=ggremesh(Fm,Vm,optionStructRemesh);

%%
switch select_amputation_case
    case 'tf'
        FT_amp{3}=Fm;
        VT_amp{3}=Vm;
        CT_amp{3}=3*ones(size(Fm,1),1);
        
    case 'tt'
        FT_amp{6}=Fm;
        VT_amp{6}=Vm;
        CT_amp{6}=6*ones(size(Fm,1),1);
end


%%
cFigure; 
subplot(1,2,1); hold on; 
% gpatch(FT,VT,'w','none',0.25);
gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
axisGeom; camlight headlight; 
colormap gjet; icolorbar; 

subplot(1,2,2); hold on; 
gpatch(FT_amp,VT_amp,'w','k',0.5);
axisGeom; camlight headlight; 
gdrawnow; 

%% Process femur taper
switch select_amputation_case
    case 'tf'
        Ft=FT_amp{1};
        Vt=VT_amp{1};             
    case 'tt'
        Ft=FT_amp{2};
        Vt=VT_amp{2};        
end

Ebs=patchBoundary(Ft,Vt);
indBs=edgeListToCurve(Ebs);
indBs=indBs(1:end-1); 

P_mid=mean(Vt(indBs,:),1);

radiusEnd=boneRoundFactor*mean(minDist(Vt(indBs,:),P_mid))/2;
P_end=P_mid-[0 0 radiusEnd];

[Fsb,Vsb,indBs,Nd,P1,P2,P3,XB,YB,ZB]=smoothBezierClose(Ft,Vt,Ebs,P_end,numBezierPoints);

%%
% Visualization
cFigure; hold on; 
gpatch(Ft,Vt,'w','none',0.25);
gpatch(Ebs,Vt,'none','b',1,3);
hp=gpatch(Fsb,Vsb,'w','k',1,1); 
% hp.FaceColor='interp';
plotV(P_end_centroid,'k.','MarkerSize',35);

quiverVec(Vt(indBs,:),Nd(indBs,:),radiusEnd/4,'k');

plotV(P1,'r.-','MarkerSize',15,'LineWidth',2);
plotV(P2,'g.-','MarkerSize',15,'LineWidth',2);
plotV(P3,'b.-','MarkerSize',15,'LineWidth',2);
for q=1:size(XB,2)
    plotV([XB(:,q) YB(:,q) ZB(:,q)],'k.-','LineWidth',0.5,'MarkerSize',5);
end

axisGeom; camlight headlight; 
gdrawnow; 

%% Merge bone components and remesh 
pointSpacing=mean(patchEdgeLengths(Ft,Vt));
[Ft,Vt]=joinElementSets({Ft,Fsb},{Vt,Vsb});
[Ft,Vt]=mergeVertices(Ft,Vt);

optionStructRemesh.pointSpacing=pointSpacing; %Set desired point spacing
optionStructRemesh.disp_on=0; % Turn off command window text display
[Ft,Vt]=ggremesh(Ft,Vt,optionStructRemesh);

%%
switch select_amputation_case
    case 'tf'
        FT_amp{1}=Ft;
        VT_amp{1}=Vt;
        CT_amp{1}=ones(size(Ft,1),1);
    case 'tt'
        FT_amp{2}=Ft;
        VT_amp{2}=Vt;
        CT_amp{2}=2*ones(size(Ft,1),1);
end

%%

cFigure; 
subplot(1,2,1); hold on; 
% gpatch(FT,VT,'w','none',0.25);
gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
axisGeom; camlight headlight; 
colormap gjet; icolorbar; 

subplot(1,2,2); hold on; 
gpatch(FT_amp,VT_amp,'w','none',0.5);
gpatch(Ft,Vt,'w','k',1);
axisGeom; camlight headlight; 
gdrawnow; 

cFigure; 
hold on; 
% gpatch(FT,VT,'w','none',0.25);
gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
axisGeom; camlight headlight; 
colormap gjet; icolorbar; 
gdrawnow; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process fibula taper      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch select_amputation_case        
    case 'tt'
        %% FIBULA
        Ft=FT_amp{3};
        Vt=VT_amp{3};
        
        Ebs=patchBoundary(Ft,Vt);
        indBs=edgeListToCurve(Ebs);
        indBs=indBs(1:end-1);
        
        P_mid=mean(Vt(indBs,:),1);
        
        radiusEnd=boneRoundFactor*mean(minDist(Vt(indBs,:),P_mid))/2;
        P_end=P_mid-[0 0 radiusEnd];
        
        [Fsb,Vsb,indBs,Nd,P1,P2,P3,XB,YB,ZB]=smoothBezierClose(Ft,Vt,Ebs,P_end,numBezierPoints);
        
        %%
        % Visualization
        cFigure; hold on;
        gpatch(Ft,Vt,'w','none',0.25);
        gpatch(Ebs,Vt,'none','b',1,3);
        hp=gpatch(Fsb,Vsb,'w','k',1,1);
        % hp.FaceColor='interp';
        %plotV(P_end_centroid,'k.','MarkerSize',35);
        quiverVec(Vt(indBs,:),Nd(indBs,:),radiusEnd/4,'k');
        plotV(P1,'r.-','MarkerSize',15,'LineWidth',2);
        plotV(P2,'g.-','MarkerSize',15,'LineWidth',2);
        plotV(P3,'b.-','MarkerSize',15,'LineWidth',2);
        for q=1:size(XB,2)
            plotV([XB(:,q) YB(:,q) ZB(:,q)],'k.-','LineWidth',0.5,'MarkerSize',5);
        end
        
        axisGeom; camlight headlight;
        gdrawnow;


        %% Merge bone components and remesh
        pointSpacing=mean(patchEdgeLengths(Ft,Vt));
        [Ft,Vt]=joinElementSets({Ft,Fsb},{Vt,Vsb});
        [Ft,Vt]=mergeVertices(Ft,Vt);
        
        optionStructRemesh.pointSpacing=pointSpacing; %Set desired point spacing
        optionStructRemesh.disp_on=0; % Turn off command window text display
        [Ft,Vt]=ggremesh(Ft,Vt,optionStructRemesh);
        
        %% TIBIA
        FT_amp{3}=Ft;
        VT_amp{3}=Vt;
        CT_amp{3}=3*ones(size(Ft,1),1);
        
        %%
        cFigure;
        subplot(1,2,1); hold on;
        % gpatch(FT,VT,'w','none',0.25);
        gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
        axisGeom; camlight headlight;
        colormap gjet; icolorbar;
        
        subplot(1,2,2); hold on;
        gpatch(FT_amp,VT_amp,'w','none',0.5);
        gpatch(Ft,Vt,'w','k',1);
        axisGeom; camlight headlight;
        gdrawnow; 
%         
%                 %%
%         cFigure;
%         hold on;
%         gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
%         axisGeom; camlight headlight;
%         colormap gjet; icolorbar;
%         gdrawnow; 
            
        %% Cut the top muscles surface
        %fileNames={'right_femur','right_tibia','right_fibula','right_patella','right_leg_skin','right_leg_muscles'};
        Fm=FT_amp{6};
        Vm=VT_amp{6};
        
        cutLevelNow=max(Vm(:,3))-topCropOffset_tf;
        %Use triSurfSlice to process cut
        snapTolerance=mean(patchEdgeLengths(Fm,Vm))/100;
        
        P=[0 0 cutLevelNow]; %Point on plane
        [Fc,Vc,~,logicSide]=triSurfSlice(Fm,Vm,[],P,-nCut,snapTolerance);
        [Fc,Vc]=patchCleanUnused(Fc(~logicSide,:),Vc); %Clean-up mesh
        
        % Remesh using ggremesh
        optionStruct3.pointSpacing=mean(patchEdgeLengths(Fm,Vm));
        optionStruct3.disp_on=0; % Turn off command window text display
        optionStruct3.pre.max_hole_area=100; %Max hole area for pre-processing step
        optionStruct3.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
        
        [Fc,Vc]=ggremesh(Fc,Vc,optionStruct3);
        
        Ebs=patchBoundary(Fc,Vc);
        indBs=edgeListToCurve(Ebs);
        indBs=indBs(1:end-1);
        
        Fc_muscle=Fc;
        Vc_muscle=Vc;
        indBs_muscle=indBs;   
        %% Cut the top femur surface
        %Get the surface of the femur
        Ff=FT_amp{1};
        Vf=VT_amp{1};
        
        %Use triSurfSlice to process cut
        snapTolerance=mean(patchEdgeLengths(Ff,Vf))/100;
        [Fc,Vc,~,logicSide]=triSurfSlice(Ff,Vf,[],P,-nCut,snapTolerance);
        [Fc,Vc]=patchCleanUnused(Fc(~logicSide,:),Vc); %Clean-up mesh
        
        % Remesh using ggremesh
        optionStruct3.pointSpacing=mean(patchEdgeLengths(Ff,Vf));
        optionStruct3.disp_on=0; % Turn off command window text display
        optionStruct3.pre.max_hole_area=100; %Max hole area for pre-processing step
        optionStruct3.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
        
        [Fc_femur,Vc_femur]=ggremesh(Fc,Vc,optionStruct3);
        Ebs_femur=patchBoundary(Fc_femur,Vc_femur);
        indBs_femur=edgeListToCurve(Ebs_femur);
        indBs_femur=indBs_femur(1:end-1);
        
        %%
        %%
        
        z=0.5*(mean(Vc_femur(indBs_femur,3))+mean(Vc_femur(indBs_femur,3)));
        
        Vc_muscle(indBs_muscle,3)=z;%muscle
        Vc_femur(indBs_femur,3)=z;%femur
        
        %%
        pointSpacing=mean(patchEdgeLengths(Fc_muscle,Vc_muscle));
        [Ftt,Vtt]=regionTriMesh3D({Vc_muscle(indBs_muscle,:),Vc_femur(indBs_femur,:)},pointSpacing,0,'linear');
        
        
        %%
        cFigure;
        gpatch(Fc_muscle,Vc_muscle,'gw','k',1);
        gpatch(Ff,Vf,'rw','k',1);
        gpatch(Ftt,Vtt,'bw','k',0.6);
        colormap gjet;
        axisGeom; camlight headlight;
        gdrawnow;
        %%
        %FT_amp{1}=fliplr(Fc_femur); %invert femur normals
        FT_amp{1}=Fc_femur; %invert femur normals
        VT_amp{1}=Vc_femur;
        CT_amp{1}=1*ones(size(Fc_femur,1),1);
        
        FT_amp{6}=Fc_muscle;
        VT_amp{6}=Vc_muscle;
        CT_amp{6}=6*ones(size(Fc_muscle,1),1);
        
        FT_amp{7}=Ftt;
        VT_amp{7}=Vtt;
        CT_amp{7}=7*ones(size(Ftt,1),1);
        
        %%
        
        cFigure; hold on;
        gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
        % patchNormPlot(FT_amp,VT_amp);
        axisGeom; camlight headlight;
        colormap gjet; icolorbar;
        gdrawnow;

        %Cut the top skin surface
        Fs=FT_amp{5};
        Vs=VT_amp{5};
        
        %Use triSurfSlice to process cut
        snapTolerance=mean(patchEdgeLengths(Fs,Vs))/100;
        P=[0 0 cutLevelNow]; %Point on plane
        [Fc,Vc,~,logicSide]=triSurfSlice(Fs,Vs,[],P,-nCut,snapTolerance);
        [Fc,Vc]=patchCleanUnused(Fc(~logicSide,:),Vc); %Clean-up mesh
        
        % Remesh using ggremesh
        optionStruct3.pointSpacing=mean(patchEdgeLengths(Fs,Vs));
        optionStruct3.disp_on=0; % Turn off command window text display
        optionStruct3.pre.max_hole_area=100; %Max hole area for pre-processing step
        optionStruct3.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
        
        [Fc,Vc]=ggremesh(Fc,Vc,optionStruct3);
        
        Ebs=patchBoundary(Fc,Vc);
        indBs=edgeListToCurve(Ebs);
        indBs=indBs(1:end-1);
        
        Fc_skin=Fc;
        Vc_skin=Vc;
        indBs_skin=indBs;
        
        
        %Get the surface of the muscle
        Fc_muscle=FT_amp{6};
        Vc_muscle=VT_amp{6};

        Ebs_muscle=patchBoundary(Fc_muscle,Vc_muscle);
        indBs_muscle=edgeListToCurve(Ebs_muscle);
        indBs_muscle=indBs_muscle(1:end-1);
        
        %%
        z=0.5*(mean(Vc_muscle(indBs_muscle,3))+mean(Vc_muscle(indBs_muscle,3)));
        
        Vc_skin(indBs_skin,3)=z;%skin
        Vc_muscle(indBs_muscle,3)=z;%muscle
        
        %%
        pointSpacing=mean(patchEdgeLengths(Fc_skin,Vc_skin));
        [Fsm,Vsm]=regionTriMesh3D({Vc_skin(indBs_skin,:),Vc_muscle(indBs_muscle,:)},pointSpacing,0,'linear');
        
        %%
        cFigure;
        gpatch(Fc_skin,Vc_skin,'gw','k',1);
        gpatch(Fc_muscle,Vc_muscle,'rw','k',1);
        gpatch(Fsm,Vsm,'bw','k',0.6);
        colormap gjet;
        axisGeom; camlight headlight;
        gdrawnow;
        
        
        %%
        FT_amp{5}=Fc_skin;
        VT_amp{5}=Vc_skin;
        CT_amp{5}=5*ones(size(Fc_skin,1),1);
        
        FT_amp{8}=Fsm;
        VT_amp{8}=Vsm;
        CT_amp{8}=8*ones(size(Fsm,1),1);
        
        %%
        cFigure; hold on;
        gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
        patchNormPlot(FT_amp{1},VT_amp{1});
        axisGeom; camlight headlight;
        colormap gjet; icolorbar;
        gdrawnow;  
        
        case 'tf'
        %fileNames={'right_femur','right_leg_skin','right_leg_muscles'};
        %Cut the top muscle surface
        Fm=FT_amp{3};
        Vm=VT_amp{3};
        
        cutLevelNow=max(Vm(:,3))-topCropOffset_tf;
        %Use triSurfSlice to process cut
        snapTolerance=mean(patchEdgeLengths(Fm,Vm))/100;
        
        P=[0 0 cutLevelNow]; %Point on plane
        [Fc,Vc,~,logicSide]=triSurfSlice(Fm,Vm,[],P,-nCut,snapTolerance);
        [Fc,Vc]=patchCleanUnused(Fc(~logicSide,:),Vc); %Clean-up mesh
        
        % Remesh using ggremesh
        optionStruct3.pointSpacing=mean(patchEdgeLengths(Fm,Vm));
        optionStruct3.disp_on=0; % Turn off command window text display
        optionStruct3.pre.max_hole_area=100; %Max hole area for pre-processing step
        optionStruct3.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
        
        [Fc,Vc]=ggremesh(Fc,Vc,optionStruct3);
        
        Ebs=patchBoundary(Fc,Vc);
        indBs=edgeListToCurve(Ebs);
        indBs=indBs(1:end-1);
        
        Fc_muscle=Fc;
        Vc_muscle=Vc;
        indBs_muscle=indBs;
        
        %Get the surface of the femur
        Ff=FT_amp{1};
        Vf=VT_amp{1};
        
        %Use triSurfSlice to process cut
        snapTolerance=mean(patchEdgeLengths(Ff,Vf))/100;
        %P=[0 0 cutLevelNow]; %Point on plane
        [Fc,Vc,~,logicSide]=triSurfSlice(Ff,Vf,[],P,-nCut,snapTolerance);
        [Fc,Vc]=patchCleanUnused(Fc(~logicSide,:),Vc); %Clean-up mesh
        
        % Remesh using ggremesh
        optionStruct3.pointSpacing=mean(patchEdgeLengths(Ff,Vf));
        optionStruct3.disp_on=0; % Turn off command window text display
        optionStruct3.pre.max_hole_area=100; %Max hole area for pre-processing step
        optionStruct3.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
        
        [Fc_femur,Vc_femur]=ggremesh(Fc,Vc,optionStruct3);
        Ebs_femur=patchBoundary(Fc_femur,Vc_femur);
        indBs_femur=edgeListToCurve(Ebs_femur);
        indBs_femur=indBs_femur(1:end-1);
        
        %%
        z=0.5*(mean(Vc_femur(indBs_femur,3))+mean(Vc_femur(indBs_femur,3)));
        
        Vc_muscle(indBs_muscle,3)=z;%muscle
        Vc_femur(indBs_femur,3)=z;%femur
        
        %%
        pointSpacing=mean(patchEdgeLengths(Fc_muscle,Vc_muscle));
        [Ftt,Vtt]=regionTriMesh3D({Vc_muscle(indBs_muscle,:),Vc_femur(indBs_femur,:)},pointSpacing,0,'linear');
        
        %%
        cFigure;
        gpatch(Fc_muscle,Vc_muscle,'gw','k',1);
        gpatch(Ff,Vf,'rw','k',1);
        gpatch(Ftt,Vtt,'bw','k',0.6);
        colormap gjet;
        axisGeom; camlight headlight;
        gdrawnow;

        %%
        FT_amp{3}=Fc_muscle;
        VT_amp{3}=Vc_muscle;
        CT_amp{3}=3*ones(size(Fc_muscle,1),1);
        
        FT_amp{4}=Ftt;
        VT_amp{4}=Vtt;
        CT_amp{4}=4*ones(size(Ftt,1),1);
        
        %%
        cFigure; hold on;
        gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
        % patchNormPlot(FT_amp,VT_amp);
        axisGeom; camlight headlight;
        colormap gjet; icolorbar;
        gdrawnow;  
       
        %Cut the top skin surface
        Fs=FT_amp{2};
        Vs=VT_amp{2};
        
        %Use triSurfSlice to process cut
        snapTolerance=mean(patchEdgeLengths(Fs,Vs))/100;
        [Fc,Vc,~,logicSide]=triSurfSlice(Fs,Vs,[],P,-nCut,snapTolerance);
        [Fc,Vc]=patchCleanUnused(Fc(~logicSide,:),Vc); %Clean-up mesh
        
        % Remesh using ggremesh
        optionStruct3.pointSpacing=mean(patchEdgeLengths(Fs,Vs));
        optionStruct3.disp_on=0; % Turn off command window text display
        optionStruct3.pre.max_hole_area=100; %Max hole area for pre-processing step
        optionStruct3.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
        [Fc,Vc]=ggremesh(Fc,Vc,optionStruct3);
        
        Ebs=patchBoundary(Fc,Vc);
        indBs=edgeListToCurve(Ebs);
        indBs=indBs(1:end-1);
        
        Fc_skin=Fc;
        Vc_skin=Vc;
        indBs_skin=indBs;
        
        %Get the surface of the muscle
        Fc_muscle=FT_amp{3};
        Vc_muscle=VT_amp{3};

        Ebm=patchBoundary(Fc_muscle,Vc_muscle);
        indBs_muscle=edgeListToCurve(Ebm);
        indBs_muscle=indBs_muscle(1:end-1);

        %%
        z=0.5*(mean(Vc_muscle(indBs_muscle,3))+mean(Vc_muscle(indBs_muscle,3)));
        
        Vc_skin(indBs_skin,3)=z;%skin
        Vc_muscle(indBs_muscle,3)=z;%muscle
        
        %%
        pointSpacing=mean(patchEdgeLengths(Fc_skin,Vc_skin));
        [Fsm,Vsm]=regionTriMesh3D({Vc_skin(indBs_skin,:),Vc_muscle(indBs_muscle,:)},pointSpacing,0,'linear');
        
        %%
        cFigure;hold on;
        gpatch(Fc_skin,Vc_skin,'gw','k',1);
        gpatch(Fc_muscle,Vc_muscle,'rw','k',1);
        gpatch(Fsm,Vsm,'bw','k',0.6);
        plotV(Vc_skin(indBs_skin,:),'r-','LineWidth',4);
        plotV(Vc_muscle(indBs_muscle,:),'g-','LineWidth',4);

        colormap gjet;
        axisGeom; camlight headlight;
        gdrawnow;

        
        %%
        FT_amp{2}=fliplr(Fc_skin);
        VT_amp{2}=Vc_skin;
        CT_amp{2}=2*ones(size(Fc_skin,1),1);
        
        FT_amp{5}=Fsm;
        VT_amp{5}=Vsm;
        CT_amp{5}=5*ones(size(Fsm,1),1);
        
        %%
        cFigure; hold on;
        gpatch(FT_amp,VT_amp,CT_amp,'none',0.5);
        patchNormPlot(FT_amp,VT_amp);
        axisGeom; camlight headlight;
        colormap gjet; icolorbar;
        gdrawnow;  

%         stlwrite(triangulation(VT_amp{1},FT_amp{1}),fullfile(saveFolder_stl,'femur_tf.stl'),'binary');
%         stlwrite(triangulation(VT_amp{2},FT_amp{2}),fullfile(saveFolder_stl,'skin_tf.stl'),'binary');
%         stlwrite(triangulation(VT_amp{3},FT_amp{3}),fullfile(saveFolder_stl,'muscle_tf.stl'),'binary');
%         stlwrite(triangulation(VT_amp{4},FT_amp{4}),fullfile(saveFolder_stl,'muscle_lid_tf.stl'),'binary');
%         stlwrite(triangulation(VT_amp{5},FT_amp{5}),fullfile(saveFolder_stl,'skin_lid_tf.stl'),'binary'); 
        stlStruct.solidNames={'Femur'};
        stlStruct.solidVertices={VT_amp{1}};
        stlStruct.solidFaces={FT_amp{1}};
        stlStruct.solidNormals={[]};
        fileName=fullfile(saveFolder_stl,'femur_tf.stl');
        export_STL_txt(fileName,stlStruct);

        stlStruct.solidNames={'Skin'};
        stlStruct.solidVertices={VT_amp{2}};
        stlStruct.solidFaces={FT_amp{2}};
        stlStruct.solidNormals={[]};
        fileName=fullfile(saveFolder_stl,'skin_tf.stl');
        export_STL_txt(fileName,stlStruct);

        stlStruct.solidNames={'Muscle'};
        stlStruct.solidVertices={VT_amp{3}};
        stlStruct.solidFaces={FT_amp{3}};
        stlStruct.solidNormals={[]};
        fileName=fullfile(saveFolder_stl,'muscle_tf.stl');
        export_STL_txt(fileName,stlStruct);

        stlStruct.solidNames={'Muscle_lid'};
        stlStruct.solidVertices={VT_amp{4}};
        stlStruct.solidFaces={FT_amp{4}};
        stlStruct.solidNormals={[]};
        fileName=fullfile(saveFolder_stl,'muscle_lid_tf.stl');
        export_STL_txt(fileName,stlStruct);
        
        stlStruct.solidNames={'Skin_lid'};
        stlStruct.solidVertices={VT_amp{5}};
        stlStruct.solidFaces={FT_amp{5}};
        stlStruct.solidNormals={[]};
        fileName=fullfile(saveFolder_stl,'skin_lid_tf.stl');
        export_STL_txt(fileName,stlStruct);
        
        
end

%% Save model
if saveOn==1
    saveName_mat=fullfile(saveFolder,[saveNameGeom,'.mat']);
    save(saveName_mat,'FT_amp','VT_amp','CT_amp');
end

%%

function [Fsb,Vsb,indBs,Nd,P1,P2,P3,XB,YB,ZB]=smoothBezierClose(Fs,Vs,Ebs,P_end,numPoints)

indBs=edgeListToCurve(Ebs);
indBs=indBs(1:end-1); 

[~,~,Ns]=patchNormal(Fs,Vs);

[~,~,NEb]=edgeVec(Ebs,Vs);
NEb=vecnormalize(NEb);
Nd=vecnormalize(cross(NEb,Ns));

cPar.n=25; 
Nd=patchSmooth(Ebs,Nd,[],cPar);

%Bezier point set 1 
P1=Vs(indBs,:); 
P1_mid=mean(P1,1);

distanceEnd=sqrt(sum((P_end-P1_mid).^2,2));

%Bezier point set 2
f=1./abs(Nd(indBs,3)); %Extend factor for direction vectors
P2=P1+distanceEnd/2*f(:,ones(1,3)).*Nd(indBs,:); %Position half-way

%Bezier point set 3
P3=P2;
P3(:,3)=P_end(:,3); %Shift to bottom
P3=P3-mean(P3,1); %Shift on own centre
[t,r]=cart2pol(P3(:,1),P3(:,2)); %Polar coordinates
[P3(:,1),P3(:,2)]=pol2cart(t,distanceEnd/2*ones(size(r))); %Force constant radius
P3=P3+P_end; %Place centered on end

%Bezier point set 4
P4=P_end;

XB=zeros(numPoints,size(P1,1));
YB=zeros(numPoints,size(P1,1));
ZB=zeros(numPoints,size(P1,1));
for q=1:1:size(P1,1)
    p=[P1(q,:); P2(q,:); P3(q,:); P4]; %Control points    
    V_bezier=bezierCurve(p,numPoints); %Compute bezier curve
    XB(:,q)=V_bezier(:,1);
    YB(:,q)=V_bezier(:,2);
    ZB(:,q)=V_bezier(:,3);
end

[Fsb,Vsb]=meshToPatch(XB,YB,ZB,1);
[Fsb,Vsb]=quad2tri(Fsb,Vsb);
[Fsb,Vsb]=mergeVertices(Fsb,Vsb);
[Fsb,Vsb]=patchCleanUnused(Fsb,Vsb);

end

%%

function [F,V]=meshToPatch(X,Y,Z,closeSection)

%Create quad patch data
[F,V] = surf2patch(X,Y,Z);

%Close section if required
if closeSection==1
    I=[(1:size(Z,1)-1)' (1:size(Z,1)-1)' (2:size(Z,1))' (2:size(Z,1))' ];
    J=[size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1)];
    F_sub=sub2ind(size(Z),I,J);
    F=[F;F_sub];    
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
