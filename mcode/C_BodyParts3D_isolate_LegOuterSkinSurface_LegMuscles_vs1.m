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
fileName_FMA_refined='skin';
saveNameMesh='right_leg_skin';
saveNameMesh1='right_leg_muscles';

pointSpacing=4; 
cutHeight=810;
saveOn=1;


Select_layer='outer_muscles_surface';%'outer_skin_surface';'outer_muscles_surface'
switch Select_layer
    case 'outer_skin_surface'
        %% Original
        fileName_mat=fullfile(loadFolder,[fileName_FMA,'.mat']);
        model=load(fileName_mat);
        F=model.faces;
        V=model.vertices;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Working with original file of the skin to cut the right leg preserving
        %everithing below the sex glands.
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
        [Fc,Vc,~,logicSide,Ec]=triSurfSlice(Fs,Vs,[],P,n,snapTolerance);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Ec,groupStruct);
        [~,indKeep]=max(groupSize);
        logicKeep=G==indKeep;
        Ec_keep=Ec(logicKeep,:);
        indCurveNow=edgeListToCurve(Ec_keep);
        VTc=Vc(indCurveNow,:);
        [V_slice_curve]=convexhull_curve(VTc,200);
        
        [Fc,Vc]=patchCleanUnused(Fc(logicSide,:),Vc);
        %Set-up orientation and location of cutting plane
        snapTolerance=mean(patchEdgeLengths(Fc,Vc))/100; %Tolerance for surface slicing
        n=vecnormalize([0 0 1]); %Normal direction to plane
        Q1=euler2DCM([0.4*pi 0 0.8*pi]);
        n=n*Q1;
        P_cut=[0 0 0]+n*306; %Point on plane
        %Slicing surface
        [Fd,Vd,Cd,logicSide,Ed]=triSurfSlice(Fc,Vc,[],P_cut,n);
        %Compose isolated cut geometry and boundary curves
        [Fe,Ve]=patchCleanUnused(Fd(logicSide==1,:),Vd);

            
        cFigure; hold on;
        title('The femur surface');
        
        gpatch(F,V,'w','none',0.25);
        gpatch(Fe,Ve,'rw','none',1);
        plotV(V_slice_curve,'k-','LineWidth',4);
        
        axisGeom;
        camlight headlight;
        drawnow;
        
        Ve=[Ve;V_slice_curve;];

        %% Construct alpha shape
        % This deteriorates surface quality (bridges concave regions) but is a
        % termporary "quick-and-dirty" work-around to obtain outer skin surface
        % only.
        shp = alphaShape(Ve,max(patchEdgeLengths(Fe,Ve)),'HoleThreshold',500,'RegionThreshold',100);
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
    case 'outer_muscles_surface'
        
        %% Refined
        fileName_mat=fullfile(saveFolder,[fileName_FMA_refined,'.mat']);
        model=load(fileName_mat);
        Fr=model.faces;
        Vr=model.vertices;


        %% Define the muscle layer. This involves working with the remesh mesh of the
        %original skin mesh. Right now it's the manual process.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Cut the right leg preserving everithing below the sex glands.
        logicVertices=Vr(:,1)<1 & Vr(:,3)<=cutHeight+mean(patchEdgeLengths(Fr,Vr));
        logicFaces=any(logicVertices(Fr),2);
        logicFaces=triSurfLogicSharpFix(Fr,logicFaces);
        Fs=Fr(logicFaces,:);
        [Fs,Vs]=patchCleanUnused(Fs,Vr);
        
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
        
        %Set-up orientation and location of cutting plane
        snapTolerance=mean(patchEdgeLengths(Fc,Vc))/100; %Tolerance for surface slicing
        n=vecnormalize([0 0 1]); %Normal direction to plane
        P_cut=[0 0 0]+n*cutHeight; %Point on plane
        
        %Slicing surface
        [Fd,Vd,Cd,logicSide,~]=triSurfSlice(Fc,Vc,[],P_cut,n);
        %Compose isolated cut geometry and boundary curves
        [Fe,Ve]=patchCleanUnused(Fd(logicSide==1,:),Vd);
        
        cFigure; hold on;
        title('The femur surface');
        
        gpatch(Fr,Vr,'w','none',0.25);
        gpatch(Fe,Ve,'rw','none',1);
        
        axisGeom;
        camlight headlight;
        drawnow;

        snapTolerance=mean(patchEdgeLengths(Fe,Ve))/100;
        numLoftGuideSlicesFat=41;
        n=vecnormalize([0 0 1]);
        numPoints=200;
        %Find lowest and highest point on the skin
        [~,indMax]=max(Ve(:,3));
        [~,indMin]=min(Ve(:,3));
        P=Ve(indMax,:);
        PA=Ve(indMin,:);
        P(:,3)=P(:,3)-110;
        P0=P;
        PA(:,3)=PA(:,3)+110;
        zz=linspacen(P(:,3),PA(:,3),numLoftGuideSlicesFat);
        e=[(1:numPoints)' [(2:numPoints)';1]];
        VC=cell(numLoftGuideSlicesFat-4,1);
        
        q=1;
        P=Ve(indMax,:);
        P(:,3)=P(:,3)-6;
        [Ff,Vf,~,logicSide,Eb]=triSurfSlice(Fe,Ve,[],P,n,snapTolerance);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        [~,indKeep]=max(groupSize);
        indKeep=2;
        logicKeep=G==indKeep;
        Eb_keep=Eb(logicKeep,:);
        indCurveNow=edgeListToCurve(Eb_keep);
        VT=Vf(indCurveNow,:);
        [V_slice_curve]=convexhull_curve(VT,numPoints);
        VC{q}=V_slice_curve;
        FL=[];
        VL=V_slice_curve;
        
          
        [Fh,Vh]=regionTriMesh3D({VL},[],0,'linear');
        Fh=fliplr(Fh);
        
        %%
        for q=2:1:numLoftGuideSlicesFat-4
            P(:,3)=zz(q);
            if (q==20)
                P(:,3)=zz(q)+6;
            end
            [Ff,Vf,~,logicSide,Eb]=triSurfSlice(Fe,Ve,[],P,n,snapTolerance);
            groupStruct.outputType='label';
            [G,~,groupSize]=tesgroup(Eb,groupStruct);
            
            if (q==2 | q==3 | q==5 | q==6 |q==7 | q==13 | q==15 | q==16 | q==18 | q==19| q==20)
                [~,indKeep]=max(groupSize);
            end
            if (q==4 | q==8 | q==9 | q==10 | q==11 | q==12 | q==14 | q==17 | q==22 | q==24 | q==25 | q==27 ...
                | q==28 | q==30 | q==31 | q==33 | q==34 | q==36 | q==37)
                indKeep=2;
            end
            if (q==23)
                indKeep=4;
            end
            if (q==26 | q==29 | q==35)
                indKeep=3;
            end
            if (q==32)
                indKeep=1;
            end
            
            
            if (q==21)
                indKeep1=2;
                indKeep2=3;
                indKeep3=4;
                logicKeep1=G==indKeep1;
                logicKeep2=G==indKeep2;
                logicKeep3=G==indKeep3;
                Eb_keep1=Eb(logicKeep1,:);
                Eb_keep2=Eb(logicKeep2,:);
                Eb_keep3=Eb(logicKeep3,:);
                indCurveNow1=edgeListToCurve(Eb_keep1);
                indCurveNow2=edgeListToCurve(Eb_keep2);
                indCurveNow3=edgeListToCurve(Eb_keep3);
                VT=[Vf(indCurveNow1,:);Vf(indCurveNow2,:);Vf(indCurveNow3,:);];
            end
            
            if (q~=21)
            logicKeep=G==indKeep;
            Eb_keep=Eb(logicKeep,:);
            indCurveNow=edgeListToCurve(Eb_keep);
            VT=Vf(indCurveNow,:);
            end
            
            [V_slice_curve]=convexhull_curve(VT,numPoints);
            f=[e+(q-2)*numPoints fliplr(e)+(q-1)*numPoints];
            
            
            FL=[FL;f];
            VL=[VL;V_slice_curve];
            VC{q}=V_slice_curve;
            
            cFigure;hold on;
            gpatch(Ff(logicSide,:),Vf,'w','none',1);
            gpatch(Ff(~logicSide,:),Vf,'w','none',0.25);
            gpatch(Eb,Vf,'none','b',1,3);
            plotV(VL,'r-','LineWidth',4);
            axisGeom; axis manual; camlight headligth;
            set(gca,'FontSize',25);
            gdrawnow;

        end

        [Fg,Vg]=regionTriMesh3D({V_slice_curve},[],0,'linear');
        ne=mean(patchNormal(Fg,Vg),1);
        if dot(ne,[0 0 1])>0
            Fg=fliplr(Fg);
        end
        
        [FL,VL]=subQuad(FL,VL,3,3);
        [FL,VL]=quad2tri(FL,VL,'a');
        FL=fliplr(FL);
        
        cFigure;
        subplot(1,2,1); hold on;
        gpatch(Fe,Ve,'rw','none',1);
        axisGeom;
        camlight headligth;
        set(gca,'FontSize',25);
        
        subplot(1,2,2); hold on;
        plotV(P0,'g.','MarkerSize',30);
        plotV(PA,'g.','MarkerSize',30);
        gpatch(Fg,Vg,'w','k',1);
        %patchNormPlot(Fg,Vg);
        
        gpatch(Fh,Vh,'w','k',1);
        %patchNormPlot(Fh,Vh);
        
        gpatch(FL,VL,'bw','k',1);
        %patchNormPlot(FL,VL);
        
        gpatch(Ff,Vf,'w','none',0.25);
        
        
        nn=size(VC);
        for q=1:1:nn(1)
            Vpp=VC{q};
            plotV(Vpp([1:size(Vpp,1) 1],:),'k-','LineWidth',4);
        end
        
        axisGeom; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        gdrawnow;
        
        % Join and merge surfaces
        [F_muscles,V_muscles,C_muscles]=joinElementSets({Fg,FL,Fh},{Vg,VL,Vh});
        [F_muscles,V_muscles]=mergeVertices(F_muscles,V_muscles);
        
        cFigure; hold on;
        title('Completed skin+fat surface');
        gpatch(Ff,Vf,'w','none',0.5);
        gpatch(F_muscles,V_muscles,C_muscles,'none',1);

        for q=1:1:nn(1)
            Vpp=VC{q};
            plotV(Vpp([1:size(Vpp,1) 1],:),'k-','LineWidth',4);
        end
        
        axisGeom; camlight headlight;
        colormap(gjet(250)); icolorbar;
        gdrawnow;


        cPar1.Method='HC'; %Smooth method
        cPar1.n=5; %Number of iterations
        [V_muscles_smoothed]=patchSmooth(F_muscles,V_muscles,[],cPar1);
    
        cFigure;hold on;
        gpatch(F_muscles,V_muscles_smoothed,'w','k',1);
        axisGeom; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        gdrawnow;

        %% Store model in structure
        modelNew.source=fileName_FMA_refined;
        modelNew.faces=F_muscles;
        modelNew.vertices=V_muscles;
        modelNew.pointSpacing=pointSpacing;
        modelNew.cutHeight=cutHeight;
        
        %% Save model
        if saveOn==1
            saveName_mat=fullfile(saveFolder,[saveNameMesh1,'.mat']);
            save(saveName_mat,'-struct','modelNew')
        end
end



function [V_slice_curve]=convexhull_curve(VT,numPoints)
P1=VT(:,[1 2]); %Point set flattened to 2D
DT=delaunayTriangulation(P1); %Delaunay triangulation of 2D set
VD=DT.Points; %Delaunay point set
VD(:,3)=min(VT(:,3)); %Set z-coord to minimum for now
indChull=DT.convexHull; %Ordered point list for conhex hull
indChull=indChull(1:end-1); %Trim away last (=start) point to avoid double
V_chull=VD(indChull,:); %Vertices for convex hull
D=max(pathLength(V_chull)); %Get length of sole curve for resampling
numResample=numPoints;%ceil(D./snapTolerance);
%V_slice_curve=evenlySampleCurve(V_chull,numResample,'pchip',1);
V_slice_curve=evenlySampleCurve(V_chull,numResample,0.01,1);


[~,indStart]=max(V_slice_curve(:,1));
    if indStart>1
        V_slice_curve=[V_slice_curve(indStart:end,:); V_slice_curve(1:indStart-1,:)];
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
