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
fontSize=40;
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
        %% Working with original file of the body surfaces wrapped around organs to cut the right side of the body preserving
        %% everithing below the cutHeight
        logicVertices=V(:,1)<0 & V(:,3)<=cutHeight+mean(patchEdgeLengths(F,V));
        logicFaces=any(logicVertices(F),2);
        logicFaces=triSurfLogicSharpFix(F,logicFaces);
        Fs=F(logicFaces,:);
        [Fs,Vs]=patchCleanUnused(Fs,V);
        
        cFigure; hold on;        
        gpatch(Fs,Vs,'w','none',0.25);
        axis off;
        axisGeom;
        camlight headlight;
        drawnow;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Cut again preserving only right leg
        optionStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Fs,optionStruct);
        [~,indMax]=max(groupSize);
        logicKeep=G==indMax;
        Fs=Fs(logicKeep,:);
        [Fs,Vs]=patchCleanUnused(Fs,Vs);
        
        cFigure; hold on;
        gpatch(Fs,Vs,'w','none',0.25);
        axis off;
        axisGeom;
        camlight headlight;
        drawnow;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Slice the surface of the right leg at cutHeight
        %% to detect all edges to form the list of curves within the slice
        snapTolerance=mean(patchEdgeLengths(Fs,Vs))/100;
        n=vecnormalize([0 0 1]); %Normal direction to plane
        P=[0 0 cutHeight]; %Point on plane
        [Fc,Vc,~,logicSide,Ec]=triSurfSlice(Fs,Vs,[],P,n,snapTolerance);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Ec,groupStruct);
        
        for i=1:size(groupSize,2)
            logicKeep=G==i;
            Ec_keep=Ec(logicKeep,:);
            indCurveNow=edgeListToCurve(Ec_keep);
            VTcc{i}=Vc(indCurveNow,:);
              
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Find in the list the curve with the min value along X axis with asociated index 
        %% which corresponds to the edge of the skin surface
        for i=1:size(groupSize,2)
            [val(i),idx(i)] = min([VTcc{i}(:,1)]);   
        end
        [~,indKeep]=min(val(:));
        logicKeep=G==indKeep;
        Ec_keep=Ec(logicKeep,:);
        indCurveNow=edgeListToCurve(Ec_keep);
        VTc=Vc(indCurveNow,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Forming the 2D convex hull around the curve
        [V_slice_curve]=convexhull_curve(VTc,200);
        [Fc,Vc]=patchCleanUnused(Fc(logicSide,:),Vc);
        
        cFigure; hold on;        
        gpatch(Fc,Vc,'w','none',0.25);
        plotV(V_slice_curve,'k-','LineWidth',4);
        axis off;
        axisGeom;
        camlight headlight;
        drawnow;

        %% Working with right leg to cut preserving
        %% everyting below the sex glands.
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
        gpatch(F,V,'w','none',0.25);
        gpatch(Fe,Ve,'rw','none',1);
        plotV(V_slice_curve,'k-','LineWidth',4);
        axis off;
        axisGeom;
        camlight headlight;
        drawnow;
        
        %% Forming the point cloud out of the processed cut right leg a
        %% and the edge of the skin (curve) formed at cutHeight
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
        hAxis(1)=subplot(1,2,1); hold on;
        %gpatch(Fs,Vs,'w','none',0.25);
        gpatch(Fe,Ve,'w','none',1);
        plotV(V_slice_curve,'k-','LineWidth',4);

        text(min(V(:,1))+50,min(V(:,2))+99,'A','FontSize',fontSize);
        axis off;
        axisGeom;
        camlight headlight;
        
        hAxis(2)=subplot(1,2,2); hold on;
        gpatch(Fn,Vn,'w','k',1,lineWidth);
        text(min(V(:,1))+40,min(V(:,2))+110,'B','FontSize',fontSize);
        axis off;
        axisGeom;
        pos = get( hAxis(2), 'Position' );
        pos(1)=pos(1)-0.3;
        set( hAxis(2), 'Position', pos ) ;
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Define the muscle.         
        %% Working with original file of the body surfaces wrapped around organs to cut the right side of the body preserving
        %% everything below the cutHeight
        logicVertices=Vr(:,1)<1 & Vr(:,3)<=cutHeight+mean(patchEdgeLengths(Fr,Vr));
        logicFaces=any(logicVertices(Fr),2);
        logicFaces=triSurfLogicSharpFix(Fr,logicFaces);
        Fs=Fr(logicFaces,:);
        [Fs,Vs]=patchCleanUnused(Fs,Vr);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Cut preserving everything from the right leg
        optionStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Fs,optionStruct);
        [~,indMax]=max(groupSize);
        logicKeep=G==indMax;
        Fs=Fs(logicKeep,:);
        [Fe,Ve]=patchCleanUnused(Fs,Vs);
        
        cFigure; hold on;  
        gpatch(Fr,Vr,'w','none',0.25);
        gpatch(Fe,Ve,'r','none',1);
        axis off;
        axisGeom;
        camlight headlight;
        drawnow;

        snapTolerance=mean(patchEdgeLengths(Fe,Ve))/100;
        numLoftGuideSlicesFat=60;
        n=vecnormalize([0 0 1]);
        numPoints=200;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Set up lower and upper limits for the muscle height
        %% Find these limits on the skin surface
        [~,indMax]=max(Ve(:,3));
        [~,indMin]=min(Ve(:,3));
        P=Ve(indMax,:);
        PA=Ve(indMin,:);
        P(:,3)=P(:,3)-122.5;%Upper limit 1
        P0=P;
        PA(:,3)=PA(:,3)+110;%Lower limit
        zz=linspacen(P(:,3),PA(:,3),numLoftGuideSlicesFat);%number of slices within the range between limits
        e=[(1:numPoints)' [(2:numPoints)';1]];
        VC=cell(numLoftGuideSlicesFat,1);
        
        q=1;
        P=Ve(indMax,:);
        P(:,3)=P(:,3)-10;%Upper limit
        [Ff,Vf,~,logicSide,Eb]=triSurfSlice(Fe,Ve,[],P,n,snapTolerance);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct); 
        for i=1:size(groupSize,2)
            logicKeep=G==i;
            Eb_keep=Eb(logicKeep,:);
            indCurveNow=edgeListToCurve(Eb_keep);
            VT=Vf(indCurveNow,:);
            [V_slice_curve]=convexhull_curve(VT,numPoints);
            Vcc{i}=V_slice_curve;
        end

        cFigure;hold on;
        gpatch(Ff(logicSide,:),Vf,'w','none',1);
        gpatch(Ff(~logicSide,:),Vf,'w','none',0.25);
        gpatch(Eb,Vf,'none','b',1,3); 
        for q=1:1:size(Vcc,2)
            Vpp=Vcc{q};
            plotV(Vpp([1:size(Vpp,1) 1],:),'k-','LineWidth',4);
        end
        axisGeom; axis off;
        axis manual; camlight headligth;
        set(gca,'FontSize',25);
        gdrawnow;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% first_min_value along X axis which corresponds to the edge of the skin surface
        for i=1:size(groupSize,2)
            [val(i),idx(i)] = min([Vcc{i}(:,1)]);   
        end
        [~,idx_a]=min(val(:));
        index=linspace(1,size(groupSize,2),size(groupSize,2));
        index1 = index(find(index~=idx_a));
        %% second_min_value along X axis which corresponds to the edge of the muscle surface
        for i=1:size(index1,2)
            [val(i),idx(i)] = min([Vcc{index1(i)}(:,1)]);   
        end
        [~,idx_a]=min(val(:));
        V_slice_curve=Vcc{idx_a};
        
        %% Reference upper curve corresponds to the beginning of the muscle layer
        VC{q}=V_slice_curve;
        FL=[];
        VL=V_slice_curve;
        [~,indStart]=max(VL(:,1));
        if indStart>1
            VL=[VL(indStart:end,:); VL(1:indStart-1,:)];
        end
        
        VL0=VL;
        VL(:,1)=0.7.*VL(:,1);
        Xshift=mean(VL,1)-mean(VL0,1);
        VL(:,1)=VL(:,1)-Xshift(1)-30; 
        
        VL0=VL;
        [Fh,Vh]=regionTriMesh3D({VL},[],0,'linear');
        Fh=fliplr(Fh);
        
        
        cFigure;hold on;
        gpatch(Ff(logicSide,:),Vf,'w','none',1);
        gpatch(Ff(~logicSide,:),Vf,'w','none',0.25);
        gpatch(Fh,Vh,'w','k',0.25);
        gpatch(Eb,Vf,'none','b',1,3);
        plotV(VL,'k-','LineWidth',4);
        axisGeom; %axis off;
        axis manual; camlight headligth;
        set(gca,'FontSize',25);
        gdrawnow;

        VC{1}=VL;
        %% Loop over number of slices along Z axis
        for i=2:1:numLoftGuideSlicesFat
            P(:,3)=zz(i);
            [Ff,Vf,~,logicSide,Eb]=triSurfSlice(Fe,Ve,[],P,n,snapTolerance);
            groupStruct.outputType='label';
            [G,~,groupSize]=tesgroup(Eb,groupStruct);
            
            %
            %             for q=1:1:max(size(groupSize))
            %             logicKeep=G==q;
            %             Eb_keep=Eb(logicKeep,:);
            %             indCurveNow=edgeListToCurve(Eb_keep);
            %             VT=Vf(indCurveNow,:);
            %             [V_slice_curve]=convexhull_curve(VT,numPoints);
            %             Vcc{i}=V_slice_curve;
            %             end
            
            VT={};
            for q=1:1:max(size(groupSize))
                logicKeep=G==q;
                Eb_keep=Eb(logicKeep,:);
                indCurveNow=edgeListToCurve(Eb_keep);
                VT{q}=Vf(indCurveNow',:);
                Xmin(q)=min(VT{q}(:,1));
            end
            
            [Xmin_value,indOut]=min(Xmin);
            count=1;
            VTT={};
            for q=1:1:max(size(groupSize))
                if q~=indOut
                    indKeep=q;
                    logicKeep=G==indKeep;
                    Eb_keep=Eb(logicKeep,:);
                    indCurveNow=edgeListToCurve(Eb_keep);
                    VTT{count}=Vf(indCurveNow,:);
                    count=count+1;
                end
            end
            
            V_points_keep=[];
            for q=1:1:max(size(VTT))
                Vpp=VTT{q};
                V_points_keep=[V_points_keep;Vpp];
            end
            
            [V_slice_curve]=convexhull_curve(V_points_keep,numPoints);
            f=[e+(i-2)*numPoints fliplr(e)+(i-1)*numPoints];
            
            
            FL=[FL;f];
            VL=[VL;V_slice_curve;];
            VC{i}=V_slice_curve;
            
%             cFigure;hold on;
%             gpatch(Ff(logicSide,:),Vf,'w','none',1);
%             gpatch(Ff(~logicSide,:),Vf,'w','none',0.25);
%             gpatch(Eb,Vf,'none','b',1,3);
%             plotV(VL,'r-','LineWidth',4);
%             axisGeom; axis off;axis manual; camlight headligth;
%             set(gca,'FontSize',25);
%             gdrawnow;

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
        gpatch(Ff(logicSide,:),Vf,'w','none',1);
        gpatch(Ff(~logicSide,:),Vf,'w','none',0.25);
        for q=1:1:size(VC,1)
            Vpp=VC{q};
           
            plotV(Vpp([1:size(Vpp,1) 1],:),'k-','LineWidth',4);
        end
        axisGeom;axis off;
        camlight headligth;
        set(gca,'FontSize',25);
        
        subplot(1,2,2); hold on;
        gpatch(FL,VL,'bw','k',1);        
        axisGeom; axis off;
        camlight headligth;
        set(gca,'FontSize',25);
        gdrawnow;

        % Join and merge surfaces
        [F_muscles,V_muscles,C_muscles]=joinElementSets({Fg,FL,Fh},{Vg,VL,Vh});
        [F_muscles,V_muscles]=mergeVertices(F_muscles,V_muscles);
        
        cFigure; hold on;
        %title('Completed skin+fat surface');
        gpatch(Ff,Vf,'w','none',0.5);
        gpatch(F_muscles,V_muscles,C_muscles,'none',1);

        for q=1:1:size(VC,1)
            Vpp=VC{q};
            plotV(Vpp([1:size(Vpp,1) 1],:),'k-','LineWidth',4);
        end
        
        axisGeom; camlight headlight;
        colormap(gjet(250)); icolorbar;
        gdrawnow;


        cPar1.Method='HC'; %Smooth method
        cPar1.n=5; %Number of iterations
        [V_muscles_smoothed]=patchSmooth(F_muscles,V_muscles,[],cPar1);
    
        %% Remesh
        optionStructRemesh.pointSpacing=4; %Set desired point spacing
        [Fn,Vn]=ggremesh(F_muscles,V_muscles_smoothed,optionStructRemesh);
    
        cFigure;hold on;
        gpatch(Fe,Ve,'w','none',0.1);
        gpatch(Fn,Vn,'r','k',1);
        axisGeom; axis off; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        gdrawnow;
        
        cFigure;        
        a1=subplot(1,2,1); hold on;
        gpatch(Ff(logicSide,:),Vf,'w','none',1);
        gpatch(Ff(~logicSide,:),Vf,'w','none',0.25);
        for q=1:1:size(VC,1)
            Vpp=VC{q};
            plotV(Vpp([1:size(Vpp,1) 1],:),'k-','LineWidth',4);
        end
        
        text(min(Vf(:,1))-100,min(Vf(:,2))+130,'A','FontSize',fontSize);
        axisGeom;axis off;
        axis manual; 
        camlight headligth;
        set(gca,'FontSize',25);

        a2=subplot(1,2,2); hold on;
        %gpatch(Fe,Ve,'w','none',0.1);
        gpatch(Fn,Vn,'w','k',1);
        text(min(Vn(:,1))+50,min(Vn(:,2))+250,'B','FontSize',fontSize);
        axisGeom; axis off; 
        axis manual; 
        pos = get( a2, 'Position' );
        pos(1)=pos(1)-0.3;
        set( a2, 'Position', pos ) ;
        camlight headligth;
        set(gca,'FontSize',25);   
        gdrawnow;
        

        %% Store model in structure
        modelNew.source=fileName_FMA_refined;
        modelNew.faces=Fn;
        modelNew.vertices=Vn;
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
