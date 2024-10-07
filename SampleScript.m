%% Guide Star method of obtaining absolute depth-resolved optic axis orientation of PS-OCT data
% This is a demonstration of guide star method

% [1] Absolute depth resolved optic axis measurement with catheter-based
% polarization sensitive optical coherence tomography
%% Load Intensity image and find sheath positions for the next step
% DESCRIPTIVE TEXT
clc
clear
addpath(genpath('DATA'));
load(fullfile('DATA','Intensity.mat'))
axialRange=1101:1800;
imagesc(Intensity(axialRange,:))
fprintf(" . \n . \n . \n -ACTION REQUIRED: Find the inner and outer sheath positions in the first A-line...")
%% Sheath and sample surface segmentation
% DESCRIPTIVE TEXT
addpath(genpath('Codes'));
Int=Intensity(axialRange,:);
innerSheath = 61;
outerSheath = 108;
TilenumberOfALines = 32;
[InnerSheathPosition, OuterSheathPosition, SampleSurface] = CatheterSegmentationClean(Int,innerSheath,outerSheath, TilenumberOfALines);
fprintf("\n -Sheath positions and sample surface are identified")
% InnerSheathPosition= InnerSheathPosition';
% OuterSheathPosition= OuterSheathPosition';

%% LOAD SYSTEM COMPENSATION
load('SystemCompensation.mat');
fprintf("\n -System compensation is loaded")
%% Reconstruction of retradance matrix M, using spectrally binned stokes vector
% DESCRIPTIVE TEXT

% load('StokesVectors.mat')
load('S1.mat')
load('S2.mat')
dopTh=0.7;
pstruct.fwx = 12;
pstruct.dopTh =dopTh;
pstruct.systemCompensation = loc;
pstruct.axial = axialRange;
out = RetardanceMatrix(S1,S2,pstruct);
fprintf("\n -Retardance Matrix MM is reconstructed")
%% Retardanc Matrix at the sheath and sample surface positions
% DESCRIPTIVE TEXT
MM = out.MM;
sheathInR = RMatrixPosition(MM,InnerSheathPosition);
sheathOutR = RMatrixPosition(MM,OuterSheathPosition);
surfR = RMatrixPosition(MM,SampleSurface);
fprintf("\n -Retardance matrices at the requested positions are reconstructed")
%% Finding correction terms
% DESCRIPTIVE TEXT
fw = 100;
depTh = 0;
[LHS,RHS,outStruct] = findRotAndSheathCorrectionGuidestar(cleanCathRetTrace(sheathInR(:,:,:),fw,depTh),cleanCathRetTrace(surfR(:,:,:),fw,depTh));
[LHS1,RHS1,outStruct1] = findRotAndSheathCorrectionApparent(cleanCathRetTrace(sheathInR(:,:,:),fw,depTh),cleanCathRetTrace(surfR(:,:,:),fw,depTh));
 
%% SECTION TITLE
% DESCRIPTIVE TEXT

MMguidestar= pagemtimes(permute(LHS,[1,2,4,3]),pagemtimes(MM,permute(RHS,[1,2,4,3])));  
MMCatheter= pagemtimes(permute(LHS1,[1,2,4,3]),pagemtimes(MM,permute(RHS1,[1,2,4,3]))); 
st.dop=out.dop;
st.dopTh=dopTh;
st.surf=SampleSurface;
st.dzRange=axialRange;
retlocGuidestar = DepthResolvedRet(MMguidestar,st);
retlocCatheter = DepthResolvedRet(MMCatheter,st);
%% SECTION TITLE
% DESCRIPTIVE TEXT Some visualization

retGuidestar = squeeze(sqrt(sum(retlocGuidestar.^2,1)));
phiGuidestar = squeeze(atan2(retlocGuidestar(2,:,:),retlocGuidestar(1,:,:)));

retCatheter = squeeze(sqrt(sum(retlocCatheter.^2,1)));
phiCatheter = squeeze(atan2(retlocCatheter(2,:,:),retlocCatheter(1,:,:)));


% retVecGuidestar = signal2isolum(phiGuidestar,retGuidestar,[-pi,pi],[0,.1],'C2');
% retVecCatheter = signal2isolum(phiCatheter,retCatheter,[-pi,pi],[0,.1],'C2');
%% SECTION TITLE
% DESCRIPTIVE TEXT

mask=Int>=75;
mask=imfill(mask,"holes");

ax1=subplot(2,2,1);
imagesc(Int,[65,115])
hold on
plot(InnerSheathPosition,'r')
hold on
plot(OuterSheathPosition,'b')
hold on
plot(SampleSurface,'w')
title("Intensity",'FontSize',14)
colormap(ax1,gray)
colorbar
 

ax2=subplot(2,2,2);
imagesc(out.dop)
title("Depolarization",'FontSize',14)
colormap(ax2,"hot")
 colorbar

 PHI1=phiGuidestar.*mask;
ax5=subplot(2,2,3);
imagesc(PHI1)
colormap(ax5,"hsv")
 colorbar
title("Guide star, Optic axis",'FontSize',14)
% colormap(ax3,cmapD)
PHI2=phiCatheter.*mask;
ax6=subplot(2,2,4);
imagesc(PHI2)
title("Catheter modelling, Optic axis",'FontSize',14)
colormap(ax6,"hsv")
colorbar
 








