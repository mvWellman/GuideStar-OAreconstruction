%% Guide Star method of obtaining absolute depth-resolved optic axis orientation of PS-OCT data
% This is a demonstration of guide star method

% [1] Absolute depth resolved optic axis measurement with catheter-based
% polarization sensitive optical coherence tomography
%% Load Intensity image and find sheath positions for the next step
% DESCRIPTIVE TEXT
clc
clear
load('Intensity.mat');
axialRange=1101:1800;
imagesc(Intensity(axialRange,:))
fprintf(" . \n . \n . \n -ACTION REQUIRED: Find the inner and outer sheath positions in the first A-line...")
%% Sheath and sample surface segmentation
% DESCRIPTIVE TEXT
Str.CatheterType='LowProfile';
Str.OuterSheath='Yes'; 
innerSheath = 61;
outerSheath = 108;
Int=Intensity(axialRange,:);
[InSheathPosition,OutSheathPosition,Surf] = CatheterSegmentation(Int,innerSheath,outerSheath,Str);
fprintf("\n -Sheath positions and sample surface are identified")
%% LOAD SYSTEM COMPENSATION
load('SystemCompensation.mat');
fprintf("\n -System compensation is loaded")
%% Reconstruction of retradance matrix M, using spectrally binned stokes vector
% DESCRIPTIVE TEXT
load('StokesVectors.mat')
dopTh=0.7;
pstruct.fwx = 12;
pstruct.dopTh =dopTh;
pstruct.systemCompensation = loc;
pstruct.axial = axialRange;
out = RetardanceMatrix(Stokes{1,1},Stokes{1,2},pstruct);
fprintf("\n -Retardance Matrix MM is reconstructed")
%% Retardanc Matrix at the sheath and sample surface positions
% DESCRIPTIVE TEXT
MM = out.MM;
sheathInR = RMatrixPosition(MM,InSheathPosition);
sheathOutR = RMatrixPosition(MM,OutSheathPosition);
surfR = RMatrixPosition(MM,Surf);
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
st.surf=Surf;
st.dzRange=axialRange;
retlocGuidestar = DepthResolvedRet(MMguidestar,st);
retlocCatheter = DepthResolvedRet(MMCatheter,st);
%% SECTION TITLE
% DESCRIPTIVE TEXT Some visualization

retGuidestar = squeeze(sqrt(sum(retlocGuidestar.^2,1)));
phiGuidestar = squeeze(atan2(retlocGuidestar(2,:,:),retlocGuidestar(1,:,:)));

retCatheter = squeeze(sqrt(sum(retlocCatheter.^2,1)));
phiCatheter = squeeze(atan2(retlocCatheter(2,:,:),retlocCatheter(1,:,:)));
%% SECTION TITLE
% DESCRIPTIVE TEXT

mask=Int>=75;
mask=imfill(mask,"holes");

ax1=subplot(2,2,1);
imagesc(Int,[65,115])
hold on
plot(InSheathPosition,'r')
hold on
plot(OutSheathPosition,'b')
hold on
plot(Surf,'w')
title("Intensity",'FontSize',14)
colormap(ax1,gray)

 

ax2=subplot(2,2,2);
imagesc(out.dop)
title("Depolarization",'FontSize',14)
colormap(ax2,"hot")
 

ax5=subplot(2,2,3);
imagesc(phiGuidestar.*mask)
colormap(ax5,"hsv")
 
title("Guide star, Optic axis",'FontSize',14)
% colormap(ax3,cmapD)
ax6=subplot(2,2,4);
imagesc(phiCatheter.*mask)
title("Catheter modelling, Optic axis",'FontSize',14)
colormap(ax6,"hsv")
 








