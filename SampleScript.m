%% Guide Star method of obtaining absolute depth-resolved optic axis orientation of PS-OCT data
% This is a demonstration of guide star method

% [1] Absolute depth resolved optic axis measurement with catheter-based
% polarization sensitive optical coherence tomography
%% 
clc
clear
% Load the tomogram data from a .mat file
load(fullfile('DATA','Tomogram.mat'));
% Compute the intensity image using the log10 of the sum of squared
% magnitudes of two channels from Tomogram.
Intensity = 10*(log10(abs(tom.ch1(:,1:2:end)).^2 + abs(tom.ch2(:,2:2:end)).^2)) ;
%% Find sheath positions for the next step
axialRange=1101:1800;
imagesc(Intensity(axialRange,:),[65,115])
disp('-ACTION REQUIRED: Find the inner and outer sheath positions in the first A-line...')
%% Sheath and sample surface segmentation
% This section segments the inner sheath, outer sheath, and sample surface
% by calling the getSheathAndSampleInterfacePosition function.
inputImage=Intensity(axialRange,:);
% Define the initial positions for the inner and outer sheath
innerSheath = 61;
outerSheath = 108;
[innerSheathPosition, outerSheathPosition, sampleSurfacePosition] = getSheathAndSampleInterfacePosition(inputImage, innerSheath, outerSheath);
disp('-Sheath positions and sample surface are identified')
%% LOAD SYSTEM COMPENSATION
% Load system compensation data
% % it is better to regenerate system compensation
systemComp = systemCompensation(fullfile('DATA','1SystemCompensation.mat'));
% % disp('-System compensation is loaded')
%% Reconstruction of retradance matrix M and spectrally binned stokes vectors
numberOfSamplesPerAline = 3264;  % Number of samples in each A-line of the tomogram
Nmapping = numberOfSamplesPerAline / 2; %mapping number

% Convert tomogram data to Stokes vectors and apply scaling factor
[S1,S2,scalingFactor] = convertFullTomToStokesBins(tom,Nmapping,5);

dopThreshhold=0.7;
pstruct = struct;
pstruct.fwx = 12;
pstruct.dopTh =dopThreshhold;
pstruct.systemCompensation = systemComp;
pstruct.axial = axialRange;
pstruct.scalingFactor = scalingFactor;
retardanceMeasurement = PSProcessGlobalSymmetricCumulative(S1,S2,pstruct);
disp('-RetardanceMeasurement and rotationMatrix "MM" is reconstructed')
%% Retardanc Matrix at the sheath and sample surface positions
% Extract rotation matrices at the positions of interest (inner sheath, outer sheath, and sample surface)
rotationMatrix = retardanceMeasurement.MM;
sheathInRotMatrix = rotationMatrixAtPosition(rotationMatrix,innerSheathPosition);
sheathOutRotMatrix = rotationMatrixAtPosition(rotationMatrix,outerSheathPosition);
surfRotMatrix = rotationMatrixAtPosition(rotationMatrix,sampleSurfacePosition);
disp('-Rotation matrices at the requested positions are extracted')
%% Finding correction terms
% Apply corrections to the rotation matrices using the guide star method and catheter modeling
fw = 100;  % Filter width
depTh = 0; % Depth threshold for filtering

% Compute the guide star corrections
[leftCorrection,rightCorrection] = GuidestarCorrection(cleanCathRetTrace(sheathInRotMatrix,fw,depTh),cleanCathRetTrace(surfRotMatrix,fw,depTh));
GuidestarRotationMatrix= pagemtimes(permute(leftCorrection,[1,2,4,3]),pagemtimes(rotationMatrix,permute(rightCorrection,[1,2,4,3])));  
% Compute the catheter modeling corrections
[leftCorrectionCath,rightCorrectionCath] = CatheterModellingCorrection(cleanCathRetTrace(sheathInRotMatrix,fw,depTh),cleanCathRetTrace(surfRotMatrix,fw,depTh));
CatheterModelRotationMatrix= pagemtimes(permute(leftCorrectionCath,[1,2,4,3]),pagemtimes(rotationMatrix,permute(rightCorrectionCath,[1,2,4,3]))); 
disp('-Rotation matrices are corrected')
%% SECTION TITLE
% Compute the depth-resolved retardance matrices using both correction methods
degreeOfPolarization = retardanceMeasurement.dop;
structure.dop=degreeOfPolarization;
structure.dopTh=dopThreshhold;
structure.surf=sampleSurfacePosition;

localRetardanceGuidestar = getDepthResolvedRet(GuidestarRotationMatrix,structure);
localRetardanceCatheterModel = getDepthResolvedRet(CatheterModelRotationMatrix,structure);
disp('-Depth resolved rotation matrices are obtained')
%% RETARDANCE AND OPTIC AXIS ORIENTATION
% Extract local optic axis orientations computed with both methods
GuidestarLocalOA = squeeze(atan2(localRetardanceGuidestar(2,:,:),localRetardanceGuidestar(1,:,:)));
CatheterModellingLocalOA = squeeze(atan2(localRetardanceCatheterModel(2,:,:),localRetardanceCatheterModel(1,:,:)));

Retardance = squeeze(sqrt(sum(localRetardanceGuidestar.^2,1)));
disp('-Local optic axis orientations are computed with both methods')
%% Visualization
ax1=subplot(2,2,1);
imagesc(inputImage,[65,115])
hold on
plot(innerSheathPosition,'r')
hold on
plot(outerSheathPosition,'b')
hold on
plot(sampleSurfacePosition,'w')
title("Intensity",'FontSize',14)
colormap(ax1,gray)
colorbar
 

ax2=subplot(2,2,2);
imagesc(degreeOfPolarization)
title("Degree of polarization",'FontSize',14)
colormap(ax2,gray)
colorbar


ax3=subplot(2,2,3);
cmap=colorcet('C6');
rgb = ind2rgb(uint8(255*mod(GuidestarLocalOA./2/pi+1/2,1)),cmap);
mask= retardanceMeasurement.dop<dopThreshhold;
rgb(repmat(mask,[1,1,3])>0) = 0.6;
imagesc(rgb)
colormap(ax3,cmap)
cb = colorbar;
cb.Ticks = [0   0.25  0.5   0.75   1];
cb.TickLabels = {  '-π/2','-π/4','0', 'π/4', 'π/2'};
title("Guide star, Optic axis",'FontSize',14)

ax4=subplot(2,2,4);
rgb = ind2rgb(uint8(255*mod(CatheterModellingLocalOA./2/pi+1/2,1)),cmap);
rgb(repmat(mask,[1,1,3])>0) = 0.6;


imagesc(rgb)
colormap(ax4,cmap)
cb = colorbar;
cb.Ticks = [0   0.25  0.5   0.75   1];
cb.TickLabels = {  '-π/2','-π/4','0', 'π/4', 'π/2'};
title("Catheter modelling, Optic axis",'FontSize',14)