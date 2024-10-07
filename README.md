# GuideStar-OAreconstruction
Absolute depth resolved optic axis reconstruction of catheter based PS-OCT system.

## Table of Contents
1. Catheter-based PS-OCT
2. GuideStar processing overview
3. Example data description

## 2.  Catheter-based PS-OCT
Catheter-based PS-OCT is an optical imaging technique which allows for minimially invasive characterization of biological tissues. Optical coherence tomography (OCT) has been used in several clinical areas including opthalmology, cardiology, gastroenterology and neurosurgery. Polarization-sensitive optical coherence tomography (PS-OCT) extends upon conventional OCT by measuring the polarization state of backscattered light in order to estimate the optical properties of the tissues, such as birefringence and optic axis. 

Birefringence is a useful contrast mechanism in biological tissues due to its increased contrast for fibrillar tissues, such as collagen, muscle and myelinated nerve fibers which are vital to both normal architecture as well as several pathological processes. The orientation of these fibers can also be found by accurately determining the depth-resolved optic axis orientation of the tissue, providing unprecidented detailed mapping of tissues in-vivo. However, this process can be hindered in catheter-based applications due to the deleterious effects of several factors including the catheter tortuosity, catheter sheath retardance and non-uniform rotational distortion (NURD). This processing pipeline overcomes these effects by capitalizing on the catheter sheath as a static "guide-star" to allow for compensation of these effects.

## 2. GuideStar Processing Overview
The GuideStar method enables accurate characterization of the optic axis of tissues as described in [1]. Briefly, we have developed a novel calibration approach to rigorously reconstruct the absolute, depth-resolved optic axis orientation by utilizing the retardance of the catheter sheath as an optic axis orientation guide star signal for each probe orientation, assuming that the sheath acts as a linear retarder with a constant optic axis orientation. This method relies on both 1) accurate segmentation of the inner and outer sheath boundaries and 2) decomposition of the polarimetric signals and compensation for the effects of the catheter. In the code, these two elements are separated into folders

- **Segmentation:** Contains sheath segmentation relevent functions
- **AbsolutLocalOA:** Contains polarization effect compensation and OA estimation relevant functions

These two elements are combined together in a test script **SampleScript.m**

## 3. Example data description

