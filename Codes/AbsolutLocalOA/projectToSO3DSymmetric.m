function Rout = projectToSO3DSymmetric(Min)
%Rout = projectToSO3(Min) takes Min, an array of 3x3 matrices, assumed to be
%D-transpose symmetric, and 'projects' them into SO3. 

%This is simply a helper function to avoid confusing with boolDSym
Rout = projectToSO3(Min,true);

