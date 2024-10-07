function [Rf,dep] = cleanCathRetTrace(Rin,filtwidth,depTh)
%function cleanRetTrace(Rin,filtwidth,depTh) takes the 3x3 SO3 matrix 
%Rin representing the cumulative roundtrip retardance of a sample or
%catheter interface as a function of catheter rotation and filters it with
%a filterwidth of filtwidth. Locations where the averaged SO3 matrices
%have depolarization below depTh are interpolated. Signal is assumed
%periodic.

% Assuming Rin is 3x3xNAlines

% build a filter
ff = shiftdim(hanning(filtwidth),-2);
ff = ff/sum(ff);

% filtered matrix; contains both retardance and depolarization
Mf = imfilter(Rin,ff,'circular');

% compute depolarization index of the averaged matrices
dep = sqrt(sum(sum(Mf.^2,1),2)/3);

% remove depolarizing part with polar decomposition
Rf = projectToSO3DSymmetric(Mf);

% identify areas below depTh
inds = dep<depTh;

%Rforig = Rf;

if sum(inds)>0
    % find continuous segments that need interpolation
    starting = find(inds&circshift(~inds,[0,0,1]));
    ending = find(inds&circshift(~inds,[0,0,-1]));

%    Rf2 = Rf;
    % linearly interpolate rotation matrices
    for ind = 1:numel(starting)
        startloc = mod(starting(ind)-2,size(Mf,3))+1;
        endloc = mod(ending(ind),size(Mf,3))+1;

        R1 = Rf(:,:,startloc);
        R2 = Rf(:,:,endloc);

        R1sq = makeRot3x3(decomposeRot(R1)/2);

%         dr = decomposeRot(R1sq.'*R2*R1sq.');
%         Rf(:,:,starting(ind):ending(ind)) = pagemtimes(R1sq,pagemtimes(makeRot3x3(dr*linspace(0,1,ending(ind)-starting(ind)+1)),R1sq));


        R2sq = makeRot3x3(decomposeRot(R2)/2);
        fractInd = 1/(ending(ind)-starting(ind)+1);
        deltaR = R2sq*R1sq.';
        % local derivatives at edges
        Rbefore = Rf(:,:,mod(startloc-2,size(Mf,3))+1);
        Rafter = Rf(:,:,mod(endloc,size(Mf,3))+1);

        fun = @(x)optVangle(x,deltaR,R1sq,fractInd,Rbefore,Rafter);

        vopt = fminsearch(fun,0);


        V = makeRot3x3([0;0;vopt]);
%        R2sq = makeRot3x3(decomposeRot(R2)/2);
        dr = decomposeRot(V*deltaR);
%        dr = decomposeRot(V*makeRot3x3(decomposeRot(R1sq.'*R2*R1sq.')/2));
        temp = pagemtimes(makeRot3x3(dr*linspace(0,1,ending(ind)-starting(ind)+1)),R1sq);
        Rf(:,:,starting(ind):ending(ind)) = pagemtimes(pagetranspose(temp).*[1,1,-1;1,1,-1;-1,-1,1],temp);


%         
% 
% 
%         V = makeRot3x3([0;0;1]);
%         drp = decomposeRot(R1sq.'*V*makeRot3x3(dr/2)*R1sq);
% 
%         (R1sq*(drp.*[1;1;-1])) + R1sq*drp
% 
%         decomposeRot(R1sq*dR1*(makeRot3x3(dr/2)*R1sq).')
%         decomposeRot(R1sq*dR2*(makeRot3x3(dr/2)*R1sq).')
% 
%         loc = makeRot3x3(drp*linspace(0,1,ending(ind)-starting(ind)+1));
%         Rf2(:,:,starting(ind):ending(ind)) = pagemtimes(pagetranspose(loc).*[1,1,-1;1,1,-1;-1,-1,1],pagemtimes(R1,loc));

    end

end


function out = optVangle(vangle,deltaR,R1sq,fractInd,Rbefore,Rafter)

    V = makeRot3x3([0;0;vangle]);
    dr = decomposeRot(V*deltaR);
    Ra = makeRot3x3(-dr*fractInd)*R1sq;
    Rb = makeRot3x3(dr*(1+fractInd))*R1sq;

    out = sum(sum(((Ra.'.*[1,1,-1;1,1,-1;-1,-1,1])*Ra-Rbefore).^2 + ((Rb.'.*[1,1,-1;1,1,-1;-1,-1,1])*Rb-Rafter).^2));

