
function [xNdm, xElm, xIntNd, xBdyN, xBdyD, xBdyDP, xNdcond]...
    = feSetHnodesQuadratic(NdmF, ElmF, IntNdF, BdyNF, BdyDF, BdyDPF, NdcondF)

  NndF = size(NdmF, 1);            % [ID, X, Y, U, NAdjEdges]'
  NelmF = size(ElmF, 1) ;          % [Node1, Node2, Node3]'

  Nint = size(IntNdF, 1);
  NbdN = size(BdyNF, 1);
  NbdD = size(BdyDF, 1);

  NbdDP = size(BdyDPF, 1);
  Ncond = size(NdcondF, 1);
  
% change the variables so that we can use this code as a function routine.

NdEltCell = cell(NndF, 2); % {Element No.}{th on the Element}
Map2New = zeros(NndF, 1);  % Map from old node to new extended node
% % Map2Old = zeros(NndF, 1);  % Map from new node to old extended node

for ie = 1:NelmF
    
    for jt = 1:3
        ndi = ElmF(ie, jt);

        % Set Elt/th on each node
        NdEltCell{ndi, 1} = [NdEltCell{ndi, 1}; ie];
        NdEltCell{ndi, 2} = [NdEltCell{ndi, 2}; jt];
    end
    
end

ElcondEl = zeros(NelmF, 3, 2);
for ie = 1:NelmF

    for th1 = 1:3

        if th1 == 3 th2 = 1;
        else th2 = th1+1;
        end

        % Two nodes which we need to seek a element on.
        NodeH = ElmF(ie, th1);
        NodeT = ElmF(ie, th2);
        
        % Array of elements that contain  NodeH. in Th position
        EltAr = NdEltCell{NodeH, 1};
        EthAr = NdEltCell{NodeH, 2};
        NuEH = size(EltAr, 1);

        % Find element which has NodeT as a node
        for in = 1:NuEH
            elt = EltAr(in, 1);
            if elt ~= ie
                for id = 1:3
                    if id ~= EthAr(in, 1)
                        if ElmF(elt, id) == NodeT

                            ElcondEl(ie, th1, 1) = elt;
                            ElcondEl(ie, th1, 2) = id;

                        end
                    end
                end
            end
        end
    end
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Depends on the order of approximation

% Quadratic basis function
NuNdElt = 6;



xElm = zeros(NelmF, NuNdElt);

cntNd = 0;

for i = 1:NelmF
    
    for j = 1:NuNdElt

        nod = xElm(i, j);
        if nod == 0
            cntNd = cntNd + 1;
            xElm(i, j) = cntNd;
        
            if j <= 3
                elold = ElmF(i, j);
                
                % Define Old/New Mappings
                Map2New(elold, 1) = cntNd;
% %                 Map2Old(cntNd, 1) = elold;
                
                Nsize = size(NdEltCell{elold, 1}, 1);

                emt = NdEltCell{elold, 1};
                eth = NdEltCell{elold, 2};

                for k = 1:Nsize 
                    qi = emt(k, 1);
                    qj = eth(k, 1);
                    xElm(qi, qj) = cntNd;
                end

            elseif j <= 6

                tEl = ElcondEl(i, j-3, 1);
                tTh = ElcondEl(i, j-3, 2);
                if tEl ~= 0
                    xElm(tEl, tTh+3) = cntNd;
                end
            end
        end
    end
end


% Map from old node to new node ordering(# is for oldnodes)
% Map2New
% % Map2Old

xNdm = zeros(cntNd, 4);


for i = 1:NelmF    % 3 
    ndId = ElmF(i, :);
    tmpVec = zeros(3,4);
    
    for t = 1:3
        tmpVec(t,:) = NdmF(ndId(1, t),:);
    end

    for j = 1:NuNdElt % 6
        if j <= 3
            ndold = ElmF(i,j);
            ndnew = xElm(i,j);
            if xNdm(ndnew, 4) == 0  % flag:->Never assigned
                xNdm(ndnew, 1) = ndnew;
                xNdm(ndnew, 2) = tmpVec(j, 2);
                xNdm(ndnew, 3) = tmpVec(j, 3);
                xNdm(ndnew, 4) = 1;
            end
        elseif j <= 6
            gNd = xElm(i,j);
            if xNdm(gNd, 4) == 0
                nd1 = j-3;
                nd2 = j-2;
                if nd2 >= 4
                    nd2 = 1;
                end
                xNdm(gNd, 1) = gNd;
                xNdm(gNd, 2) = 0.5 * (tmpVec(nd1, 2) + tmpVec(nd2, 2));
                xNdm(gNd, 3) = 0.5 * (tmpVec(nd1, 3) + tmpVec(nd2, 3));
                xNdm(gNd, 4) = 1;
            end
        end
    end
end


% Edge Information #Positive:Internal Edge, #0:Bd-N Edge, #-1:Bd-D Edge
Egm = zeros(NelmF, 3);
Egm(:,:) = ElcondEl(:,:,1);    % This shows Internal Edge/Boundary Edge, doesn't show Bd-D/Bd-N

% Update Egm and set BdyNF/BdyDF
for i = 1:NbdD
    
    NmNode = BdyDF(i, 1);
    emt = NdEltCell{NmNode, 1};
    eth = NdEltCell{NmNode, 2};
    
    Nnbd = size(emt, 1);
    for j = 1:Nnbd
        p1 = emt(j, 1);
        p2 = eth(j, 1);
        NbElt = ElcondEl(emt(j, 1), eth(j, 1));
        if NbElt == 0
            Egm(p1, p2) = -1; % This edge is Dirichlet Boundary Edge
            break;
        end
    end
    
end


% Extended Edge/order information for each node.
NdqEltCell = cell(cntNd, 3); % {Element No.}{th on the Element}

% Geometric Nodes processing
for in = 1:NndF    % copy the content of NdEltCell into NdqElltCell

    newIn = Map2New(in);
    NdqEltCell{newIn, 1} = NdEltCell{in, 1};
    NdqEltCell{newIn, 2} = NdEltCell{in, 2};
    NdqEltCell{newIn, 3} = zeros(size(NdEltCell{in, 1})); % This means the node is on geometric node

end
for je = 1:NelmF
    for ke = 1:3
        
        nod = xElm(je, ke+3);
        
        if size(NdqEltCell{nod, 1}, 1) == 0
            %Self info
            NdqEltCell{nod, 1} = [NdqEltCell{nod, 1}; je];
            NdqEltCell{nod, 2} = [NdqEltCell{nod, 2}; ke];
            NdqEltCell{nod, 3} = [NdqEltCell{nod, 3}; 0.5];
 
            %Nbd info
            nbElt = ElcondEl(je, ke, 1);
            if nbElt ~= 0
                NdqEltCell{nod, 1} = [NdqEltCell{nod, 1}; nbElt];
                NdqEltCell{nod, 2} = [NdqEltCell{nod, 2}; ElcondEl(je, ke, 2)];
                NdqEltCell{nod, 3} = [NdqEltCell{nod, 3}; 0.5];
            end

        end
    end
end


% Table for updating nodes in partitioning
locTbl = zeros(cntNd, 1);
NbdDq = 2*NbdD;
bdDTbl = zeros(NbdDq, 1);  % Num of all bdy is 2*origial bdyD

xIntNd = IntNdF;
xBdyN = BdyNF;
xBdyD = BdyDF;

% Update old node indices to new node 
for i = 1:Nint
    newNode = Map2New(IntNdF(i, 1), 1);
    xIntNd(i, 1) = newNode;
    locTbl(newNode, 1) = 1;
end

for i = 1:NbdN
    newNode = Map2New(BdyNF(i, 1), 1);
    xBdyN(i, 1) = newNode;
    locTbl(newNode, 1) = 1;
end

for i = 1:NbdD
    newNode = Map2New(BdyDF(i, 1), 1);
    xBdyD(i, 1) = newNode;
    bdDTbl(i, 1) = 1;
    locTbl(newNode, 1) = 1;
end

condTbl = locTbl; % Copy locTbl which original node are set to 1

for i = 1:NelmF
    for j = 1:3
        newNode = xElm(i, j+3);
        edgeInfo = Egm(i,j);
        if edgeInfo > 0         % Interior Edge
            if locTbl(newNode, 1) == 0
                xIntNd = [xIntNd;newNode];
                locTbl(newNode, 1) = 1;
            end

        elseif edgeInfo == -1   % Dirichlet Bdy Edge
            if locTbl(newNode, 1) == 0
                xBdyD = [xBdyD;newNode];
                locTbl(newNode, 1) = 1;
            end
        else                    % Newmann Bdy Edge
            if locTbl(newNode, 1) == 0
                xBdyN = [xBdyN;newNode];
                locTbl(newNode, 1) = 1;
            end
        end
    end
end


% Copy content of BdyDPF to front positions of xBdyDP
xBdyDP = BdyDPF;

%Total Solution Set including D-Bdy conditions
Solvec = zeros(cntNd, 1);
for ip = 1:NbdD
    Solvec(xBdyD(ip, 1), 1) = BdyDPF(ip, 1);
end


% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % Module to be updated in dealing with unknown solutions
% % % % % % 
% % % % % % for iq = NbdD+1:NbdDq
% % % % % %     dnd = xBdyD(iq, 1);
% % % % % %     xv = xNdm(dnd, 2);
% % % % % %     yv = xNdm(dnd, 3);
% % % % % %     
% % % % % %     xBdyDP(iq, 1) = feInSolution(xv, yv);
% % % % % %     Solvec(xBdyD(iq, 1), 1) = xBdyDP(iq, 1);
% % % % % % end
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% In case we don't know exact solution:
for iq = NbdD+1:NbdDq
    dnd = xBdyD(iq, 1);
    ee = NdqEltCell{dnd, 1};
    e1 = NdqEltCell{dnd, 2};
    am = NdqEltCell{dnd, 3};
    
    if e1 == 3
        e2 = 1;
    else
        e2 = e1+1;
    end
    nd1 = xElm(ee, e1);
    nd2 = xElm(ee, e2);
    
    xBdyDP(iq, 1) = (1-am) * Solvec(nd1, 1) + am * Solvec(nd2, 1);
    Solvec(xBdyD(iq, 1), 1) = xBdyDP(iq, 1);
end

% In regards to conductivity, we just maintain it
xNdcond = NdcondF;

return
