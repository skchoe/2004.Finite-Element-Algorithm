
function [NdMat, Elm, BdyD, BdyN, IntNd, BdyDP, NdcondF] = ...
    feMeshTorsoTrs37(NdmF, ElmF, BdyDF, BdyNF, IntNdF, BdyDPF, NdcondF, drawFlag)

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
NuNdElt = 10;   % Number of nodes on an element
Ntri = 3;       % Number of edges on a triangle
Ndeg = 2;       % Number of nodes on the interior of an edge

xElm = zeros(NelmF, NuNdElt);

cntNd = 0;

for i = 1:NelmF

    % Ordering of node index in given Element from 1-3
    for j = 1:Ntri

        nod = xElm(i, j);
        if nod == 0
            cntNd = cntNd + 1;
            xElm(i, j) = cntNd;
        
            elold = ElmF(i, j);
                
            % Define Old/New Mappings
            Map2New(elold, 1) = cntNd;
                
            Nsize = size(NdEltCell{elold, 1}, 1);

            emt = NdEltCell{elold, 1};
            eth = NdEltCell{elold, 2};

            for k = 1:Nsize 
                 qi = emt(k, 1);
                 qj = eth(k, 1);
                 xElm(qi, qj) = cntNd;
            end
        end
    end

    % Ordering of node index in given Element from 4-9
    for j = 1:Ntri
        
        th = Ntri + 2*j-1;  % Index 4,6,8
        nod1 = xElm(i, th);
        nod2 = xElm(i, th+1);
        
        if nod1 == 0 & nod2 == 0

            cN1 = cntNd+1;
            cN2 = cntNd+2;
            xElm(i, th) = cN1;
            xElm(i, th+1) = cN2;
            cntNd = cntNd + 2;

            tEl = ElcondEl(i, j, 1);
            tTh = ElcondEl(i, j, 2);
            if tEl ~= 0
                thN = Ntri + 2*tTh - 1;
                xElm(tEl, thN) = cN2;
                xElm(tEl, thN+1) = cN1;
            end
        end
    end
    
    % Ordering of node index in given Element from 10
    cntNd = cntNd+1;
    xElm(i, NuNdElt) = cntNd;
end



% Map from old node to new node ordering(# is for oldnodes)
% Map2New
% % Map2Old

xNdm = zeros(cntNd, 3);


for i = 1:NelmF    % 3 
    ndId = ElmF(i, :);  % Get old id for retreiving node geometry
    tmpVec = zeros(3,3);
    
    for t = 1:3
        tmpVec(t,:) = NdmF(ndId(1, t),:);
    end

    for j = 1:Ntri
        ndold = ElmF(i,j);
        ndnew = xElm(i,j);
        if xNdm(ndnew, 3) == 0  % flag:->Never assigned
            xNdm(ndnew, 1) = tmpVec(j, 1);
            xNdm(ndnew, 2) = tmpVec(j, 2);
            xNdm(ndnew, 3) = 1;  % set flag to 'on'
        end
    end
  

    for j = 1:Ntri
        th = Ntri + 2*j-1;  % Index 4,6,8
        gNd1 = xElm(i, th);
        gNd2 = xElm(i, th+1);
        if xNdm(gNd1, 3) == 0 & xNdm(gNd2, 3) == 0
            nd1 = j;
            nd2 = j+1;
            if nd2 >= 4
                nd2 = 1;
            end
            xNdm(gNd1, 1) = 2/3 * tmpVec(nd1, 1) + 1/3 * tmpVec(nd2, 1);
            xNdm(gNd1, 2) = 2/3 * tmpVec(nd1, 2) + 1/3 * tmpVec(nd2, 2);
            xNdm(gNd1, 3) = 1;

            xNdm(gNd2, 1) = 1/3 * tmpVec(nd1, 1) + 2/3 * tmpVec(nd2, 1);
            xNdm(gNd2, 2) = 1/3 * tmpVec(nd1, 2) + 2/3 * tmpVec(nd2, 2);
            xNdm(gNd2, 3) = 1;
        end
    end
    
    gNd = xElm(i, NuNdElt);
    xNdm(gNd, 1) = 1/3 * (sum(tmpVec(:,1)));
    xNdm(gNd, 2) = 1/3 * (sum(tmpVec(:,2)));
    xNdm(gNd, 3) = 1;
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
NdqEltCell = cell(cntNd, 3); % {Element No.}{th on the Element}{flag}

% Geometric Nodes processing
for in = 1:NndF    % copy the content of NdEltCell into NdqElltCell for each geometric node
    
    newIn = Map2New(in);
    NdqEltCell{newIn, 1} = NdEltCell{in, 1};
    NdqEltCell{newIn, 2} = NdEltCell{in, 2};
    NdqEltCell{newIn, 3} = zeros(size(NdEltCell{in, 1})); % This means the node is on geometric node
    
end

% Approximation Nodes processing
for je = 1:NelmF
    
    for ke = 1:3 % Number of edges
        
        th = Ntri + 2*ke-1;
        nod1 = xElm(je, th);
        
        if size(NdqEltCell{nod1, 1}, 1) == 0
            %Self info
            NdqEltCell{nod1, 1} = [NdqEltCell{nod1, 1}; je];
            NdqEltCell{nod1, 2} = [NdqEltCell{nod1, 2}; ke];
            NdqEltCell{nod1, 3} = [NdqEltCell{nod1, 3}; 1/3];
 
            %Nbd info
            nbElt = ElcondEl(je, ke, 1);
            if nbElt ~= 0
                NdqEltCell{nod1, 1} = [NdqEltCell{nod1, 1}; nbElt];
                NdqEltCell{nod1, 2} = [NdqEltCell{nod1, 2}; ElcondEl(je, ke, 2)];
                NdqEltCell{nod1, 3} = [NdqEltCell{nod1, 3}; 2/3];
            end

        end
        
        nod2 = xElm(je, th+1);
        
        if size(NdqEltCell{nod2, 1}, 1) == 0
            %Self info
            NdqEltCell{nod2, 1} = [NdqEltCell{nod2, 1}; je];
            NdqEltCell{nod2, 2} = [NdqEltCell{nod2, 2}; ke];
            NdqEltCell{nod2, 3} = [NdqEltCell{nod2, 3}; 2/3];
 
            %Nbd info
            nbElt = ElcondEl(je, ke, 1);
            if nbElt ~= 0
                NdqEltCell{nod2, 1} = [NdqEltCell{nod2, 1}; nbElt];
                NdqEltCell{nod2, 2} = [NdqEltCell{nod2, 2}; ElcondEl(je, ke, 2)];
                NdqEltCell{nod2, 3} = [NdqEltCell{nod2, 3}; 1/3];
            end

        end
    end
end

% Testing NdqEltCell 
% % % node = 38;
% % % a = NdqEltCell{node,1}
% % % b = NdqEltCell{node,2}
% % % c = NdqEltCell{node,3}

% Table for updating nodes in partitioning
locTbl = zeros(cntNd, 1);
NbdDq = 3*NbdD; % Since we have 3 D-bdy nodes in each edge 
bdDTbl = zeros(NbdDq, 1);  % Num of all bdy is 3*origial bdyD

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
        
        th = Ntri + 2*j-1;
        newNode1 = xElm(i, th);
        newNode2 = xElm(i, th+1);
        
        edgeInfo = Egm(i,j);
        
        if edgeInfo > 0         % Interior Edge
            if locTbl(newNode1, 1) == 0
                xIntNd = [xIntNd;newNode1];
                locTbl(newNode1, 1) = 1;
            end
            if locTbl(newNode2, 1) == 0
                xIntNd = [xIntNd;newNode2];
                locTbl(newNode2, 1) = 1;
            end
            
        elseif edgeInfo == -1   % Dirichlet Bdy Edge
            if locTbl(newNode1, 1) == 0
                xBdyD = [xBdyD;newNode1];
                locTbl(newNode1, 1) = 1;
            end
            if locTbl(newNode2, 1) == 0
                xBdyD = [xBdyD;newNode2];
                locTbl(newNode2, 1) = 1;
            end
       else                    % Newmann Bdy Edge
            if locTbl(newNode1, 1) == 0
                xBdyN = [xBdyN;newNode1];
                locTbl(newNode1, 1) = 1;
            end
            if locTbl(newNode2, 1) == 0
                xBdyN = [xBdyN;newNode2];
                locTbl(newNode2, 1) = 1;
            end
       end
    end
end

% Insert center node to interior set xIntNd
for i = 1:NelmF
    xIntNd = [xIntNd;xElm(i, NuNdElt)];
end


% Copy content of BdyDPF to front positions of xBdyDP
xBdyDP = BdyDPF;

%Total Solution Set including D-Bdy conditions
Solvec = zeros(cntNd, 1);
for ip = 1:NbdD
    Solvec(xBdyD(ip, 1), 1) = BdyDPF(ip, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Module to be updated in dealing with unknown solutions

for iq = NbdD+1:NbdDq
    dnd = xBdyD(iq, 1);
    xv = xNdm(dnd, 1);
    yv = xNdm(dnd, 2);
    
    xBdyDP(iq, 1) = feInSolution(xv, yv);
    Solvec(xBdyD(iq, 1), 1) = xBdyDP(iq, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % % % In case we don't know exact solution:
% % % % % % 
% % % % % % for iq = NbdD+1:NbdDq
% % % % % %     dnd = xBdyD(iq, 1);
% % % % % %     ee = NdqEltCell{dnd, 1};
% % % % % %     e1 = NdqEltCell{dnd, 2};
% % % % % %     am = NdqEltCell{dnd, 3};
% % % % % %     
% % % % % %     if e1 == 3
% % % % % %         e2 = 1;
% % % % % %     else
% % % % % %         e2 = e1+1;
% % % % % %     end
% % % % % %     nd1 = xElm(ee, e1);
% % % % % %     nd2 = xElm(ee, e2);
% % % % % %     
% % % % % %     xBdyDP(iq, 1) = (1-am) * Solvec(nd1, 1) + am * Solvec(nd2, 1);
% % % % % %     Solvec(xBdyD(iq, 1), 1) = xBdyDP(iq, 1);
% % % % % % end
% % % % % % 

% Assigning conductivity matrix to new node set xNdm
xNdcond = zeros(cntNd, 4);
conFlg = zeros(cntNd, 1);

for ic = 1:NndF
    newNd = Map2New(ic, 1);
    xNdcond(newNd, :) = NdcondF(ic, :);
    conFlg(newNd, 1) = 1;
end

for ie = 1:NelmF
    for je = 1:3
        th = Ntri + 2*je-1;
        
        nd1 = xElm(ie, th);
        nd2 = xElm(ie, th+1);
        if conFlg(nd1, 1) == 0 & conFlg(nd2, 1) == 0
            conEltVec = NdqEltCell{nd1, 1};
            conThVec = NdqEltCell{nd1, 2};
            conFcVec = NdqEltCell{nd1, 3};
            
            ee = conEltVec(1, 1);
            e1 = conThVec(1, 1);
            am = conFcVec(1, 1);
            
            if e1 == 3
                e2 = 1;
            else
                e2 = e1+1;
            end
            node1 = xElm(ee, e1);
            node2 = xElm(ee, e2);

%             % Finding conductivity in interpolation nodes(by interpolating between known values)
%             for im = 1:4 % [2x2] matrix form of conductivity
%                 xNdcond(nd1, im) = (1-am) * xNdcond(node1, im) + am * xNdcond(node2, im);
%                 xNdcond(nd2, im) = am * xNdcond(node1, im) + (1-am) * xNdcond(node2, im);
%             end

            % Finding conductivity in interpolation nodes(by exact evaluating)
            cx1 = xNdm(nd1, 1); cy1 = xNdm(nd1, 2);
            cx2 = xNdm(nd2, 1); cy2 = xNdm(nd2, 2);
            
            cdtmat1 = feInConduct(cx1, cy1);
            cdtmat2 = feInConduct(cx2, cy2);
            for im = 1:4 % [2x2] matrix form of conductivity
                xNdcond(nd1, im) = cdtmat1(1,1,im);
                xNdcond(nd2, im) = cdtmat2(1,1,im);
            end
            
            
            conFlg(nd1, 1) = 1;
            conFlg(nd2, 1) = 1;
        end
    end
    % Conductivity on the center node of each element
    nd0 = xElm(ie, NuNdElt);
    nodes = zeros(Ntri, 1);
    for i = 1:Ntri
        nodes(i, 1) = xElm(ie, i);
    end
    for im = 1:4
        xNdcond(nd0, im) = ( xNdcond(nodes(1, 1), im) +  xNdcond(nodes(2, 1), im) + xNdcond(nodes(3, 1), im) )/3;
    end
end

NdMat = xNdm;


Elm = zeros(9*NelmF, 3);

for iel = 1:NelmF
    cMl = 9*(iel-1);

    cT = cMl + 1;
    Elm(cT, 1) = xElm(iel, 1); 
    Elm(cT, 2) = xElm(iel, 4);        
    Elm(cT, 3) = xElm(iel, 9);

    
    cT = cMl + 2;
    Elm(cT, 1) = xElm(iel, 2); 
    Elm(cT, 2) = xElm(iel, 6);        
    Elm(cT, 3) = xElm(iel, 5);

   
    cT = cMl + 3;
    Elm(cT, 1) = xElm(iel, 3); 
    Elm(cT, 2) = xElm(iel, 8);        
    Elm(cT, 3) = xElm(iel, 7);


    cT = cMl + 4;
    Elm(cT, 1) = xElm(iel, 4); 
    Elm(cT, 2) = xElm(iel, 5);        
    Elm(cT, 3) = xElm(iel, 10);


    cT = cMl + 5;
    Elm(cT, 1) = xElm(iel, 6); 
    Elm(cT, 2) = xElm(iel, 7);        
    Elm(cT, 3) = xElm(iel, 10);
 
    
    cT = cMl + 6;
    Elm(cT, 1) = xElm(iel, 8); 
    Elm(cT, 2) = xElm(iel, 9);        
    Elm(cT, 3) = xElm(iel, 10);
 
   
    cT = cMl + 7;
    Elm(cT, 1) = xElm(iel, 4); 
    Elm(cT, 2) = xElm(iel, 10);        
    Elm(cT, 3) = xElm(iel, 9);


    cT = cMl + 8;
    Elm(cT, 1) = xElm(iel, 6); 
    Elm(cT, 2) = xElm(iel, 10);        
    Elm(cT, 3) = xElm(iel, 5);


    cT = cMl + 9;
    Elm(cT, 1) = xElm(iel, 8); 
    Elm(cT, 2) = xElm(iel, 10);        
    Elm(cT, 3) = xElm(iel, 7);

end



BdyD    = xBdyD;
BdyN    = xBdyN;
IntNd   = xIntNd;
BdyDP   = xBdyDP;
NdcondF = xNdcond;


% % Coordinates
% save mat/data_pts.mat NdMat % -ASCII
% 
% % Elements
% save mat/data_elt.mat Elm % -ASCII
% 
% % Dirichlet Boundary Coordinate Indices
% save mat/data_bdyD.mat BdyD  %-ASCII
% 
% % Neumann Boundary Coordinate Indices
% save mat/data_bdyN.mat BdyN  %-ASCII
% 
% % Interior Coordinate Indices
% save mat/data_intNd.mat IntNd  %-ASCII
% 
% % Boundary Potential on Dirichlet Boundary
% save mat/data_bdyDP.mat BdyDP %-ASCII
% 
% % Conductivity data for each node
% save mat/data_condc.mat NdcondF %-ASCII

% % % % % % Visualization
if drawFlag == 1
  figno = 24;
  subPlotData(figno, cntNd, NdMat, Elm, BdyN, BdyD);
end 

return


function subPlotData(figno, Nnd, NdmF, ElmF, BdyN, BdyD) 

  Xvec = NdmF(:,1)';
  Yvec = NdmF(:,2)';

  NelF = size(ElmF, 1);

  figure(figno);
    triplot(ElmF, Xvec, Yvec, 'Color', [0.5 0.2 0.5]); hold on;

% %     Plot global nodes
%     for i = 1:Nnd
%         str = num2str(i);
%         text(Xvec(1, i), Yvec(1, i), str, 'Color', [0.5 0 0], 'Fontsize', 8);
%     end

%     hold on;
%     
% %     Plot Element No.
%     for i = 1:NelF
%         one = ElmF(i,1); two = ElmF(i,2); thr = ElmF(i,3);
%         xp = (Xvec(1, one) + Xvec(1, two) + Xvec(1, thr))/3;
%         yp = (Yvec(1, one) + Yvec(1, two) + Yvec(1, thr))/3;
%         str = num2str(i);
%         text(xp, yp, str, 'Color', [0 0.5 0], 'Fontsize', 7);
%     end
% 
%     hold on;


  Nbd = size(BdyD, 1);
  Nbn = size(BdyN, 1);

  for i = 1:Nbd
      px = NdmF(BdyD(i, 1), 1);
      py = NdmF(BdyD(i, 1), 2);
      
      plot(px, py, 'd', 'Color', [1 0 0]);
  end

return
