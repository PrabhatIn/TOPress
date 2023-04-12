function TOPress_external_arch(nelx,nely,volfrac,penal,rmin,etaf,betaf,lst,maxit)
%% ___PART 1.____________________________MATERIAL AND FLOW PARAMETERS
E1 = 1;
Emin = E1*1e-6;
nu = 0.30;
[Kv,epsf,r,Dels] = deal(1,1e-7,0.1,2);              % flow parameters
[Ds, kvs]= deal((log(r)/Dels)^2*epsf,Kv*(1 - epsf));                    % flow parameters
%% ____PART 2._______________FINITE ELEMENT ANALYSIS PREPARATION and NON-DESIGN DOMAIN
[nel,nno] = deal(nelx*nely, (nelx+1)*(nely+1));
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
Udofs = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
[Lnode,Rnode,]= deal(1:nely+1, (nno-nely):nno);
[Bnode,Tnode]= deal((nely+1):(nely+1):nno, 1:(nely+1):(nno-nely));
[Pdofs,allPdofs, allUdofs] = deal(Udofs(:,2:2:end)/2,1:nno,1:2*nno);
iP = reshape(kron(Pdofs,ones(4,1))',16*nel,1);
jP = reshape(kron(Pdofs,ones(1,4))',16*nel,1);
iT = reshape(kron(Udofs,ones(4,1))',32*nel,1);
jT = reshape(kron(Pdofs,ones(1,8))',32*nel,1);
iK = reshape(kron(Udofs,ones(8,1))',64*nel,1);
jK = reshape(kron(Udofs,ones(1,8))',64*nel,1);
Kp = 1/6*[4 -1 -2 -1;-1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4]; % flow matrix: Darcy Law
KDp = 1/36*[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]; % Drainage matrix
Te = 1/12*[-2 2 1 -1;-2 -1 1 2;-2 2 1 -1;-1 -2 2 1;-1 1 2 -2; -1 -2 2 1; -1 1 2 -2; -2 -1 1 2]; % transformation matrix
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
ke = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);             %stiffness matrix
IFprj=@(xv,etaf,betaf)((tanh(betaf*etaf) + tanh(betaf*(xv-etaf)))/...            %projection function
    (tanh(betaf*etaf) + tanh(betaf*(1 - etaf))));
dIFprj=@(xv,etaf,betaf) betaf*(1-tanh(betaf*(xv-etaf)).^2)...
    /(tanh(betaf*etaf)+tanh(betaf*(1-etaf)));                    % derivative of the projection function
[NDS, NDV ] = deal( [], [] );
act = setdiff((1 : nel)', union( NDS, NDV ));
%% ____PART 3.______PRESSURE & STRUCTURE B.C's, LOADs
[PF, Pin] =deal(0.00001*ones(nno,1),1); %pressure-field preparation
PF([Tnode,Lnode,Rnode, Bnode(1:nelx/10), Bnode(end-nelx/10:end)]) = Pin; PF(Bnode(nelx/10+1:end-nelx/10-1)) = 0; % applying pressure load
fixedPdofs = allPdofs(PF~=0.00001);
freePdofs  = setdiff(allPdofs,fixedPdofs);
pfixeddofsv = [fixedPdofs' PF(fixedPdofs)]; % p-fixed and its value
fixedUdofs = [2*Bnode(nelx/10+1)-1  2*Bnode(nelx/10+1)  2*Bnode(end-nelx/10-1)-1 2*Bnode(end-nelx/10-1)]; %fixed displ.
freeUdofs = setdiff(allUdofs,fixedUdofs);
[U, lam1] = deal(zeros(2*nno,1),zeros(nno,1));
%% ___PART 4._________________________________________FILTER PREPARATION
[dy,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
h = max( 0, rmin-sqrt( dx.^2 + dy.^2 ) );                                % conv. kernel
Hs = imfilter(ones(nely, nelx ), h);                               % matrix of weights (filter)
%% ___PART 5.__________________________MMA OPTIMIZATION PREPARATION & INITIALIZATION
[x,dVol0] = deal(zeros(nel,1),ones(nel,1)/(nel*volfrac));
 x(act) = (volfrac*(nel-length(NDV))-length(NDS) )/length(act); x(NDS) = 1; % updating the variable
[nMMA,mMMA,xphys,xMMA,mvLt] = deal(length(act),1,x,x(act),0.1);
[xminvec,xmaxvec] = deal(zeros(nMMA,1),ones(nMMA,1)); %Min. & Max
[low, upp] = deal(xminvec,xmaxvec); % Low and Upp limits MMA
[cMMA,dMMA, a0, aMMA] = deal(1000*ones(mMMA,1),zeros(mMMA,1),1,zeros(mMMA,1));
[xold1,xold2] = deal(xMMA);
[loop, change] =deal(0,1);
dVol = imfilter(reshape(dVol0, nely, nelx)./Hs,h);  %filtered volume sensitivity
%% ____PART 6._____________________________________MMA OPTIMIZATION LOOP
while(loop<maxit && change>0.001)
    loop = loop + 1;                                                   % Updating the opt. iteration
    %___PART 6.1__________SOLVING FLOW BALANCE EQUATION
    Kc = Kv*(1-(1-epsf)*IFprj(xphys,etaf,betaf));         %Flow coefficient
    Dc = Ds*IFprj(xphys,etaf,betaf);                            %Drainage coefficient
    Ae = reshape(Kp(:)*Kc' + KDp(:)*Dc',16*nel,1);  %Elemental flow matrix
    AG = (sparse(iP,jP,Ae)+ sparse(iP,jP,Ae)')/2;      %Global flow matrix
    Aff = AG(freePdofs,freePdofs);                           %AG for free pressure dofs
    PF(freePdofs) = decomposition(Aff,'ldl')\(-AG(freePdofs,fixedPdofs)*pfixeddofsv(:,2));
    PF(pfixeddofsv(:,1)) = pfixeddofsv(:,2);              % Final P-field
    %__PART 6.2_DETERMINING CONSISTENT NODAL LOADS and GLOBAL Disp. Vector
    Ts = reshape(Te(:)*ones(1,nel), 32*nel, 1);        %Elemental transformation matrix
    TG = sparse(iT, jT, Ts);                                       %Global transformation matrix
    F = -TG*PF;                                                       % Dertmining nodal forces
    E = Emin + xphys.^penal*(E1 - Emin);                %Material interpolation
    Ks = reshape(ke(:)*E',64*nel,1);                         %Elemental stiffness matrix
    KG = (sparse(iK,jK,Ks) + sparse(iK,jK,Ks)')/2;    %Global stiffnes matrix
    U(freeUdofs) = decomposition(KG(freeUdofs,freeUdofs),'chol','lower')\F(freeUdofs); %Global Disp. Vect.
    %__PART 6.3__OBJECTIVE, CONSTRAINT and THEIR SENSITIVITIES COMPUTATION
    obj = U'*KG*U; % determining objective
    lam1(freePdofs) = (2*U(freeUdofs)'*TG(freeUdofs,freePdofs))/Aff;
    objsT1 = -(E1 - Emin)*penal*xphys.^(penal - 1).*sum(([U(Udofs)]*ke).*[U(Udofs)],2);
    dC1k = -dIFprj(xphys,etaf,betaf).* sum((lam1(Pdofs)*(kvs*Kp)) .* PF(Pdofs),2);
    dC1h =  dIFprj(xphys,etaf,betaf).* sum((lam1(Pdofs)*(Ds*KDp)) .* PF(Pdofs),2);
    objsT2 = dC1k + dC1h; objsens = (objsT1 + lst*objsT2); % final sensitivities
    Vol = sum(xphys)/(nel*volfrac)-1;                                          % volume fraction
    if(loop ==1), normf = 1000/(obj);save normf normf;end, load normf normf;
    objsens = imfilter(reshape(objsens*normf, nely, nelx)./Hs,h);
    %___PART 6.4______________________SETTING and CALLING MMA OPTIMIZATION
    xval = xMMA;
    [xminvec, xmaxvec]= deal(max(0, xval - mvLt),min(1, xval + mvLt));
    [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(mMMA,nMMA,loop,xval,xminvec,xmaxvec,xold1,xold2, ...
        obj*normf,objsens(act),objsens(act)*0,Vol,dVol(act),dVol(act)*0,low,upp,a0,aMMA,cMMA,dMMA);
    [xold2,xold1, xnew]= deal(xold1, xval,xmma);
    change = max(abs(xnew-xMMA)); % Calculating change
    xMMA = xnew;   xphys(act) = xnew;xphys = imfilter(reshape(xphys, nely, nelx),h)./Hs; 
    xphys = xphys(:); xphys(NDS) = 1; xphys(NDV) = 0; 
    %___PART 6.5_____________________________Printing and plotting results
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,obj*normf,mean(xphys),change);
    colormap(gray); imagesc(1-reshape(xphys, nely, nelx));caxis([0 1]);axis equal off;drawnow;
end
%%____________________________Printing the final result with pressure field
[xx,yy,nod_dV] = deal(0:nelx, 0:-1:-nely,zeros(nno,1));
convert =[1; 1/2*ones(nely-1,1);1;reshape([1/2*ones(nelx-1,1),1/4*ones(nelx-1,nely-1)...
    ,1/2*ones(nelx-1,1)]',[],1);1; 1/2*ones(nely-1,1); 1];
PFP = figure(2); set(PFP,'color','w');
for i = 1:nel, nod_dV(Pdofs(i,:)) = nod_dV(Pdofs(i,:)) + xphys(i); end
xphysnodal = full(reshape(nod_dV.*convert,nely+1,nelx+1)); % nodal density grid
PnodMat   = -full(reshape((PF)/Pin, nely+1, nelx+1));% nodal pressure field
winter1 = [ linspace(0,0.85)' linspace(0.8078,0.65)' linspace(0.9,0.15)']; % colorbar style
colormap([winter1; 0. 0. 0. ]); black = 0.01* ones(nely+1, nelx+1);
PrDenHand = surf(xx,yy,PnodMat, 'FaceColor', 'interp'); caxis([-1 0.0]);hold on
Hbl = imagesc(xx, yy, black, [-1 0.01]);  shading interp;xlim ([0 nelx]); ylim([-nely, 0]); view(0,90);
PwDenHand = Hbl; axis equal off; set(PwDenHand, 'AlphaData', xphysnodal);
set(PrDenHand, 'CData', PnodMat); set(PrDenHand,'ZData',PnodMat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    The code, TOPress,  is provided for  pedagogical purposes. A  detailed                           %
%    description is presented in the paper:"TOPress: a MATLAB implementation for               %
%    topology optimization of structures subjected to design-dependent pressure                    %
%    loads"   Structural and Mutltidisciplinary Optimization, 2023.                                             %
%                                                                                                                                               %
%    One can download the code and its extensions for the different problems                        %
%    from the online supplementary material and also from:                                                     %
%                                      https://github.com/PrabhatIn/TOPress                                             %
%                                                                                                                                               %
%    Please send your comment to: pkumar@mae.iith.ac.in                                                    %
%                                                                                                                                               %
%    One may also refer to the following two papers for more detail:                                        % 
%                                                                                                                                               %
%    1. Kumar P, Frouws JS, Langelaar M (2020) Topology optimization of fluidic                    %
%    pressure-loaded structures and compliant mechanisms using the Darcy method.            %
%    Structural and Multidisciplinary Optimization 61(4):1637-1655                                          %
%    2. Kumar P, Langelaar M (2021) On topology optimization of design-dependent              % 
%    pressure-loaded three-dimensional structures and compliant mechanisms.                     %
%    International Journal for Numerical Methods in Engineering 122(9):2205-2220                %
%                                                                                                                                               %
%    Disclaimer:                                                                                                                         %
%    The author does not guarantee that the code is free from erros but reserves                   %
%    all rights. Further, the author shall not be liable in any event caused by                           % 
%    use of the above 100-line code and its extensions                                                            %
%                                                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%