function Top_mfse2D_note(nptx,npty,refine,volfra,corlencoe)
%% 定义参数 %%
E0 = 1; Emin = 1e-9*E0; nu = 0.3; ptdist = 1;
elelen = ptdist/refine; nelx = refine*nptx; nely = refine*npty;
tolne = nelx*nely; tolnd = (nelx+1)*(nely+1); tolnpt = nptx*npty;
tolvol = tolne*elelen^2;
fprintf([' Number of material-field  points:%10i \n'...
    ' Number of finite elements:%10i\n'],tolnpt,tolne);

%% 准备有限元分析 %%
nodenrs = reshape(1:tolnd,1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,tolne,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1], tolne,1);
iK = reshape(kron(edofMat,ones(8,1))', 64*tolne,1);
jK = reshape(kron(edofMat,ones(1,8))', 64*tolne,1);
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE0 = 1/(1-nu^2)/24*([A11 A12; A12' A11]+nu*[B11 B12; B12' B11]);
% 定义边界条件
F = sparse(2, 1, -1, 2*tolnd, 1); U = zeros(2*tolnd,1);
fixeddofs = [1:2:2*(nely + 1),2*(nelx + 1)*(nely + 1),2*(nelx)*(nely + 1)];
freedofs = setdiff(1:2*tolnd,fixeddofs);

%% 材料场级数展开 %%
[eIntopMat,ptIntopMat] = MFSE2D(nptx,npty,refine,corlencoe,1);
save('MFSE2D.mat','eIntopMat','ptIntopMat');
load('MFSE2D.mat');

%% 初始化设计变量 %%
beta = 0.5; penal = 3;
x = (-log(1/volfra-1)/beta)*ones(1,tolnpt)/ptIntopMat; x = x';

%% 初始化迭代 %%
loop = 0; obj = 0.; neig = length(x);
change = 1.; ichange = 1; n = neig;
xmin = -1000*ones(n,1); xmax = 1000*ones(n,1);
low = xmin; upp = xmax;
xold1 = x;  xold2 = x; clf;
Obj = []; Volf = [];
[Xe, Ye] = meshgrid((0.5:1:nelx)*elelen, (nely-0.5:-1:0.5)*elelen);

%% 开始迭代 %%
while (change >= 0.005 || beta < 20)
    loop = loop + 1; objold = obj;
    %% 投影到密度场 %%
    ePhi = eIntopMat'*x;
    [ePhiProj,edproj] = THRESHOLD(ePhi,beta);

    %% 做有限元分析 %%
    sK = reshape(KE0(:)*(Emin + ePhiProj'.^penal*(E0-Emin)), 64*tolne, 1);
    K = sparse(iK,jK,sK); K=(K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    
    %% 目标函数和灵敏度评估 %%
    obj = F'*U;
    ce = sum((U(edofMat)*KE0).*U(edofMat), 2);
    dcdx = eIntopMat*(-penal*(E0-Emin)*ePhiProj(:).^(penal-1).*ce.*edproj(:));
    vol = sum(ePhiProj(:)*elelen^2);
    voldgdx = eIntopMat*(edproj(:)*elelen^2);
    
    %% 更新设计变量 %%
    m = 1; cc = 10000*ones(m,1); d = zeros(m,1); a0 = 1; a = zeros(m,1);
    fval = zeros(m, 1); dfdx = zeros(m, n);
    fval(1) = 100*(vol/tolvol-volfra);
    dfdx(1,:) = 100*voldgdx/tolvol;
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp]=...
        mmasub(m,n,loop,x,xmin,xmax,xold1,xold2, ...
        obj,dcdx,fval,dfdx,low,upp,a0,a,cc,d);
    xold2 = xold1; xold1 = x; x = xmma;
    
    %% 投影参数更新 %%
    change = abs(obj-objold)/obj;
    if change < 0.005 && loop > 30
        ichange = ichange+1;
    else
        ichange = 1;
    end
    if mod(ichange,3) == 0
        beta = min(beta * 1.1,20);
    end

    %% 输出结果并绘制密度图 %%
    fprintf([' It.:%5i Obj.:%9.4f Vol:%7.4f numdesvars :%5i' ...
        ' beta:%5.1f ch.:%6.3f\n'],...
        loop,obj,vol/tolvol,neig,beta,change);  
    figure(1); clf;
    displayx = zeros(nely, 2*nelx);
    displayx(:, 1:nelx) = flip(reshape(ePhiProj, nely, nelx),2);
    displayx(:, nelx+1:end) = displayx(:, nelx:-1:1);
    colormap(gray); clims=[-1 0];imagesc(-displayx,clims); 
    axis equal; axis tight;
    title('Elemental density distribution');
    set(gca,'XTick',[0 1e5]);set(gca,'YTick',[0 1e5]);
    figure(2); clf;
    ePhi = reshape(ePhi,size(Xe));
    contourf([Xe],[Ye],[fliplr(ePhi)],[0 0]);
    colormap([0 0 0; 0 0 1;1 0 0; 0 1 0; 1 1 1]);
    axis equal; axis tight; set(gca,'XTick',[]); set(gca,'YTick',[]);
    figure(3); clf;
    Obj = cat(2,Obj,obj); Volf = cat(2,Volf,vol/tolvol);
    plotConvergence(Obj,Volf);
end
%% 保存模型信息
save('result.mat','ePhiProj');
end
