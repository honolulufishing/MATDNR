function [bus,gen,trsfm,shtcr,loss,elapsed_time,Record_KKT,success] = VARIPM(bus,gen,Lnbr,trsfm,shtc,shtr,vctr)


% ---------------- 发电机功率数据初始化 -------------- %
bus_s = bus;
% reference bus index
refgen = find(bus(:, 2) == 3);
ref = bus(refgen,1);
% PV bus indices
pvgen = find(bus(:, 2) == 2);
pv = bus(pvgen,1);
gen(:,3) = 0;
gen(find(gen(:,1)==ref),2) = 0;
% -------------------- END -------------------- %

% -------------- 变压器变比数据初始化 ------------- %
j = sqrt(-1);
if ~isempty(trsfm)
    temp = find(trsfm(:,10) == 1 );
    trsfm(temp,6) = 1;  % 将变压器的变比初始化值定义为1，以后优化中会发生变化   
    ntb = size(find(trsfm(:,10)==1),1);
else
     ntb =0;
end
% --------------------- END -------------------- %

% 生成虚拟节点
[bus,trsfm,shtcr,tbnode,tbgt,tbbt,tbkmax,tbkmin] = initialIPM(bus,gen,Lnbr,trsfm,shtc,shtr);

% 生成海森矩阵
[save_H,GYbus,BYbus ,GYbus_non,BYbus_non,Ybus_non] = hessianIPM(bus,gen,Lnbr,trsfm,shtc,shtr,tbnode,tbgt,tbbt);

% -------------- 优化变量个数信息初始化 ------------- %
nb = size(bus_s,1);

refgen = find(bus(:, 2) == 3);
ref = bus(refgen,1);                   % reference bus index
npvg = size(gen,1);                    % PV加slack节点
ncr = size(shtcr,1);                   %  shtcr矩阵的变量个数
n = nb+ntb+nb-1+ntb+ntb+npvg+ncr;      %  变量个数
nb_add = nb + ntb;                     %  新模型的节点个数
% --------------------- END -------------------- %

% -------------- 不等式约束的上下限lb和ub ------------- %
% 不等式约束的上下限
if ~isempty(trsfm)
    lb = zeros(nb+ntb+npvg+ncr,1);   % 不等式下限约束
    ub = zeros(nb+ntb+npvg+ncr,1);   % 不等式上限约束
    lb(1:nb,1) = bus(1:nb,9).^2;     % 电压幅值下限
    ub(1:nb,1) = bus(1:nb,10).^2;    % 电压幅值上限
    lb(nb+1:nb+ntb,1) = tbkmin;      % 变比下限值约束
    ub(nb+1:nb+ntb,1) = tbkmax;      % 变比上限值约束
else
    lb = zeros(nb+npvg,1);   % 不等式下限约束
    ub = zeros(nb+npvg,1);   % 不等式上限约束
    lb(1:nb,1) = bus(1:nb,9).^2;     % 电压幅值下限
    ub(1:nb,1) = bus(1:nb,10).^2;    % 电压幅值上限
end

temp_ref = find(gen(:,1)==ref);
temp_pv = setdiff(1:size(gen,1),temp_ref)';
temp =[ temp_ref; temp_pv ];   % 将平衡节点放在gen第1位，发电机节点中没有PQ节点

lb(nb+ntb+1:nb+ntb+npvg,1) = gen(temp,5); % 发电机无功出力下限
ub(nb+ntb+1:nb+ntb+npvg,1) = gen(temp,4); % 发电机无功出力上限

if ~isempty(shtcr)    
lb(nb+ntb+npvg+1:nb+ntb+npvg+ncr,1) = shtcr(:,4); % 无功补偿下限
ub(nb+ntb+npvg+1:nb+ntb+npvg+ncr,1) = shtcr(:,3); % 无功补偿上限
end


% --------------------- END -------------------- %

% -------------- 内点法求解初始迭代点x0 ------------- %
% x初始化：x=[e f k QG Qcri]
x = zeros(n,1);               % 迭代变量

idgen = find(gen(:,7)==0);
nidgen = size(idgen,1);

y = zeros(2*nb-1+ntb+ntb+nidgen,1);  % 等式约束的拉格朗日向量
h = zeros(nb+ntb+npvg+ncr,1); % 不等式约束函数向量
x(1:nb,1) = (bus(1:nb,9)+bus(1:nb,10))/2;  % 普通节点电压实部的初始值
if ~isempty(trsfm)
    x(nb+1:nb+ntb,1) = x(tbnode(:,1),1);            % 引入的节点电压实部的初始值
end
x(nb+ntb+1:2*(nb+ntb)-1,1) = 0;                 % 所有节点（除了平衡节点的外）的电压相角的初值均为0，所以对应虚部也为0
x(2*(nb+ntb):n,1) = (lb(nb+1:nb+ntb+npvg+ncr,1)+ub(nb+1:nb+ntb+npvg+ncr,1))/2; % 发电机注入无功、补偿无功的的初始值

% ------------------ END ----------------- %

% -------------- 内点法求解对偶变量y,l,u,z,w的初始点 ------------- %
y(1:nb-1,1) = -10;    % y 拉格朗日对偶变量（对应功率等式约束）
y(2*nb-1+ntb+ntb+1:2*nb-1+ntb+ntb+nidgen) = 10; % y 拉格朗日对偶变量（对应PV节点等式约束）
e = x(1:nb+ntb,1);                                               % 所有节点的实部值
f = [x(nb+ntb+1:nb+ntb+ref-1,1);0;x(nb+ntb+ref:2*nb+2*ntb-1,1)]; % 所有节点的虚部值
e_old = x(1:nb,1);                                                 % 所有普通节点的实部值
f_old = [x(nb+ntb+1:nb+ntb+ref-1,1);0;x(nb+ntb+ref:2*nb+ntb-1,1)]; % 所有普通节点的虚部值
% 不等式处h(x)初值求取
h = [ e_old.^2+f_old.^2 ; x(2*(nb+ntb):n,1)];
% 计算不等式约束函数值 h
l = h-lb;
u = ub-h;
z = 10*ones(nb+ntb+npvg+ncr,1);
w = z;                % w,z 拉格朗日对偶变量（对应不等式约束）
% ------------------ END ----------------- %

% -------------- 内点法生成雅克比矩阵 ------------- %
% [Jf,Jg,Jh] = jacobiIPM(bus,gen,trsfm,shtcr,x,e,f,e_old,f_old,GYbus,BYbus,GYbus_non,BYbus_non,Ybus_non,tbnode,tbgt,tbbt);
% ------------------ END ----------------- %

% -------------- 设置内点法求解的参数信息 ------------- %
Gap = l'*z+u'*w;   % 计算互补间隙的初始值
MaxIter = 200;    % 最大迭代次数
tol_1 = 10^(-6);   % 互补间隙收敛精度
tol_2 = 10^(-6);   % KKT条件最大范数收敛精度
% ------------------ END ----------------- %

C = sparse(2*nb-1+ntb+ntb+nidgen,2*nb-1+ntb+ntb+nidgen);
P = [];

success = 0;

Record_goalf = zeros(MaxIter,1);
Record_KKT = zeros(MaxIter,4);
Record_Gap = zeros(MaxIter,1); % Complementary Gap of PCPDIPM

ii = 1;
tStart = cputime;
while ii<= MaxIter

    % 计算等式约束函数值g
    g = geqIPM(bus,gen,trsfm,shtcr,x,GYbus,BYbus,e,f,e_old,f_old,tbnode,tbgt,tbbt);
    
    % 计算目标函数goalf
    goalf = g(ref,:);    % 目标函数的值
    g(ref,:) = [];       % 删掉平衡节点有功注入所对应的行
    
    
    % -------------- 内点法生成雅克比矩阵 ------------- %
    [Jf,Jg,Jh] = jacobiIPM(bus,gen,trsfm,shtcr,x,e,f,e_old,f_old,GYbus,BYbus,GYbus_non,BYbus_non,Ybus_non,tbnode,tbgt,tbbt);
    % ------------------ END ----------------- %
    
    
    % --------- 计算KKT方程 -------- %
    Lx = Jf - Jg'*y - Jh'*(z-w);  % 计算Lx
    Ly = g;                       % 计算Ly
    Lz = h - l - lb;              % 计算Lz
    Lw = h + u - ub;              % 计算Lw
    % ------------- END ------------- %
    
    % 记录第ii次迭代的目标函数值
    Record_goalf(ii,1) = goalf;
    % 记录第ii次迭代的KKT方程不平衡量
    Record_KKT(ii,:) = [max(abs(Lx)) max(abs(Ly))  max(abs(Lz))  max(abs(Lw))];
    Record_Gap(ii,:) = Gap;
    
    % Gap与KT条件的收敛判据
    if Gap < tol_1                                      % Gap的收敛判据
        if max([abs(Lx);abs(Ly);abs(Lz);abs(Lw)])<tol_2 % KKT的收敛判据
            
            success = 1;
            
            % 有功网损
            j = sqrt(-1);
            gen(find(gen(:,1)==ref),2) = 0;
            loss = goalf + sum(gen(:,2),1) + sum(bus(:,3),1);
            % 电压幅值和相角
            V = abs(e+j*f);
            bus(:,7) =  V;
            bus(:,8) = angle(e+j.*f);
            % 发电机机端电压和无功功率
            gen(:,6) = bus(gen(:,1),7);
            gen(temp,3) = x(2*(nb+ntb)+ntb:2*(nb+ntb)+ntb-1+npvg,1);
            gen(temp_ref,2) = goalf;
            if ~isempty(shtc)
            % 电容电抗的注入无功
            idshtc = find(shtc(:,6)==1);
            shtcr(idshtc ,2) = x(2*(nb+ntb)+ntb+npvg:n,1);
            end
            if ~isempty(trsfm)
            % 变压器变比
            idtrsfm = find(trsfm(:,10)==1);
            trsfm(idtrsfm,6) = x(2*(nb+ntb):2*(nb+ntb)+ntb-1,1);
            end
            
            % ----------------------- Output Solution ------------------%
            
            disp('   ')
            fprintf('Success:Optimal Power Flow using IPM converged in  %g iterations.\n', ii)
            fprintf('The System Total Real Loss = %g MW\n\n', loss)

%             % Output Bus Solution
%             disp('                                                                            ')
%             disp('   Optimal Power Flow Solution by IPM With a Second-order Reactive Optimization Model')
%             fprintf('                         No. of Iterations = %g \n', ii)
%             fprintf('                 The System Total Real Loss = %g MW\n\n', loss)
%             head =['    Bus  Bus    ------Load------    ---Injected---     Voltage    Angle  '
%                 '    No.  Type    MW       Mvar       GS       BS         Mag.      Deg.  '
%                 '                                                                         '];
%             disp(head)
%             for i = 1:size(bus,1) - size(trsfm,1)
%                 fprintf(' %5g', bus(i,1)), fprintf(' %5g', bus(i,2)),
%                 fprintf(' %8.3f', bus(i,3)), fprintf(' %9.3f', bus(i,4)),
%                 fprintf(' %9.3f', bus(i,5)),  fprintf(' %9.3f', bus(i,6)),
%                 fprintf(' %9.3f ', bus(i,7)), fprintf(' %8.3f\n',bus(i,8)*180/pi)
%             end
%             fprintf('      \n'), fprintf('  Total     ')
%             fprintf(' %8.3f', sum(bus(:,3))), fprintf(' %9.3f\n\n', sum(bus(:,4)))
%             
%             
%             % Output Generator Solution
%             disp('                                     ')
%             disp('                   Generator Solution')
%             
%             head =['                                                     '
%                 '    Bus     ----Gen----       -----Q-----       Vol  '
%                 '    No.     MW      MVA       Max     Min       Mag. '
%                 '                                                     '];
%             disp(head)
%             for i = 1:size(gen,1)
%                 fprintf(' %5g', gen(i,1)), fprintf(' %8.3f', gen(i,2)),
%                 fprintf(' %8.3f', gen(i,3)), fprintf(' %9.3f',gen(i,4)),
%                 fprintf(' %9.3f', gen(i,5)),  fprintf(' %8.3f\n', gen(i,6)),
%             end
%             fprintf('      \n'), fprintf('  Total ')
%             fprintf(' %6.3f', sum(gen(:,2))), fprintf(' %8.3f\n\n', sum(gen(:,3)))
%             
%             
%             % Output Transformer Branch Solution
%             disp('                                       ')
%             disp('                   Transformer Solution')
%             
%             head =['                                                                        '
%                 '   TrsBR.  Bus   Bus     R      X     ----- Tap -----   Series   State  '
%                 '    NO.    from   to    p.u.   p.u.   Ratio  Max  Min      #       #    '
%                 '                                                                        '];
%             disp(head)
%             
%             for i = 1:size(trsfm,1)
%                 fprintf(' %5g', trsfm(i,1)), fprintf(' %5.0f', trsfm(i,2)),
%                 fprintf(' %5g', trsfm(i,3)), fprintf(' %8.3f', trsfm(i,4)),
%                 fprintf(' %6.3f', trsfm(i,5)),fprintf(' %5.5f',trsfm(i,6)),
%                 fprintf(' %5.1f', trsfm(i,7)),fprintf(' %5.1f',trsfm(i,8)),
%                 fprintf(' %5g',trsfm(i,9)),fprintf(' %5g\n', trsfm(i,10)),
%              end
%            
%             % Output shtcr Solution
%             disp('                                       ')
%             disp('                   Shunt   Solution')
%             head =['                                 '
%                    '   Bus.   shtcr    Max     Min   '
%                    '    NO.    Mvar    Mvar    Mvar  '
%                    '                                 '];
%             disp(head)
%             
%             for i = 1:size(shtcr,1)
%                 fprintf(' %5g',shtcr(i,1)), fprintf('   %5.3f', shtcr(i,2)),
%                 fprintf('   %5.1f', shtcr(i,3)), fprintf('   %5.1f\n', shtcr(i,4)),
%             end
%             fprintf('      \n'), fprintf('  Total ')
%             fprintf(' %6.3f\n', sum(shtcr(:,2)))
            
            % ----------------------- END ------------------%
            
            break;
        end
    end
 
    % --------------- 求解仿射方程 ------------------ %
    % 把目标函数、等式约束、不等式约束的海森矩阵的拉咯朗系数合并在一起
    yy = [-y(1:ref-1,1);1;-y(ref:2*nb-1+ntb+ntb+nidgen,1);-(z-w)];
    H = sparse(save_H(:,1),save_H(:,2),save_H(:,3).*yy(save_H(:,4)),n,n);
    ln1 = reciprocal(l);
    un1 = reciprocal(u);
    vmid = z.*ln1 + w.*un1;
    H = H+Jh'*diag(vmid,0)*Jh;     % 求到了H
  
    vmid = (z.*l+z.*Lz).*ln1-(w.*u-w.*Lw).*un1;
    kesai_aff = Lx + Jh'*vmid;
    % 求到了kesai_aff
    A = [H ,-Jg'
        Jg, C];
    B = [kesai_aff;Ly];
    
    if ii==1||ii==2
        P = symrcm(A);
    end
    [L,U] = lu(A(P,P));
    dxy_aff(P,1) = -U\(L\B(P));
    
    dx_aff = dxy_aff(1:n,1);
    dy_aff = dxy_aff(n+1:n+2*nb-1+ntb+ntb+nidgen,1);
    
    dl_aff = Jh*dx_aff + Lz;
    du_aff = -dl_aff + Lz - Lw;
    dz_aff = -(z.*l + z.*dl_aff).*ln1;
    dw_aff = -(w.*u + w.*du_aff).*un1;
    % ------------END----------- %
    
    % 确定仿射步长：stepp_aff,stepd_aff,互补间隙：Gap_aff
    templ = find(dl_aff<0);tempu = find(du_aff<0);
    stepp_aff = min([-l(templ,1)./dl_aff(templ,1);-u(tempu,1)./du_aff(tempu,1);1]);
    tempz = find(dz_aff<0);tempw=find(dw_aff<0);
    stepd_aff = min([-z(tempz,1)./dz_aff(tempz,1);-w(tempw,1)./dw_aff(tempw,1);1]);
    
    if isempty(stepp_aff)
        stepp_aff=1;
    end
    if isempty(stepd_aff)
        stepd_aff=1;
    end
    Gap_aff=(z+stepd_aff*dz_aff)'*(l+stepp_aff*dl_aff)+(w+stepd_aff*dw_aff)'*(u+stepp_aff*du_aff);
    % 中心参数sigema,预测障碍因子mu
    %     sigema=(Gap_aff/Gap)^3;
    %     mu=sigema*Gap/2/(nb+2*ntb+npvg+ncr);
    mu=min((Gap_aff/Gap)^2,0.2)*Gap_aff/2/(nb+ntb+npvg+ncr);
    
    % -----------间接求解校正方程------------- %
    vmid=(dl_aff.*dz_aff-mu).*ln1-(dw_aff.*du_aff-mu).*un1;
    kesai_cor=Jh'*vmid;
    B = [kesai_cor;sparse(2*nb-1+ntb+ntb+nidgen,1)];
  
    dxy_cor(P,1)=-U\(L\full(B(P)));
    dx_cor=dxy_cor(1:n,1);
    dy_cor=dxy_cor(n+1:n+2*nb-1+ntb+ntb+nidgen,1);
    dl_cor=Jh*dx_cor;
    du_cor=-dl_cor;
    dz_cor=(mu-dl_aff.*dz_aff-z.*dl_cor).*ln1;
    dw_cor=(mu-du_aff.*dw_aff-w.*du_cor).*un1;
    
    dx=dx_cor+dx_aff;
    dy=dy_cor+dy_aff;
    dl=dl_cor+dl_aff;
    du=du_cor+du_aff;
    dz=dz_cor+dz_aff;
    dw=dw_cor+dw_aff;
    % ------------END----------- %
   
    % ----------直接求解校正方程---------- %
    %     vmid=(z.*l-mu+dz_aff.*dl_aff+z.*Lz).*ln1-(w.*u-mu+dw_aff.*du_aff-w.*Lw).*un1;
    %     kesai=Lx+Jh'*vmid;
    %     B=[kesai;Ly];
    %     dxy=-(A\B);
    %     dx=dxy(1:n,1);
    %     dy=dxy(n+1:n+2*nb-1+ntb,1);
    %     dl=Jh*dx+Lz;
    %     du=-dl+Lz-Lw;
    %     dz=-(z.*l-mu+dl_aff.*dz_aff+z.*dl).*ln1;
    %     dw=-(w.*u-mu+du_aff.*dw_aff+w.*du).*un1;
    % ------------END----------- %
    
    % 确定迭代步长：stepp,stepd
    templ=find(dl<0);
    tempu=find(du<0);
    stepp=0.9995*min([-l(templ,1)./dl(templ,1);-u(tempu,1)./du(tempu,1);1]);
    tempz=find(dz<0);tempw=find(dw<0);
    stepd=0.9995*min([-z(tempz,1)./dz(tempz,1);-w(tempw,1)./dw(tempw,1);1]);
    if isempty(stepp)
        stepp=1;
    end
    if isempty(stepd)
        stepd=1;
    end
    
    % 更新原变量和对偶变量x,y,l,u,z,w
    x = x + stepp*dx;
    y = y + stepd*dy;
    l = l + stepp*dl;
    u = u + stepp*du;
    z = z + stepd*dz;
    w = w + stepd*dw;
    
    % 得到了新一次的变量值x,y,l,u,z,w
    e = x(1:nb+ntb,1);                            % 所有节点的实部值
    f = [x(nb+ntb+1:nb+ntb+ref-1,1);0;x(nb+ntb+ref:2*nb+2*ntb-1,1)];% 所有节点的虚部值
    e_old = x(1:nb,1);                            % 所有普通节点的实部值
    f_old = [x(nb+ntb+1:nb+ntb+ref-1,1);0;x(nb+ntb+ref:2*nb+ntb-1,1)];% 所有普通节点的虚部值
    h = [e_old.^2+f_old.^2;x(2*(nb+ntb):n,1)];    % 计算不等式约束函数值

    Gap = l'*z+u'*w;     % 更新Gap

    % -------------- 内点法生成雅克比矩阵 ------------- %
    % [Jf,Jg,Jh] = jacobiIPM(bus,gen,trsfm,shtcr,x,e,f,e_old,f_old,GYbus,BYbus,GYbus_non,BYbus_non,Ybus_non,tbnode,tbgt,tbbt);
    % ------------------ END ----------------- %

    if ii == MaxIter%最大迭代次数的收敛判据
        disp('Warning:内点法达到了最大循环次数，无法收敛！');
        break;
    end
    
    ii = ii + 1;    % 迭代次数ii增加1
end

elapsed_time =  cputime- tStart;



% 虚拟节点剔除
bus = bus(1:nb,:);
% % 便于用于潮流计算
% bus (:,3:6) = bus (:,3:6);
if ~isempty(trsfm)
idtrsfm = find(trsfm(:,10)==1);
trsfm(idtrsfm,2) = tbnode(:,1);
end


% Figure
% figure(1)
% xk = 1:ii;
% plot(xk,Record_Gap(1:ii,1),'ro-')
% xlabel('Iteration','fontsize',10,'fontweight','b');
% ylabel('Gap','fontsize',10,'fontweight','b');
% title('Complementary Gap of PCPDIPM','fontsize',10,'fontweight','b')

