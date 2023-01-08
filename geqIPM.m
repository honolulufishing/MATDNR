
function g = geqIPM(bus,gen,trsfm,shtcr,x,GYbus,BYbus,e,f,e_old,f_old,tbnode,tbgt,tbbt)


% reference bus index
refgen = find(bus(:, 2) == 3);
ref = bus(refgen,1);
% PV bus indices
pvgen = find(bus(:, 2) == 2);
pv = bus(pvgen,1);
% PQ bus indices
Nbus([ref;pv]) = 0;
pq = find(Nbus);
if ~isempty(trsfm)
ntb = size(find(trsfm(:,10)==1),1);
else
    ntb =0;
end
nb = size(bus,1) - ntb;                % 原始系统节点个数
npvg = size(gen,1);
nb_add = nb + ntb;                     %  新模型的节点个数

ncr = size(shtcr,1);                   %  shtcr矩阵的变量个数

n = nb+ntb+nb-1+ntb+ntb+npvg+ncr;      %  变量个数


% % % ------------- 潮流功率平衡方程计算方法(1) -------------- %
% Pbus = -bus(1:nb,3);
% Qbus = -bus(1:nb,4);
% % 节点注入有功、无功(记入了各个节点的负荷有功、无功，计及了发电机的注入有功，所以Pbus,Qbus是不变的)
% Pbus(gen(:,1),1) = Pbus(gen(:,1),1) + gen(:,2);  % 把发电机中的有功注入增加到Pbus中来
% 
%     %----------计算等式约束函数值g-----------%
%     g = zeros(2*nb+ntb+ntb,1);%对等式约束函数值进行初始化
%     % 计算各个节点注入有功等式约束的函数值
%     gP = (GYbus*e-BYbus*f).*e_old+(BYbus*e+GYbus*f).*f_old - Pbus;
%     % 计算各个节点注入无功等式约束的函数值
%     gQ = (GYbus*e-BYbus*f).*f_old-(BYbus*e+GYbus*f).*e_old-Qbus;
% 
%     gQ([ref;pv],1) = gQ([ref;pv],1) - x(2*(nb+ntb)+ntb:2*(nb+ntb)-1+ntb+npvg,1); % 用发电机注入无功来更新gQ
%     gQ(shtcr(:,1),1) = gQ(shtcr(:,1),1) - x(2*(nb+ntb)+ntb+npvg:n,1);
% 
%     % 用电容电抗注入无功来更新gQ
%     % 计算各个节点注入无功等式约束的函数值
%     % 以下是对原可调变压器高压端节点的注入有、无功进行修正
%     eii = e(tbnode(:,2),1);fii=f(tbnode(:,2),1);
%     ej = e(tbnode(:,3),1);fj=f(tbnode(:,3),1);
%     % 分别给引入节点电压的实虚部、原可调变压器低压端的节点电压实虚部赋值
%     tbgP = (eii.^2+fii.^2).*tbgt-(eii.*ej+fii.*fj).*tbgt+(eii.*fj-ej.*fii).*tbbt;%计算可调变压器原高压端节点注入有功应该增加的部分
%     tbgQ = -(eii.^2+fii.^2).*tbbt+(eii.*ej+fii.*fj).*tbbt+(eii.*fj-ej.*fii).*tbgt;%计算可调变压器原高压端节点注入无功应该增加的部分
%   for m = 1:ntb
%         t = tbnode(m,1);%取出可调变压器原高压端节点
%         gP(t,1) = gP(t,1) + tbgP(m,1);%对原可调变压器高压端节点的注入有功进行修正
%         gQ(t,1) = gQ(t,1) + tbgQ(m,1);
%   end
%  % ------------- END -------------- %

% % ------------- 潮流功率平衡方程计算方法二 -------------- %
Pbus = bus(1:nb,3);
Qbus = bus(1:nb,4);
% 节点注入有功、无功(记入了各个节点的负荷有功、无功，计及了发电机的注入有功，所以Pbus,Qbus是不变的)
Pbus(gen(:,1),1) = Pbus(gen(:,1),1) + gen(:,2);  % 把发电机中的有功注入增加到Pbus中来

% g = zeros(2*nb+ntb+ntb,1);   % 对等式约束函数值进行初始化
j = sqrt(-1);

U =  sqrt(e.^2+f.^2).*(cos( atan(f./e))+j*sin( atan(f./e)));
Sbus = Pbus + j*Qbus ;

Sbus([ref;pv],:) = Sbus([ref;pv],:) + j*x(2*(nb+ntb)+ntb:2*(nb+ntb)-1+ntb+npvg,1); % 用发电机注入无功来更新gQ
if ~isempty(shtcr)
Sbus(shtcr(:,1),:) = Sbus(shtcr(:,1),:)  + j*x(2*(nb+ntb)+ntb+npvg:n,1);
end

Y = GYbus+j.*BYbus;

mis = U(1:nb,1) .* conj(Y * U) - Sbus;

delta_P = real(mis(1:nb,:));
delta_Q = imag(mis(1:nb,:));

gP = delta_P;
gQ = delta_Q;

    % 用电容电抗注入无功来更新gQ
    % 计算各个节点注入无功等式约束的函数值
    
    if ~isempty(trsfm)
    % 以下是对原可调变压器高压端节点的注入有、无功进行修正
    eii = e(tbnode(:,2),1);fii=f(tbnode(:,2),1);
    ej = e(tbnode(:,3),1);fj=f(tbnode(:,3),1);
    % 分别给引入节点电压的实虚部、原可调变压器低压端的节点电压实虚部赋值
    tbgP = (eii.^2+fii.^2).*tbgt-(eii.*ej+fii.*fj).*tbgt+(eii.*fj-ej.*fii).*tbbt;%计算可调变压器原高压端节点注入有功应该增加的部分
    tbgQ = -(eii.^2+fii.^2).*tbbt+(eii.*ej+fii.*fj).*tbbt+(eii.*fj-ej.*fii).*tbgt;%计算可调变压器原高压端节点注入无功应该增加的部分
  for m = 1:ntb
        t = tbnode(m,1);%取出可调变压器原高压端节点
        gP(t,1) = gP(t,1) + tbgP(m,1);%对原可调变压器高压端节点的注入有功进行修正
        gQ(t,1) = gQ(t,1) + tbgQ(m,1);
  end
 
% ------------- END -------------- %


% 对原可调变压器高压端节点的注入有、无功进行修正
eii = e(tbnode(:,2),1);fii=f(tbnode(:,2),1);
ei = e(tbnode(:,1),1);fi=f(tbnode(:,1),1);%分别给原可调变压器高压端的节点电压实虚部赋值
ki = x((nb+ntb)*2:(nb+ntb)*2+ntb-1,1);
g_theta = ei.*fii-eii.*fi;%计算可调变压器原高压端节点和引入节点之间的电压相角相等的等式约束的函数值
gK = ei-ki.*eii;
    else
        gK=[];g_theta=[];
   end

idgen = find(gen(:,7)==0);
gPV = zeros(size(idgen,1),1);
for i = 1: size(idgen,1)
    
    if gen(idgen(i,1),1) > ref
        gPV(i,1) = x(gen(idgen(i,1),1),1).^2  +  x(gen(idgen(i,1),1)+nb+ntb-1,1).^2 -  gen(idgen(i,1),6).^2;
    elseif gen(idgen(i,1),1) < ref
        gPV(i,1) = x(gen(idgen(i,1),1),1).^2  +  x(gen(idgen(i,1),1)+nb+ntb,1).^2 -  gen(idgen(i,1),6).^2;
    elseif gen(idgen(i,1),1) == ref
        gPV(i,1) = x(gen(idgen(i,1),1),1).^2   -  gen(idgen(i,1),6).^2;
    end
end

g = [gP;gQ;g_theta;gK; gPV ];