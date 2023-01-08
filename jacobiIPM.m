function [Jf,Jg,Jh] = jacobiIPM(bus,gen,trsfm,shtcr,x,e,f,e_old,f_old,GYbus,BYbus,GYbus_non,BYbus_non,Ybus_non,tbnode,tbgt,tbbt)

% reference bus index
refgen = find(bus(:, 2) == 3);
ref = bus(refgen,1);
% PV bus indices
pvgen = find(bus(:, 2) == 2);
pv = bus(pvgen,1);
 if ~isempty(trsfm)
ntb = size(find(trsfm(:,10)==1),1);
 else
     ntb = 0;
 end
nb = size(bus,1) - ntb;                % 原始系统节点个数

npvg = size(gen,1); 
nb_add = nb + ntb;                     %  新模型的节点个数 

ncr = size(shtcr,1);                   %  shtcr矩阵的变量个数

n = nb+ntb+nb-1+ntb+ntb+npvg+ncr;      %  变量个数

% -------------- 不等式约束的雅克比矩阵 ------------- %
% 雅克比矩阵形成思路：先按含虚拟节点的全变量分析问题，再在后面去除虚拟节点部分和进行节点移号
% 行（不等式排序数值）- 列 - 值
Jh_1 = [[1:nb]',[1:nb]',2*e_old
       [1:nb]',[1:nb]'+nb_add,2*f_old
        nb+[1:ntb+npvg+ncr]',2*nb_add+[1:ntb+npvg+ncr]',ones(ntb+npvg+ncr,1)
    ];
% ------------------ END ----------------- %

% -------------- 相角等式等式约束的雅克比矩阵 ------------- %
if ~isempty(trsfm)
Jg_1 = [[1:ntb]',tbnode(:,1),f(tbnode(:,2))
    [1:ntb]',tbnode(:,2)+nb_add,e(tbnode(:,1))
    [1:ntb]',tbnode(:,2),-f(tbnode(:,1))
    [1:ntb]',tbnode(:,1)+nb_add,-e(tbnode(:,2))
    [ntb+1:ntb+ntb]',tbnode(:,1),ones(ntb,1)
    [ntb+1:ntb+ntb]',tbnode(:,2),-x((nb+ntb)*2:(nb+ntb)*2+ntb-1)
    [ntb+1:ntb+ntb]',[nb_add*2+1:nb_add*2+ntb]',-e(tbnode(:,2))
    ];
Jg_1(:,1) = Jg_1(:,1) + 2*nb;  % 还是调整顺序，和原始问题的序列方式一致
else
    Jg_1=[];
end
% ------------------ END ----------------- %

% -------------- 有功无功的等式约束的雅克比矩阵 ------------- %
if ~isempty(trsfm)
    Jfg_P = [[1:nb]',[1:nb]',GYbus*e-BYbus*f
        [1:nb]',[1:nb]'+nb_add,BYbus*e+GYbus*f
        Ybus_non(:,1),Ybus_non(:,2),GYbus_non(:,3).*e(GYbus_non(:,1))+BYbus_non(:,3).*f(GYbus_non(:,1))
        Ybus_non(:,1),Ybus_non(:,2)+nb_add,-BYbus_non(:,3).*e(BYbus_non(:,1))+GYbus_non(:,3).*f(GYbus_non(:,1))
        tbnode(:,1),tbnode(:,2),tbgt.*(2*e(tbnode(:,2),1)-e(tbnode(:,3),1))+tbbt.*f(tbnode(:,3),1)
        tbnode(:,1),tbnode(:,3),-tbgt.*e(tbnode(:,2),1)-tbbt.*f(tbnode(:,2),1)
        tbnode(:,1),tbnode(:,2)+nb_add,tbgt.*(2*f(tbnode(:,2),1)-f(tbnode(:,3),1))-tbbt.*e(tbnode(:,3),1)
        tbnode(:,1),tbnode(:,3)+nb_add,-tbgt.*f(tbnode(:,2),1)+tbbt.*e(tbnode(:,2),1)
        ];
else
    Jfg_P = [[1:nb]',[1:nb]',GYbus*e-BYbus*f
        [1:nb]',[1:nb]'+nb_add,BYbus*e+GYbus*f
        Ybus_non(:,1),Ybus_non(:,2),GYbus_non(:,3).*e(GYbus_non(:,1))+BYbus_non(:,3).*f(GYbus_non(:,1))
        Ybus_non(:,1),Ybus_non(:,2)+nb_add,-BYbus_non(:,3).*e(BYbus_non(:,1))+GYbus_non(:,3).*f(GYbus_non(:,1))
        ];
end
% tbbt - tbbt_1
% Ybus_non - Ybus_non_1
% tbgt-tbgt_1
if ~isempty(trsfm)
Jg_Q = [[1:nb]',[1:nb]',-GYbus*f-BYbus*e
    [1:nb]',[1:nb]'+nb_add,GYbus*e-BYbus*f
    Ybus_non(:,1),Ybus_non(:,2),-BYbus_non(:,3).*e(BYbus_non(:,1))+GYbus_non(:,3).*f(BYbus_non(:,1))
    Ybus_non(:,1),Ybus_non(:,2)+nb_add,-GYbus_non(:,3).*e(GYbus_non(:,1))-BYbus_non(:,3).*f(BYbus_non(:,1))
    tbnode(:,1),tbnode(:,2),tbbt.*(-2*e(tbnode(:,2),1)+e(tbnode(:,3),1))+tbgt.*f(tbnode(:,3),1)
    tbnode(:,1),tbnode(:,3),tbbt.*e(tbnode(:,2),1)-tbgt.*f(tbnode(:,2),1)
    tbnode(:,1),tbnode(:,2)+nb_add,tbbt.*(-2*f(tbnode(:,2),1)+f(tbnode(:,3),1))-tbgt.*e(tbnode(:,3),1)
    tbnode(:,1),tbnode(:,3)+nb_add,tbbt.*f(tbnode(:,2),1)+tbgt.*e(tbnode(:,2),1)
    [ref;pv;shtcr(:,1)],2*nb_add+ntb+[1:npvg+ncr]',-ones(ncr+npvg,1)
    ];
elseif isempty(trsfm) && ~isempty(shtcr) 
    Jg_Q = [[1:nb]',[1:nb]',-GYbus*f-BYbus*e
    [1:nb]',[1:nb]'+nb_add,GYbus*e-BYbus*f
    Ybus_non(:,1),Ybus_non(:,2),-BYbus_non(:,3).*e(BYbus_non(:,1))+GYbus_non(:,3).*f(BYbus_non(:,1))
    Ybus_non(:,1),Ybus_non(:,2)+nb_add,-GYbus_non(:,3).*e(GYbus_non(:,1))-BYbus_non(:,3).*f(BYbus_non(:,1))
    [ref;pv;shtcr(:,1)],2*nb_add+ntb+[1:npvg+ncr]',-ones(ncr+npvg,1)
    ];
else
        Jg_Q = [[1:nb]',[1:nb]',-GYbus*f-BYbus*e
    [1:nb]',[1:nb]'+nb_add,GYbus*e-BYbus*f
    Ybus_non(:,1),Ybus_non(:,2),-BYbus_non(:,3).*e(BYbus_non(:,1))+GYbus_non(:,3).*f(BYbus_non(:,1))
    Ybus_non(:,1),Ybus_non(:,2)+nb_add,-GYbus_non(:,3).*e(GYbus_non(:,1))-BYbus_non(:,3).*f(BYbus_non(:,1))
    [ref;pv],2*nb_add+ntb+[1:npvg]',-ones(npvg,1)
    ];
    
end
Jg_Q(:,1) = Jg_Q(:,1)+nb;
% ------------------ END ----------------- %

% -------------- PV节点的等式约束的雅克比矩阵 ------------- %
idgen = find(gen(:,7)==0);
nidgen = size(idgen,1);
if isempty(find(gen(idgen,1)==ref))
    Jg_PV = [ [1:nidgen]',gen(idgen,1),2.*e(gen(idgen,1),1)
        [1:nidgen]',gen(idgen,1)+nb_add,2.*f(gen(idgen,1),1) ];
else
    Jg_PV = [];
    for i = 1:nidgen
        if gen(idgen(i,1),1)==ref
            temp = [ i,gen(idgen(i,1),1),2.*e(gen(idgen(i,1),1),1)];
            Jg_PV = [Jg_PV; temp];
        else
            temp = [ i,gen(idgen(i,1),1),2.*e(gen(idgen(i,1),1),1)
                i,gen(idgen(i,1),1)+nb_add,2.*f(gen(idgen(i,1),1),1)];
            Jg_PV = [Jg_PV; temp];
        end
    end
end
    
Jg_PV(:,1) = Jg_PV(:,1) + 2*nb + 2*ntb ;  % 还是调整顺序，和原始问题的序列方式一致

% ------------------ END ----------------- %



% 形成目标函数的一阶偏导数Jf
% 等式约束的一阶偏导数Jg（去除平衡节点有功平衡方程的一阶偏导数）
% 不等式约束的一阶偏导数Jh
Jh = Jh_1;                 % 不等式约束群的一阶偏导
Jg = [ Jfg_P; Jg_Q; Jg_1; Jg_PV]; % 等式约束群的一阶偏导



[temp1,temp2] = find(Jg(:,2) == ref+nb_add);
Jg(temp1,:) = [];          % 去掉了等式约束雅克比矩阵的对ref节点电压虚部f的一阶导数
[temp1,temp2] = find(Jh(:,2)==ref+nb_add);
Jh(temp1,:) = [];          % 去掉了不等式约束雅克比矩阵的对ref节点电压虚部f的一阶导数
[temp1,temp2] = find(Jg(:,1)==ref);
Jf = Jg(temp1,:);          % 目标函数的一阶偏导数
[temp1,temp2] = find(Jg(:,1)~=ref);
Jg = Jg(temp1,:);         % 等式约束群中去掉对目标函数的一阶偏导数


[temp1,temp2] = find(Jf(:,2)>ref+nb_add);
Jf(temp1,2) = Jf(temp1,2)-1; % 使得fref后的变量的列数均减去1
[temp1,temp2] = find(Jg(:,1)>ref);
Jg(temp1,1) = Jg(temp1,1)-1; % 使得Pref后的行数均减去1
[temp1,temp2] = find(Jg(:,2)>ref+nb_add); % 使得fref后的变量的列数均减去1
Jg(temp1,2) = Jg(temp1,2)-1;
[temp1,temp2] = find(Jh(:,2)>ref+nb_add);
Jh(temp1,2) = Jh(temp1,2)-1;

Jf = sparse(Jf(:,2),1,Jf(:,3),n,1);
Jg = sparse(Jg(:,1),Jg(:,2),Jg(:,3),2*nb-1+ntb+ntb+nidgen,n);
Jh = sparse(Jh(:,1),Jh(:,2),Jh(:,3),nb+ntb+npvg+ncr,n);

% ------------------ END ----------------- %