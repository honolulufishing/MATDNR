function [save_H,GYbus,BYbus ,GYbus_non,BYbus_non,Ybus_non] = hessianIPM(bus,gen,Lnbr,trsfm,shtc,shtr,tbnode,tbgt,tbbt)

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
nb = size(bus,1) - ntb;
% -------------- 优化变量个数信息初始化 ------------- %
%ncr = size(shtcr,1);                 %  shtcr矩阵的变量个数
%n = nb+ntb+nb-1+ntb+ntb+npvg+ncr;    %  变量个数 
nb_add = nb + ntb;                     %  新模型的节点个数 
% --------------------- END -------------------- %

% -------------- 生成新模型下的节点导纳矩阵Ybus------------- %
[Ybus,Yf,Yt] = makeYbus(bus,Lnbr,trsfm,shtc,shtr,pq,pv);

GYbus = real(Ybus(1:nb,:));   % Ybus实部: nb×(nb+ntb)维矩阵
BYbus = imag(Ybus(1:nb,:));   % Ybus虚部: nb×(nb+ntb)维矩阵

% 先寻找Ybus中的非零元素  begin
GYbus_non = [];
BYbus_non = [];
% GYbus_non=[非零元素的行号，非零元素的列号，非零元素];
% BYbus_non=[非零元素的行号，非零元素的列号，非零元素];
for t = 1: nb
    [row,column,non] = find(Ybus(t,:));
    GYbus_non = [ GYbus_non;[t*row',column',non.']];
end
Ybus_non = GYbus_non;
BYbus_non =[ GYbus_non(:,1:2),imag(GYbus_non(:,3)) ]; % 找出虚部非零元，维数等于nnz(BYbus)
GYbus_non(:,3) = real(GYbus_non(:,3));                % 找出实部非零元，维数等于nnz(GYbus)
% --------------------- END -------------------- %

% ------------------ 生成海森矩阵 ----------------- %
% -------------- 相角相等的等式约束的海森矩阵 ------------- %
if ~isempty(trsfm)
% save=[非零元素的行号，非零元素的列号，非零元素，所属的哪一个海森矩阵]
% H1（ei,fm）=1,H1(em,fi)=-1 ，等式约束处于2*nb+1:2*nb+ntb段 
save_1=[ tbnode(:,1),tbnode(:,2)+nb_add,ones(ntb,1),[2*nb+1:2*nb+ntb]'
    tbnode(:,2)+nb_add,tbnode(:,1),ones(ntb,1),[2*nb+1:2*nb+ntb]'
    tbnode(:,2),tbnode(:,1)+nb_add,-ones(ntb,1),[2*nb+1:2*nb+ntb]'
    tbnode(:,1)+nb_add,tbnode(:,2),-ones(ntb,1),[2*nb+1:2*nb+ntb]'
    tbnode(:,2),nb_add*2+[1:ntb]',-ones(ntb,1),[2*nb+ntb+1:2*nb+ntb+ntb]'
    nb_add*2+[1:ntb]',tbnode(:,2),-ones(ntb,1),[2*nb+ntb+1:2*nb+ntb+ntb]'
    ];
else
    save_1=[];
end
% --------------------- END -------------------- %

% -------------- 节点电压幅值上下限不等式约束的海森矩阵 ------------- %
% H2(ei,ei)=2,H2(fi,fi)=2,i=1,2...nb
save_2=[[1:nb]',[1:nb]',2*ones(nb,1),[2*nb+ntb+ntb+1:3*nb+ntb+ntb]'
    [1:nb]'+nb_add,[1:nb]'+nb_add,2*ones(nb,1),[2*nb+ntb+ntb+1:3*nb+ntb+ntb]'];
% --------------------- END -------------------- %

% -------------- 节点注入有功和无功等式约束的海森矩阵 ------------- %
% save_7为有功功率,save_8为无功功率方程的二次偏导数
if ~isempty(trsfm)
save_7 = [ GYbus_non(:,1),GYbus_non(:,2),GYbus_non(:,3),GYbus_non(:,1)
    GYbus_non(:,2),GYbus_non(:,1),GYbus_non(:,3),GYbus_non(:,1)
    GYbus_non(:,1)+nb_add,GYbus_non(:,2)+nb_add,GYbus_non(:,3),GYbus_non(:,1)
    GYbus_non(:,2)+nb_add,GYbus_non(:,1)+nb_add,GYbus_non(:,3),GYbus_non(:,1)
    BYbus_non(:,1),BYbus_non(:,2)+nb_add,-BYbus_non(:,3),BYbus_non(:,1)
    BYbus_non(:,2)+nb_add,BYbus_non(:,1),-BYbus_non(:,3),BYbus_non(:,1)
    BYbus_non(:,1)+nb_add,BYbus_non(:,2),BYbus_non(:,3),BYbus_non(:,1)
    BYbus_non(:,2),BYbus_non(:,1)+nb_add,BYbus_non(:,3),BYbus_non(:,1)
    tbnode(:,2),tbnode(:,2),2*tbgt,tbnode(:,1)
    tbnode(:,2)+nb_add,tbnode(:,2)+nb_add,2*tbgt,tbnode(:,1)
    tbnode(:,2),tbnode(:,3),-tbgt,tbnode(:,1)
    tbnode(:,3),tbnode(:,2),-tbgt,tbnode(:,1)
    tbnode(:,2),tbnode(:,3)+nb_add,tbbt,tbnode(:,1)
    tbnode(:,3)+nb_add,tbnode(:,2),tbbt,tbnode(:,1)
    tbnode(:,2)+nb_add,tbnode(:,3),-tbbt,tbnode(:,1)
    tbnode(:,3),tbnode(:,2)+nb_add,-tbbt,tbnode(:,1)
    tbnode(:,2)+nb_add,tbnode(:,3)+nb_add,-tbgt,tbnode(:,1)
    tbnode(:,3)+nb_add,tbnode(:,2)+nb_add,-tbgt,tbnode(:,1)
     ];
 
save_8 = [BYbus_non(:,1),BYbus_non(:,2),-BYbus_non(:,3),BYbus_non(:,1)+nb
    BYbus_non(:,2),BYbus_non(:,1),-BYbus_non(:,3),BYbus_non(:,1)+nb
    BYbus_non(:,1)+nb_add,BYbus_non(:,2)+nb_add,-BYbus_non(:,3),BYbus_non(:,1)+nb
    BYbus_non(:,2)+nb_add,BYbus_non(:,1)+nb_add,-BYbus_non(:,3),BYbus_non(:,1)+nb
    GYbus_non(:,1),GYbus_non(:,2)+nb_add,-GYbus_non(:,3),GYbus_non(:,1)+nb
    GYbus_non(:,2)+nb_add,GYbus_non(:,1),-GYbus_non(:,3),GYbus_non(:,1)+nb
    GYbus_non(:,1)+nb_add,GYbus_non(:,2),GYbus_non(:,3),GYbus_non(:,1)+nb
    GYbus_non(:,2),GYbus_non(:,1)+nb_add,GYbus_non(:,3),GYbus_non(:,1)+nb
    tbnode(:,2),tbnode(:,2),-2*tbbt,tbnode(:,1)+nb
    tbnode(:,2)+nb_add,tbnode(:,2)+nb_add,-2*tbbt,tbnode(:,1)+nb
    tbnode(:,2),tbnode(:,3),tbbt,tbnode(:,1)+nb
    tbnode(:,3),tbnode(:,2),tbbt,tbnode(:,1)+nb
    tbnode(:,2),tbnode(:,3)+nb_add,tbgt,tbnode(:,1)+nb
    tbnode(:,3)+nb_add,tbnode(:,2),tbgt,tbnode(:,1)+nb
    tbnode(:,2)+nb_add,tbnode(:,3),-tbgt,tbnode(:,1)+nb
    tbnode(:,3),tbnode(:,2)+nb_add,-tbgt,tbnode(:,1)+nb
    tbnode(:,2)+nb_add,tbnode(:,3)+nb_add,tbbt,tbnode(:,1)+nb
    tbnode(:,3)+nb_add,tbnode(:,2)+nb_add,tbbt,tbnode(:,1)+nb
    ];
else
    save_7 = [ GYbus_non(:,1),GYbus_non(:,2),GYbus_non(:,3),GYbus_non(:,1)
    GYbus_non(:,2),GYbus_non(:,1),GYbus_non(:,3),GYbus_non(:,1)
    GYbus_non(:,1)+nb_add,GYbus_non(:,2)+nb_add,GYbus_non(:,3),GYbus_non(:,1)
    GYbus_non(:,2)+nb_add,GYbus_non(:,1)+nb_add,GYbus_non(:,3),GYbus_non(:,1)
    BYbus_non(:,1),BYbus_non(:,2)+nb_add,-BYbus_non(:,3),BYbus_non(:,1)
    BYbus_non(:,2)+nb_add,BYbus_non(:,1),-BYbus_non(:,3),BYbus_non(:,1)
    BYbus_non(:,1)+nb_add,BYbus_non(:,2),BYbus_non(:,3),BYbus_non(:,1)
    BYbus_non(:,2),BYbus_non(:,1)+nb_add,BYbus_non(:,3),BYbus_non(:,1)
     ];
 
 save_8 = [BYbus_non(:,1),BYbus_non(:,2),-BYbus_non(:,3),BYbus_non(:,1)+nb
    BYbus_non(:,2),BYbus_non(:,1),-BYbus_non(:,3),BYbus_non(:,1)+nb
    BYbus_non(:,1)+nb_add,BYbus_non(:,2)+nb_add,-BYbus_non(:,3),BYbus_non(:,1)+nb
    BYbus_non(:,2)+nb_add,BYbus_non(:,1)+nb_add,-BYbus_non(:,3),BYbus_non(:,1)+nb
    GYbus_non(:,1),GYbus_non(:,2)+nb_add,-GYbus_non(:,3),GYbus_non(:,1)+nb
    GYbus_non(:,2)+nb_add,GYbus_non(:,1),-GYbus_non(:,3),GYbus_non(:,1)+nb
    GYbus_non(:,1)+nb_add,GYbus_non(:,2),GYbus_non(:,3),GYbus_non(:,1)+nb
    GYbus_non(:,2),GYbus_non(:,1)+nb_add,GYbus_non(:,3),GYbus_non(:,1)+nb
    ];
end

% --------------------- END -------------------- %

idgen = find(gen(:,7)==0);

% -------------- PV节点电压幅值等式约束的海森矩阵 ------------- %
% H2(ei,ei)=2,H2(fi,fi)=2,i=1,2...nb,等式约束处于2*nb+ntb+1:2*nb+ntb+nidgen段
save_9 = [[1:nb]',[1:nb]',2*ones(nb,1),[2*nb+ntb+ntb+1:2*nb+ntb+ntb+nb]'
    [1:nb]'+nb_add,[1:nb]'+nb_add,2*ones(nb,1),[2*nb+ntb+ntb+1:2*nb+ntb+ntb+nb]'];
% --------------------- END -------------------- %
save_9_temp = [];
for i = 1: size(idgen,1)  
if gen(idgen(i,1),1)~=ref
    [temp1,temp2] = find(save_9(:,2) == gen(idgen(i,1),1) );
    save_9_temp =[ save_9_temp;save_9(temp1,:)];
   [temp1,temp2] = find(save_9(:,2) == gen(idgen(i,1),1)++nb_add );
    save_9_temp =[ save_9_temp;save_9(temp1,:)];
else
    [temp1,temp2] = find(save_9(:,2) == gen(idgen(i,1),1) );
    save_9_temp =[ save_9_temp;save_9(temp1,:)];
end
end
save_9_final = save_9_temp;


save_H = [ save_1; save_2; save_7; save_8 ; save_9_final ];


[temp1,temp2] = find(save_H(:,3)==0); 
save_H(temp1,:)=[]; % 海森矩阵中的二阶偏导数生成完毕

% -------------- 处理平衡节点的虚部f（f=0） ------------- %
% 消去save矩阵中平衡节点的虚部f所在行列
[temp1,temp2]=find(save_H(:,1)==ref+nb_add);
save_H(temp1,:)=[];
[temp1,temp2]=find(save_H(:,2)==ref+nb_add);
save_H(temp1,:)=[];
[temp1,temp2]=find(save_H(:,1)>ref+nb_add);
save_H(temp1,1)=save_H(temp1,1)-1;% 行数前进一位
[temp1,temp2]=find(save_H(:,2)>ref+nb_add);
save_H(temp1,2)=save_H(temp1,2)-1; % 列数前进一位
[temp1,temp2]=find(save_H(:,3)==0);
save_H(temp1,:)=[];



% --------------------- END -------------------- %
% save_H保存了海森矩阵的所有非零元素
% ------------------ END ----------------- %