

function [bus,trsfm,shtcr,tbnode,tbgt,tbbt,tbkmax,tbkmin] = initialIPM(bus,gen,Lnbr,trsfm,shtc,shtr)


% ---------------- 系统节点类型、各类节点个数信息 -------------- %
nb = size(bus,1);
ng = size(gen,1);
nt0 = size(trsfm,1);
if ~isempty(trsfm)
    ntb = size(find(trsfm(:,10)==1),1);
else
    ntb =0;tbnode=[];tbgt=[];tbbt=[];tbkmax=[];tbkmin=[];
end
nl = size(Lnbr,1);
nlt = nl + nt0;
if ~isempty(shtr)
    nr = size(shtr,1);
else
    shtcr=shtc;
end
% ---------------- END -------------- %


% ----------- 构造变压器支路的虚拟节点 ---------- %
if ~isempty(trsfm)
    temp = find(trsfm(:,10) == 1 );
    tbnode = [trsfm(temp,2),[nb+1:nb+ntb]',trsfm(temp,3)];
    trsfm(temp,2) = [nb+1:nb+ntb]'; % 首节点变化转为虚拟节点编号
    tbgtbt = 1./(trsfm(temp,4)+j*trsfm(temp,5));  % 变压器支路阻导纳值处理
    tbgt = real(tbgtbt);
    tbbt = imag(tbgtbt);
    tbkmax = trsfm(temp,7);
    tbkmin = trsfm(temp,8);  % 对应变压器变比的不等式约束
    
    % 引入虚拟节点
    bus(nb+1:nb+ntb,1) = tbnode(:,2);
    bus(nb+1:nb+ntb,8) = bus(tbnode(:,3),8);
end

bus(find(bus(:,9)==0),9)= 1.1;     % 虚拟节点的电压下限
bus(find(bus(:,10)==0),10)= 1.1 ;  % 虚拟节点的电压上限
% --------------------- END -------------------- %

% -------------- 并联电容和并联电抗数据初始化 ------------- %
if ~isempty(shtc)
    % shtc、shtr组合成shtcr矩阵，并定义电容为正方向
    idshtc = find(shtc(:,6)==1);
    shtcr = shtc(idshtc,:);
    if ~isempty(shtr)
    for t = 1: nr
        temp = find(shtcr(:,1)==shtr(t,1));
        if ~isempty(temp) % 并联电容和并联电抗接在同一个位置
            shtcr(temp,4) = - shtr(t,3);              % 跟新shtcr的下限
            shtcr(temp,2) = shtc(temp,2) - shtr(t,2); % 电容器和电抗器的代数和
        else            % 并联电容和并联电抗没有接在同一个位置
            shtr(t,2) = - shtr(t,2);%定义为电容方向，电抗投切量相当于电容投切量的负数；
            shtr(t,4) = - shtr(t,3);
            shtr(t,3) = 0;
            shtcr = [ shtcr; shtr(t,:)];%定义为电容方向的意义正在于此
        end
    end
    end
end


% --------------------- END -------------------- %