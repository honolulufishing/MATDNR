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
% -------------- �Ż�����������Ϣ��ʼ�� ------------- %
%ncr = size(shtcr,1);                 %  shtcr����ı�������
%n = nb+ntb+nb-1+ntb+ntb+npvg+ncr;    %  �������� 
nb_add = nb + ntb;                     %  ��ģ�͵Ľڵ���� 
% --------------------- END -------------------- %

% -------------- ������ģ���µĽڵ㵼�ɾ���Ybus------------- %
[Ybus,Yf,Yt] = makeYbus(bus,Lnbr,trsfm,shtc,shtr,pq,pv);

GYbus = real(Ybus(1:nb,:));   % Ybusʵ��: nb��(nb+ntb)ά����
BYbus = imag(Ybus(1:nb,:));   % Ybus�鲿: nb��(nb+ntb)ά����

% ��Ѱ��Ybus�еķ���Ԫ��  begin
GYbus_non = [];
BYbus_non = [];
% GYbus_non=[����Ԫ�ص��кţ�����Ԫ�ص��кţ�����Ԫ��];
% BYbus_non=[����Ԫ�ص��кţ�����Ԫ�ص��кţ�����Ԫ��];
for t = 1: nb
    [row,column,non] = find(Ybus(t,:));
    GYbus_non = [ GYbus_non;[t*row',column',non.']];
end
Ybus_non = GYbus_non;
BYbus_non =[ GYbus_non(:,1:2),imag(GYbus_non(:,3)) ]; % �ҳ��鲿����Ԫ��ά������nnz(BYbus)
GYbus_non(:,3) = real(GYbus_non(:,3));                % �ҳ�ʵ������Ԫ��ά������nnz(GYbus)
% --------------------- END -------------------- %

% ------------------ ���ɺ�ɭ���� ----------------- %
% -------------- �����ȵĵ�ʽԼ���ĺ�ɭ���� ------------- %
if ~isempty(trsfm)
% save=[����Ԫ�ص��кţ�����Ԫ�ص��кţ�����Ԫ�أ���������һ����ɭ����]
% H1��ei,fm��=1,H1(em,fi)=-1 ����ʽԼ������2*nb+1:2*nb+ntb�� 
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

% -------------- �ڵ��ѹ��ֵ�����޲���ʽԼ���ĺ�ɭ���� ------------- %
% H2(ei,ei)=2,H2(fi,fi)=2,i=1,2...nb
save_2=[[1:nb]',[1:nb]',2*ones(nb,1),[2*nb+ntb+ntb+1:3*nb+ntb+ntb]'
    [1:nb]'+nb_add,[1:nb]'+nb_add,2*ones(nb,1),[2*nb+ntb+ntb+1:3*nb+ntb+ntb]'];
% --------------------- END -------------------- %

% -------------- �ڵ�ע���й����޹���ʽԼ���ĺ�ɭ���� ------------- %
% save_7Ϊ�й�����,save_8Ϊ�޹����ʷ��̵Ķ���ƫ����
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

% -------------- PV�ڵ��ѹ��ֵ��ʽԼ���ĺ�ɭ���� ------------- %
% H2(ei,ei)=2,H2(fi,fi)=2,i=1,2...nb,��ʽԼ������2*nb+ntb+1:2*nb+ntb+nidgen��
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
save_H(temp1,:)=[]; % ��ɭ�����еĶ���ƫ�����������

% -------------- ����ƽ��ڵ���鲿f��f=0�� ------------- %
% ��ȥsave������ƽ��ڵ���鲿f��������
[temp1,temp2]=find(save_H(:,1)==ref+nb_add);
save_H(temp1,:)=[];
[temp1,temp2]=find(save_H(:,2)==ref+nb_add);
save_H(temp1,:)=[];
[temp1,temp2]=find(save_H(:,1)>ref+nb_add);
save_H(temp1,1)=save_H(temp1,1)-1;% ����ǰ��һλ
[temp1,temp2]=find(save_H(:,2)>ref+nb_add);
save_H(temp1,2)=save_H(temp1,2)-1; % ����ǰ��һλ
[temp1,temp2]=find(save_H(:,3)==0);
save_H(temp1,:)=[];



% --------------------- END -------------------- %
% save_H�����˺�ɭ��������з���Ԫ��
% ------------------ END ----------------- %