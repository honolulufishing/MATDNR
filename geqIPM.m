
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
nb = size(bus,1) - ntb;                % ԭʼϵͳ�ڵ����
npvg = size(gen,1);
nb_add = nb + ntb;                     %  ��ģ�͵Ľڵ����

ncr = size(shtcr,1);                   %  shtcr����ı�������

n = nb+ntb+nb-1+ntb+ntb+npvg+ncr;      %  ��������


% % % ------------- ��������ƽ�ⷽ�̼��㷽��(1) -------------- %
% Pbus = -bus(1:nb,3);
% Qbus = -bus(1:nb,4);
% % �ڵ�ע���й����޹�(�����˸����ڵ�ĸ����й����޹����Ƽ��˷������ע���й�������Pbus,Qbus�ǲ����)
% Pbus(gen(:,1),1) = Pbus(gen(:,1),1) + gen(:,2);  % �ѷ�����е��й�ע�����ӵ�Pbus����
% 
%     %----------�����ʽԼ������ֵg-----------%
%     g = zeros(2*nb+ntb+ntb,1);%�Ե�ʽԼ������ֵ���г�ʼ��
%     % ��������ڵ�ע���й���ʽԼ���ĺ���ֵ
%     gP = (GYbus*e-BYbus*f).*e_old+(BYbus*e+GYbus*f).*f_old - Pbus;
%     % ��������ڵ�ע���޹���ʽԼ���ĺ���ֵ
%     gQ = (GYbus*e-BYbus*f).*f_old-(BYbus*e+GYbus*f).*e_old-Qbus;
% 
%     gQ([ref;pv],1) = gQ([ref;pv],1) - x(2*(nb+ntb)+ntb:2*(nb+ntb)-1+ntb+npvg,1); % �÷����ע���޹�������gQ
%     gQ(shtcr(:,1),1) = gQ(shtcr(:,1),1) - x(2*(nb+ntb)+ntb+npvg:n,1);
% 
%     % �õ��ݵ翹ע���޹�������gQ
%     % ��������ڵ�ע���޹���ʽԼ���ĺ���ֵ
%     % �����Ƕ�ԭ�ɵ���ѹ����ѹ�˽ڵ��ע���С��޹���������
%     eii = e(tbnode(:,2),1);fii=f(tbnode(:,2),1);
%     ej = e(tbnode(:,3),1);fj=f(tbnode(:,3),1);
%     % �ֱ������ڵ��ѹ��ʵ�鲿��ԭ�ɵ���ѹ����ѹ�˵Ľڵ��ѹʵ�鲿��ֵ
%     tbgP = (eii.^2+fii.^2).*tbgt-(eii.*ej+fii.*fj).*tbgt+(eii.*fj-ej.*fii).*tbbt;%����ɵ���ѹ��ԭ��ѹ�˽ڵ�ע���й�Ӧ�����ӵĲ���
%     tbgQ = -(eii.^2+fii.^2).*tbbt+(eii.*ej+fii.*fj).*tbbt+(eii.*fj-ej.*fii).*tbgt;%����ɵ���ѹ��ԭ��ѹ�˽ڵ�ע���޹�Ӧ�����ӵĲ���
%   for m = 1:ntb
%         t = tbnode(m,1);%ȡ���ɵ���ѹ��ԭ��ѹ�˽ڵ�
%         gP(t,1) = gP(t,1) + tbgP(m,1);%��ԭ�ɵ���ѹ����ѹ�˽ڵ��ע���й���������
%         gQ(t,1) = gQ(t,1) + tbgQ(m,1);
%   end
%  % ------------- END -------------- %

% % ------------- ��������ƽ�ⷽ�̼��㷽���� -------------- %
Pbus = bus(1:nb,3);
Qbus = bus(1:nb,4);
% �ڵ�ע���й����޹�(�����˸����ڵ�ĸ����й����޹����Ƽ��˷������ע���й�������Pbus,Qbus�ǲ����)
Pbus(gen(:,1),1) = Pbus(gen(:,1),1) + gen(:,2);  % �ѷ�����е��й�ע�����ӵ�Pbus����

% g = zeros(2*nb+ntb+ntb,1);   % �Ե�ʽԼ������ֵ���г�ʼ��
j = sqrt(-1);

U =  sqrt(e.^2+f.^2).*(cos( atan(f./e))+j*sin( atan(f./e)));
Sbus = Pbus + j*Qbus ;

Sbus([ref;pv],:) = Sbus([ref;pv],:) + j*x(2*(nb+ntb)+ntb:2*(nb+ntb)-1+ntb+npvg,1); % �÷����ע���޹�������gQ
if ~isempty(shtcr)
Sbus(shtcr(:,1),:) = Sbus(shtcr(:,1),:)  + j*x(2*(nb+ntb)+ntb+npvg:n,1);
end

Y = GYbus+j.*BYbus;

mis = U(1:nb,1) .* conj(Y * U) - Sbus;

delta_P = real(mis(1:nb,:));
delta_Q = imag(mis(1:nb,:));

gP = delta_P;
gQ = delta_Q;

    % �õ��ݵ翹ע���޹�������gQ
    % ��������ڵ�ע���޹���ʽԼ���ĺ���ֵ
    
    if ~isempty(trsfm)
    % �����Ƕ�ԭ�ɵ���ѹ����ѹ�˽ڵ��ע���С��޹���������
    eii = e(tbnode(:,2),1);fii=f(tbnode(:,2),1);
    ej = e(tbnode(:,3),1);fj=f(tbnode(:,3),1);
    % �ֱ������ڵ��ѹ��ʵ�鲿��ԭ�ɵ���ѹ����ѹ�˵Ľڵ��ѹʵ�鲿��ֵ
    tbgP = (eii.^2+fii.^2).*tbgt-(eii.*ej+fii.*fj).*tbgt+(eii.*fj-ej.*fii).*tbbt;%����ɵ���ѹ��ԭ��ѹ�˽ڵ�ע���й�Ӧ�����ӵĲ���
    tbgQ = -(eii.^2+fii.^2).*tbbt+(eii.*ej+fii.*fj).*tbbt+(eii.*fj-ej.*fii).*tbgt;%����ɵ���ѹ��ԭ��ѹ�˽ڵ�ע���޹�Ӧ�����ӵĲ���
  for m = 1:ntb
        t = tbnode(m,1);%ȡ���ɵ���ѹ��ԭ��ѹ�˽ڵ�
        gP(t,1) = gP(t,1) + tbgP(m,1);%��ԭ�ɵ���ѹ����ѹ�˽ڵ��ע���й���������
        gQ(t,1) = gQ(t,1) + tbgQ(m,1);
  end
 
% ------------- END -------------- %


% ��ԭ�ɵ���ѹ����ѹ�˽ڵ��ע���С��޹���������
eii = e(tbnode(:,2),1);fii=f(tbnode(:,2),1);
ei = e(tbnode(:,1),1);fi=f(tbnode(:,1),1);%�ֱ��ԭ�ɵ���ѹ����ѹ�˵Ľڵ��ѹʵ�鲿��ֵ
ki = x((nb+ntb)*2:(nb+ntb)*2+ntb-1,1);
g_theta = ei.*fii-eii.*fi;%����ɵ���ѹ��ԭ��ѹ�˽ڵ������ڵ�֮��ĵ�ѹ�����ȵĵ�ʽԼ���ĺ���ֵ
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