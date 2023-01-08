

function [bus,trsfm,shtcr,tbnode,tbgt,tbbt,tbkmax,tbkmin] = initialIPM(bus,gen,Lnbr,trsfm,shtc,shtr)


% ---------------- ϵͳ�ڵ����͡�����ڵ������Ϣ -------------- %
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


% ----------- �����ѹ��֧·������ڵ� ---------- %
if ~isempty(trsfm)
    temp = find(trsfm(:,10) == 1 );
    tbnode = [trsfm(temp,2),[nb+1:nb+ntb]',trsfm(temp,3)];
    trsfm(temp,2) = [nb+1:nb+ntb]'; % �׽ڵ�仯תΪ����ڵ���
    tbgtbt = 1./(trsfm(temp,4)+j*trsfm(temp,5));  % ��ѹ��֧·�赼��ֵ����
    tbgt = real(tbgtbt);
    tbbt = imag(tbgtbt);
    tbkmax = trsfm(temp,7);
    tbkmin = trsfm(temp,8);  % ��Ӧ��ѹ����ȵĲ���ʽԼ��
    
    % ��������ڵ�
    bus(nb+1:nb+ntb,1) = tbnode(:,2);
    bus(nb+1:nb+ntb,8) = bus(tbnode(:,3),8);
end

bus(find(bus(:,9)==0),9)= 1.1;     % ����ڵ�ĵ�ѹ����
bus(find(bus(:,10)==0),10)= 1.1 ;  % ����ڵ�ĵ�ѹ����
% --------------------- END -------------------- %

% -------------- �������ݺͲ����翹���ݳ�ʼ�� ------------- %
if ~isempty(shtc)
    % shtc��shtr��ϳ�shtcr���󣬲��������Ϊ������
    idshtc = find(shtc(:,6)==1);
    shtcr = shtc(idshtc,:);
    if ~isempty(shtr)
    for t = 1: nr
        temp = find(shtcr(:,1)==shtr(t,1));
        if ~isempty(temp) % �������ݺͲ����翹����ͬһ��λ��
            shtcr(temp,4) = - shtr(t,3);              % ����shtcr������
            shtcr(temp,2) = shtc(temp,2) - shtr(t,2); % �������͵翹���Ĵ�����
        else            % �������ݺͲ����翹û�н���ͬһ��λ��
            shtr(t,2) = - shtr(t,2);%����Ϊ���ݷ��򣬵翹Ͷ�����൱�ڵ���Ͷ�����ĸ�����
            shtr(t,4) = - shtr(t,3);
            shtr(t,3) = 0;
            shtcr = [ shtcr; shtr(t,:)];%����Ϊ���ݷ�������������ڴ�
        end
    end
    end
end


% --------------------- END -------------------- %