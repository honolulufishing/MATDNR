function [bus,gen,trsfm,shtcr,loss,elapsed_time,Record_KKT,success] = VARIPM(bus,gen,Lnbr,trsfm,shtc,shtr,vctr)


% ---------------- ������������ݳ�ʼ�� -------------- %
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

% -------------- ��ѹ��������ݳ�ʼ�� ------------- %
j = sqrt(-1);
if ~isempty(trsfm)
    temp = find(trsfm(:,10) == 1 );
    trsfm(temp,6) = 1;  % ����ѹ���ı�ȳ�ʼ��ֵ����Ϊ1���Ժ��Ż��лᷢ���仯   
    ntb = size(find(trsfm(:,10)==1),1);
else
     ntb =0;
end
% --------------------- END -------------------- %

% ��������ڵ�
[bus,trsfm,shtcr,tbnode,tbgt,tbbt,tbkmax,tbkmin] = initialIPM(bus,gen,Lnbr,trsfm,shtc,shtr);

% ���ɺ�ɭ����
[save_H,GYbus,BYbus ,GYbus_non,BYbus_non,Ybus_non] = hessianIPM(bus,gen,Lnbr,trsfm,shtc,shtr,tbnode,tbgt,tbbt);

% -------------- �Ż�����������Ϣ��ʼ�� ------------- %
nb = size(bus_s,1);

refgen = find(bus(:, 2) == 3);
ref = bus(refgen,1);                   % reference bus index
npvg = size(gen,1);                    % PV��slack�ڵ�
ncr = size(shtcr,1);                   %  shtcr����ı�������
n = nb+ntb+nb-1+ntb+ntb+npvg+ncr;      %  ��������
nb_add = nb + ntb;                     %  ��ģ�͵Ľڵ����
% --------------------- END -------------------- %

% -------------- ����ʽԼ����������lb��ub ------------- %
% ����ʽԼ����������
if ~isempty(trsfm)
    lb = zeros(nb+ntb+npvg+ncr,1);   % ����ʽ����Լ��
    ub = zeros(nb+ntb+npvg+ncr,1);   % ����ʽ����Լ��
    lb(1:nb,1) = bus(1:nb,9).^2;     % ��ѹ��ֵ����
    ub(1:nb,1) = bus(1:nb,10).^2;    % ��ѹ��ֵ����
    lb(nb+1:nb+ntb,1) = tbkmin;      % �������ֵԼ��
    ub(nb+1:nb+ntb,1) = tbkmax;      % �������ֵԼ��
else
    lb = zeros(nb+npvg,1);   % ����ʽ����Լ��
    ub = zeros(nb+npvg,1);   % ����ʽ����Լ��
    lb(1:nb,1) = bus(1:nb,9).^2;     % ��ѹ��ֵ����
    ub(1:nb,1) = bus(1:nb,10).^2;    % ��ѹ��ֵ����
end

temp_ref = find(gen(:,1)==ref);
temp_pv = setdiff(1:size(gen,1),temp_ref)';
temp =[ temp_ref; temp_pv ];   % ��ƽ��ڵ����gen��1λ��������ڵ���û��PQ�ڵ�

lb(nb+ntb+1:nb+ntb+npvg,1) = gen(temp,5); % ������޹���������
ub(nb+ntb+1:nb+ntb+npvg,1) = gen(temp,4); % ������޹���������

if ~isempty(shtcr)    
lb(nb+ntb+npvg+1:nb+ntb+npvg+ncr,1) = shtcr(:,4); % �޹���������
ub(nb+ntb+npvg+1:nb+ntb+npvg+ncr,1) = shtcr(:,3); % �޹���������
end


% --------------------- END -------------------- %

% -------------- �ڵ㷨����ʼ������x0 ------------- %
% x��ʼ����x=[e f k QG Qcri]
x = zeros(n,1);               % ��������

idgen = find(gen(:,7)==0);
nidgen = size(idgen,1);

y = zeros(2*nb-1+ntb+ntb+nidgen,1);  % ��ʽԼ����������������
h = zeros(nb+ntb+npvg+ncr,1); % ����ʽԼ����������
x(1:nb,1) = (bus(1:nb,9)+bus(1:nb,10))/2;  % ��ͨ�ڵ��ѹʵ���ĳ�ʼֵ
if ~isempty(trsfm)
    x(nb+1:nb+ntb,1) = x(tbnode(:,1),1);            % ����Ľڵ��ѹʵ���ĳ�ʼֵ
end
x(nb+ntb+1:2*(nb+ntb)-1,1) = 0;                 % ���нڵ㣨����ƽ��ڵ���⣩�ĵ�ѹ��ǵĳ�ֵ��Ϊ0�����Զ�Ӧ�鲿ҲΪ0
x(2*(nb+ntb):n,1) = (lb(nb+1:nb+ntb+npvg+ncr,1)+ub(nb+1:nb+ntb+npvg+ncr,1))/2; % �����ע���޹��������޹��ĵĳ�ʼֵ

% ------------------ END ----------------- %

% -------------- �ڵ㷨����ż����y,l,u,z,w�ĳ�ʼ�� ------------- %
y(1:nb-1,1) = -10;    % y �������ն�ż��������Ӧ���ʵ�ʽԼ����
y(2*nb-1+ntb+ntb+1:2*nb-1+ntb+ntb+nidgen) = 10; % y �������ն�ż��������ӦPV�ڵ��ʽԼ����
e = x(1:nb+ntb,1);                                               % ���нڵ��ʵ��ֵ
f = [x(nb+ntb+1:nb+ntb+ref-1,1);0;x(nb+ntb+ref:2*nb+2*ntb-1,1)]; % ���нڵ���鲿ֵ
e_old = x(1:nb,1);                                                 % ������ͨ�ڵ��ʵ��ֵ
f_old = [x(nb+ntb+1:nb+ntb+ref-1,1);0;x(nb+ntb+ref:2*nb+ntb-1,1)]; % ������ͨ�ڵ���鲿ֵ
% ����ʽ��h(x)��ֵ��ȡ
h = [ e_old.^2+f_old.^2 ; x(2*(nb+ntb):n,1)];
% ���㲻��ʽԼ������ֵ h
l = h-lb;
u = ub-h;
z = 10*ones(nb+ntb+npvg+ncr,1);
w = z;                % w,z �������ն�ż��������Ӧ����ʽԼ����
% ------------------ END ----------------- %

% -------------- �ڵ㷨�����ſ˱Ⱦ��� ------------- %
% [Jf,Jg,Jh] = jacobiIPM(bus,gen,trsfm,shtcr,x,e,f,e_old,f_old,GYbus,BYbus,GYbus_non,BYbus_non,Ybus_non,tbnode,tbgt,tbbt);
% ------------------ END ----------------- %

% -------------- �����ڵ㷨���Ĳ�����Ϣ ------------- %
Gap = l'*z+u'*w;   % ���㻥����϶�ĳ�ʼֵ
MaxIter = 200;    % ����������
tol_1 = 10^(-6);   % ������϶��������
tol_2 = 10^(-6);   % KKT�����������������
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

    % �����ʽԼ������ֵg
    g = geqIPM(bus,gen,trsfm,shtcr,x,GYbus,BYbus,e,f,e_old,f_old,tbnode,tbgt,tbbt);
    
    % ����Ŀ�꺯��goalf
    goalf = g(ref,:);    % Ŀ�꺯����ֵ
    g(ref,:) = [];       % ɾ��ƽ��ڵ��й�ע������Ӧ����
    
    
    % -------------- �ڵ㷨�����ſ˱Ⱦ��� ------------- %
    [Jf,Jg,Jh] = jacobiIPM(bus,gen,trsfm,shtcr,x,e,f,e_old,f_old,GYbus,BYbus,GYbus_non,BYbus_non,Ybus_non,tbnode,tbgt,tbbt);
    % ------------------ END ----------------- %
    
    
    % --------- ����KKT���� -------- %
    Lx = Jf - Jg'*y - Jh'*(z-w);  % ����Lx
    Ly = g;                       % ����Ly
    Lz = h - l - lb;              % ����Lz
    Lw = h + u - ub;              % ����Lw
    % ------------- END ------------- %
    
    % ��¼��ii�ε�����Ŀ�꺯��ֵ
    Record_goalf(ii,1) = goalf;
    % ��¼��ii�ε�����KKT���̲�ƽ����
    Record_KKT(ii,:) = [max(abs(Lx)) max(abs(Ly))  max(abs(Lz))  max(abs(Lw))];
    Record_Gap(ii,:) = Gap;
    
    % Gap��KT�����������о�
    if Gap < tol_1                                      % Gap�������о�
        if max([abs(Lx);abs(Ly);abs(Lz);abs(Lw)])<tol_2 % KKT�������о�
            
            success = 1;
            
            % �й�����
            j = sqrt(-1);
            gen(find(gen(:,1)==ref),2) = 0;
            loss = goalf + sum(gen(:,2),1) + sum(bus(:,3),1);
            % ��ѹ��ֵ�����
            V = abs(e+j*f);
            bus(:,7) =  V;
            bus(:,8) = angle(e+j.*f);
            % ��������˵�ѹ���޹�����
            gen(:,6) = bus(gen(:,1),7);
            gen(temp,3) = x(2*(nb+ntb)+ntb:2*(nb+ntb)+ntb-1+npvg,1);
            gen(temp_ref,2) = goalf;
            if ~isempty(shtc)
            % ���ݵ翹��ע���޹�
            idshtc = find(shtc(:,6)==1);
            shtcr(idshtc ,2) = x(2*(nb+ntb)+ntb+npvg:n,1);
            end
            if ~isempty(trsfm)
            % ��ѹ�����
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
 
    % --------------- �����䷽�� ------------------ %
    % ��Ŀ�꺯������ʽԼ��������ʽԼ���ĺ�ɭ�����������ϵ���ϲ���һ��
    yy = [-y(1:ref-1,1);1;-y(ref:2*nb-1+ntb+ntb+nidgen,1);-(z-w)];
    H = sparse(save_H(:,1),save_H(:,2),save_H(:,3).*yy(save_H(:,4)),n,n);
    ln1 = reciprocal(l);
    un1 = reciprocal(u);
    vmid = z.*ln1 + w.*un1;
    H = H+Jh'*diag(vmid,0)*Jh;     % ����H
  
    vmid = (z.*l+z.*Lz).*ln1-(w.*u-w.*Lw).*un1;
    kesai_aff = Lx + Jh'*vmid;
    % ����kesai_aff
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
    
    % ȷ�����䲽����stepp_aff,stepd_aff,������϶��Gap_aff
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
    % ���Ĳ���sigema,Ԥ���ϰ�����mu
    %     sigema=(Gap_aff/Gap)^3;
    %     mu=sigema*Gap/2/(nb+2*ntb+npvg+ncr);
    mu=min((Gap_aff/Gap)^2,0.2)*Gap_aff/2/(nb+ntb+npvg+ncr);
    
    % -----------������У������------------- %
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
   
    % ----------ֱ�����У������---------- %
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
    
    % ȷ������������stepp,stepd
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
    
    % ����ԭ�����Ͷ�ż����x,y,l,u,z,w
    x = x + stepp*dx;
    y = y + stepd*dy;
    l = l + stepp*dl;
    u = u + stepp*du;
    z = z + stepd*dz;
    w = w + stepd*dw;
    
    % �õ�����һ�εı���ֵx,y,l,u,z,w
    e = x(1:nb+ntb,1);                            % ���нڵ��ʵ��ֵ
    f = [x(nb+ntb+1:nb+ntb+ref-1,1);0;x(nb+ntb+ref:2*nb+2*ntb-1,1)];% ���нڵ���鲿ֵ
    e_old = x(1:nb,1);                            % ������ͨ�ڵ��ʵ��ֵ
    f_old = [x(nb+ntb+1:nb+ntb+ref-1,1);0;x(nb+ntb+ref:2*nb+ntb-1,1)];% ������ͨ�ڵ���鲿ֵ
    h = [e_old.^2+f_old.^2;x(2*(nb+ntb):n,1)];    % ���㲻��ʽԼ������ֵ

    Gap = l'*z+u'*w;     % ����Gap

    % -------------- �ڵ㷨�����ſ˱Ⱦ��� ------------- %
    % [Jf,Jg,Jh] = jacobiIPM(bus,gen,trsfm,shtcr,x,e,f,e_old,f_old,GYbus,BYbus,GYbus_non,BYbus_non,Ybus_non,tbnode,tbgt,tbbt);
    % ------------------ END ----------------- %

    if ii == MaxIter%�����������������о�
        disp('Warning:�ڵ㷨�ﵽ�����ѭ���������޷�������');
        break;
    end
    
    ii = ii + 1;    % ��������ii����1
end

elapsed_time =  cputime- tStart;



% ����ڵ��޳�
bus = bus(1:nb,:);
% % �������ڳ�������
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

