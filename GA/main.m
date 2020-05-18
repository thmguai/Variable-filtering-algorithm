%clear all;
%clc;
tic;
NP=50;%群体个数
NG=200;%迭代次数
Pc=0.85;
Pm=0.1;
L=35; %码长
x=zeros(NP,L);
fx=zeros(1,NP);
X=xlsread('C:\Users\c430\Desktop\蛋白质含量测试\波长变量筛选\GA\蛋白质校正集.xlsx');%导入数据
Y = X(:,1);
X2=xlsread('C:\Users\c430\Desktop\蛋白质含量测试\波长变量筛选\GA\蛋白质预测集.xlsx');%导入数据
Y1 = X2(:,1);
K1 = 50;
T1 = 0;
T2 = 0;
T3 = 0;
T4 = 0;
for k1 = 1:K1
for i=1:NP
    x(i,:)=Initial(L);
    fx(i)=fitness(x(i,:),X,L);
end
BEST = 0;
for k=1:NG
    if mod(k,10) == 0
        disp(k)
    end
    sumfx=sum(fx);
    Px=fx/sumfx;
%     NP1 = randsample(NP,NP,true,Px);
%     for i=1:NP
%         x(i,:)=x(NP1(i),:);
%         fx(i)=fitness(x(i,:));
%     end
    for i=1:NP
        SelFather = randsample(NP,1,true,Px);
        Selmother = randsample(NP,1,true,Px);
        %Selmother=floor(rand()*(NP-1))+1;   %随机选择母亲
        posCut=floor(rand()*(L-2))+1;   %随机确定交叉点
        r1=rand();
        if r1<=Pc
            nx(i,1:posCut)=x(SelFather,1:posCut);
            nx(i,(posCut+1):L)=x(Selmother,(posCut+1):L);
            r2=rand();
            if r2<=Pm    %变异
                posMut=round(rand()*(L-1)+1);
                nx(i,posMut)=~nx(i,posMut);
            end
        else
            nx(i,:)=x(SelFather,:);
        end
    end
    x=nx;
    for i=1:NP
        fx(i)=fitness(x(i,:),X,L);
    end
    [best(k),bestGeneIndex(k)]= max(fx);
    if BEST < best(k)
        BEST = best(k);
        BEST_GENE = x(bestGeneIndex(k),:);
    end
    bestValue(k) = BEST;
    %bestGene(k) = BEST_GENE;
end
figure(3);
plot(1:NG,bestValue);
j=1;
nir_num1(k1) = 0;
for i=1:L
    if BEST_GENE(i)==1;
        X_end1(:,j:j+19) = X(:,20*i-18:20*i+1);
        j = j + 20;
        nir_num1(k1) = nir_num1(k1) + 1;
    end
end
[Rc,RMSEC,beta,yc]= fitaaa(BEST_GENE,X,L);
% figure(1);
% plot(Y,yc,'ro');
% xlabel('Actual Value');
% ylabel('Predictive Value');
% text(7.7,9.8,['Rc = ',num2str(Rc)]);
% text(7.7,9.65,['RMSEC = ',num2str(RMSEC)]);
% hold on;
% plot([7.5,10],[7.5,10],'linewidth',1.5);
[Rp,RMSEP,yp]= fitbbb(BEST_GENE,X2,L,beta);
Rp1(k1) = RMSEP;
T1 = T1 + Rc;
T2 = T2 + RMSEC;
T3 = T3 + Rp;
T4 = T4 + RMSEP;
end
disp(T1/K1);
disp(T2/K1);
disp(T3/K1);
disp(T4/K1);
toc;
% figure(2);
% plot(Y1,yp,'ro');
% xlabel('Actual Value');
% ylabel('Predictive Value');
% text(7.7,9.8,['Rp = ',num2str(Rp)]);
% text(7.7,9.65,['RMSEP = ',num2str(RMSEP)]);
% hold on;
% plot([7.5,10],[7.5,10],'linewidth',1.5);