function IMGA
clear all;
clc;
%x1属于[-2.048,2.048],x2属于[-2.048,2.048],fitness:J=100*(x1^2-x2)^2+(1-x1)^2;
popsize=50;% the value of population 种群取复数
%CodeL=2;%基因个数
len=35;
G=200;  % the max generation
mnum=8; % mnum---记忆细胞数目
Pc=0.85;
Pm=0.05;
Tacl=0.8;     %计算抗体浓度
Best_result=3905.926;
X=xlsread('C:\Users\c430\Desktop\蛋白质含量测试\波长变量筛选\IGA\蛋白质校正集.xlsx');%导入数据
Y = X(:,1);
X2=xlsread('C:\Users\c430\Desktop\蛋白质含量测试\波长变量筛选\IGA\蛋白质预测集.xlsx');%导入数据
Y1 = X2(:,1);
K1 = 1;
T1 = 0;
T2 = 0;
T3 = 0;
T4 = 0;
for k1 = 1:K1
%---------------------初始化种群-------------------------%
for i=1:popsize
    x(i,:)=Initial(len);
    fx(i)=fitness(x(i,:),X,len);
end

%--------------------------------------------------------%
[Ag,RealValue]=affinity1(x,popsize,X,len);%抗体与抗原的亲和力
%-------------------产生初始记忆细胞----------------------%
for i=1:mnum 
    mcell(i,:)=Initial(len);
    mfx(i)=fitness(mcell(i,:),X,len);
end
%--------------------------------------------------------%
[mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len);%更新初始记忆细胞
[Order_Ag,Index_Ag]=sort(Ag);%初始记忆细胞更新初始种群
x(Index_Ag(1:mnum),:)=mx;
[Ag,RealValue]=affinity1(x,popsize,X,len);%抗体与抗原的亲和力
[Ab]=affinity2(x,popsize,len);%抗体与抗体之间的亲和力
C=(sum(Ab>Tacl,2)/popsize)';%计算抗体浓度
[Best_Ag,Index]=max(Ag);%求初始代精英抗体
Best_Value(1)=RealValue(Index);
Temp_Value=Best_Value(1);
Best_gene=x(Index,:);
[Worst_Ag,Index1]=min(Ag);%求初始代最差抗体
Worst_Value(1)=RealValue(Index1);
Worst_gene=x(Index1,:);
avg_Ag=0;
for k=1:popsize
    avg_Ag=avg_Ag+RealValue(k);
end
avg_Ag=avg_Ag/popsize;%求初始代抗体与抗原亲和度平均值
fprintf(1,'\n1--->Best : %f Worst : %f Avg : %f \n',Best_Value(1),Worst_Value(1),avg_Ag)
[mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len);%更新初始记忆细胞
%-------------------------------开始进化----------------------------------%
time(1)=1;
for gen=2:G
  time(gen)=gen;  
%---------------------1:Select Cross Mutation Operation----------------%
[E]=Select_Antibody(C,Ag,popsize);
[newench2select]=Select(x,E,popsize);%复制个体
[newench2cross]=Cross(newench2select,Pc);%交叉
[ench2]=Mutate(newench2cross,Pm);%变异
%[chrom]=decode(ench2,len,MinX,MaxX);%解码
%-----------------------------2:免疫操作-------------------------------%
[Ag,RealValue]=affinity1(ench2,popsize,X,len);%抗体与抗原的亲和力
[Ab]=affinity2(ench2,popsize,len);%抗体与抗体之间的亲和力
C=(sum(Ab>Tacl,2)/popsize)';%计算抗体浓度
[Best_cur_Ag,Index]=max(Ag);
Best_cur_gene=ench2(Index,:);
[Worst_cur_Ag,index]=min(Ag);
Worst_Value=RealValue(index);
Worst_cur_gene=ench2(index);
avg_Value=0;
for k=1:popsize
    avg_Value=avg_Value+RealValue(k);
end
avg_Value=avg_Value/popsize;
%精英抗体保留策略
if Best_cur_Ag < Best_Ag
    Best_Value(gen)=Temp_Value;
    ench2(index,:)=Best_gene;
elseif Best_cur_Ag >= Best_Ag
    Best_Value(gen)=RealValue(Index);
    Best_Ag=Best_cur_Ag;
    Temp_Value=RealValue(Index);
    Best_gene=Best_cur_gene;
end
fprintf(1,'\n%d--->BEST : %f WORST : %f AVG : %f \n',gen,Best_Value(gen),Worst_Value,avg_Value);
fprintf(1,'--->Value_F = %f',Best_Value(gen));
disp(Best_gene)
% if abs(Best_Value(gen)-Best_result)<0.001
%     fprintf(1,'在第%d代收敛到最优解 ',gen);
%     return;
% end
[mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len);%更新初始记忆细胞
[~,Index_Ag]=sort(Ag);%更新本代的基因
ench2(Index_Ag(1:mnum),:)=mx;
x=ench2;
Ag(Index_Ag(1:mnum))=mAg;
end
j=1;
for i=1:len
    if Best_gene(i)==1;
        X_end(:,j:j+19) = X(:,20*i-18:20*i+1);
        j = j + 20;
    end
end
[Rc,RMSEC,beta,yc]= fitaaa(Best_gene,X,len);
% figure(1);
% plot(Y,yc,'ro');
% xlabel('Actual Value');
% ylabel('Predictive Value');
% text(7.7,9.8,['Rc = ',num2str(Rc)]);
% text(7.7,9.65,['RMSEC = ',num2str(RMSEC)]);
% hold on;
% plot([7.5,10],[7.5,10],'linewidth',1.5);
[Rp,RMSEP,yp]= fitbbb(Best_gene,X2,len,beta);
T1 = T1 + Rc;
T2 = T2 + RMSEC;
T3 = T3 + Rp;
T4 = T4 + RMSEP;
% figure(2);
% plot(Y1,yp,'ro');
% xlabel('Actual Value');
% ylabel('Predictive Value');
% text(7.7,9.8,['Rp = ',num2str(Rp)]);
% text(7.7,9.65,['RMSEP = ',num2str(RMSEP)]);
% hold on;
% plot([7.5,10],[7.5,10],'linewidth',1.5);
end
disp(T1/K1);
disp(T2/K1);
disp(T3/K1);
disp(T4/K1);
%------------------------------画图--------------------------------%
%save IMGA_MAX_F1 time Best_Value
figure(3);
plot(time,Best_Value,'b');
xlabel('进化代数');ylabel('目标函数值');% the value of objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%免疫遗传算法子函数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
%---------------------计算抗体与抗原的亲和度----------------%
function [Ag,Real]=affinity1(x,popsize,X,len) 
for i=1:popsize
    x1=x(i,:);
    Real(i) = fitness(x1,X,len);
    Ag(i)=Real(i)/1;%在求取极大值时，两者相等；求极小值时，转为求最大值问题，并归一化
end
%---------------------------------------------------------%
%---------------------计算抗体与抗体的亲和度----------------%
function [Ab]=affinity2(ench2,popsize,len)
for i=1:popsize
    for j=i:popsize
        gene1=ench2(i,:);
        gene2=ench2(j,:);
        s=0;
        for k=1:len
            if gene1(k)==gene2(k)
                Hj(k)=0;
            else
                Hj(k)=log10(2);
            end
            s=s+Hj(k);
        end
        H_2=s/len;
        Ab(i,j)=1/(1+H_2);
        Ab(j,i)=Ab(i,j);
    end
end
%------------------------------------------------------------------------%
%----------------------------更新记忆细胞---------------------------------%
function [mx,mAg]=memorycell(x,mcell,mnum,popsize,X,len)
snum=mnum+popsize;
for i=1:popsize
    x1=x(i,:);
    GxAg(i)=fitness(x1,X,len);
end
for i=1:mnum
    GxmAg(i)=fitness(mcell(i,:),X,len);
end
TAg(1:popsize)=GxAg;
Tch(1:popsize,:)=x;
TAg(popsize+1:popsize+mnum)=GxmAg;
Tch(popsize+1:popsize+mnum,:)=mcell;
[OrderTAg,IndexTAg]=sort(TAg);
for i=1:mnum
    mx(i,:)=Tch(IndexTAg(snum+1-i),:);
    mAg(i)=OrderTAg(snum+1-i);
end

%[mgene,sumlen]=encode(mchrom,len,mnum,MaxX,MinX);
%------------------------------------------------------------------------%
%-------------------------------量化抗体的标准--------------------------------%
function [E]=Select_Antibody(C,Ag,popsize)
lanm=0.7;
miu=1.25;
for i=1:popsize
    E(i)=lanm*Ag(i)+(1-lanm)*exp(-miu*C(i));
end
%------------------------------------------------------------------------%
%-----------------------------轮盘赌策略选择抗体--------------------------%
function [newchselect]=Select(ench2,E,popsize)
sumE=sum(E);
pE=E/sumE;
psE=0;
psE(1)=pE(1);
for i=2:popsize
    psE(i)=psE(i-1)+pE(i);
end
for i=1:popsize
    sita=rand;
    for g=1:popsize
        if sita<=psE(g)
            n=g;
            break;
        end
    end
    newchselect(i,:)=ench2(n,:);
end
%------------------------------------------------------------------------%
%-----------------------------单点交叉法--------------------------------%
function [newench2cross]=Cross(newench2select,Pc)
[Ordersj,Indexsj]=sort(rand(size(newench2select,1),1));%随机交叉
newench2select=newench2select(Indexsj,:);
lchrom=size(newench2select,2);
poscut=ceil(rand(size(newench2select,1)/2,1)*(lchrom-1));%选好随机交叉点
poscut=poscut.*(rand(size(poscut))<Pc); %根据交叉概率进行交叉
for i=1:length(poscut),%length（A）表示A的最大行/列数
     newench2cross([2*i-1 2*i],:)=[newench2select([2*i-1 2*i],1:poscut(i)) newench2select([2*i 2*i-1],poscut(i)+1:lchrom)];%单点交叉
end
%------------------------------------------------------------------------%
%----------------------------------变异----------------------------------%
function [newchmutate]=Mutate(newchcross,Pm)
point=find(rand(size(newchcross))<Pm);%发现变异点
newchmutate=newchcross;
newchmutate(point)=1-newchcross(point);   
%------------------------------------------------------------------------%
function z = fitness(x,X,L)
    j = 1;
    f = 7;
    for i=1:L
        if x(i) == 1
            XX(:,j:j+19) = X(:,20*i-18:20*i+1);
            j = j+20;
        end
    end
    YY = X(:,1);
    [xl,yl,xs,ys,beta,pctvar,mse]=plsregress(XX,YY,f);%对xr和Y进行pls回归
    RMSEC = sqrt(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2)/50);

    R = sqrt(1-(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2))/(sum((YY-mean(YY)).^2)));

    z = R/(1+RMSEC);
    
function result=Initial(length)     %初始化函数
    for i=1:length
        r=rand();
        result(i)=round(r);
    end
    
