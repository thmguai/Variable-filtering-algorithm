%%无信息变量消除法
x1=xlsread('C:\Users\c430\Desktop\蛋白质含量测试\波长变量筛选\UVE\蛋白质含量预测.xlsx');
x=x1(:,2:end);%光谱矩阵
y=x1(:,1);%浓度矩阵
f=7;%ncomp最佳主因子数
[n,m]=size(x);
R=normrnd(0.4,0.242,n,m);%生成一个随机变量矩阵
XR=[x,R];%x与R进行组合
B=[];
for i=1:n
    xr=XR;
    xr(i,:)=[];
    Y=y;
    Y(i)=[];
    [xl,yl,xs,ys,beta,pctvar,mse]=plsregress(xr,Y,f);%对xr和Y进行pls回归
    B=[B,beta];%回归系数矩阵B
end
me=mean(B,2);%求均值
s=std(B')';%求标准偏差
h=me./s;
h(1)=[];
hmax=max(abs(h(1+m:2*m)));%在区间[1+m,2m]上得到最大值，并设置为阈值
index=find(abs(h) > hmax); %找到h大于阈值的波段的位置并标上索引
plot(1:m,h(1:m),'k',m+1:(2*m),h(m+1:end),'r')
hold on
plot((1:2*m),hmax,'b')
hold on
plot((1:2*m),-1*hmax,'b')
hold on
plot(index,h(index),'g+')
xlabel('real variables - index - random variables')
ylabel('系数均值/标偏')
X=x(:,index);
index = index';

[m,n] = size(index);
for i = 1:n
    y = [y,x(:,index(i))]; 
end
z = fitness_2(y);
disp(z);



