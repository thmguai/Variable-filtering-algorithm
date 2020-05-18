clear all;
X1=xlsread('C:\Users\c430\Desktop\蛋白质含量测试\波长变量筛选\PCA\蛋白质含量预测.xlsx');
X = X1(:,2:end);
Y = X1(:,1);
z=zscore(X);               %数据标准化
M=cov(z);                  %协方差
[V,D]=eig(M);             %求出协方差矩阵的特征向量、特征根
d=diag(D);                %取出特征根矩阵列向量（提取出每一主成分的贡献率）
eig1=sort(d,'descend');     %将贡献率按从大到小元素排列
v=fliplr(V);              %依照D重新排列特征向量
S=0;
i=20;
gxl = (sum(eig1(1:i)))/(sum(eig1));
% while S/sum(eig1)<0.99997
%     i=i+1;
%     S=S+eig1(i);
% end                         %求出累积贡献率大于85%的主成分
NEW=z*v(:,1:i);              %输出产生的新坐标下的数据
XX = [Y,NEW];
rd(1:78) = randperm(78);
x_train = XX(rd(1:50),:);
x_test = XX(rd(51:end),:);
[RMSEC,Rc,RMSEP,Rp] = fitness_1(x_train,x_test);
W=100*eig1/sum(eig1);
figure(1)
pareto(W);                  %画出贡献率的直方图
bar(W);