clear all;
X1=xlsread('C:\Users\c430\Desktop\蛋白质含量测试\波长变量筛选\全谱\淀粉校正集.xlsx');%导入数据
X2=xlsread('C:\Users\c430\Desktop\蛋白质含量测试\波长变量筛选\全谱\淀粉预测集.xlsx');%导入数据
[Rc,RMSEC,beta,yc] = fitaaa(X1);
[Rp,RMSEP,yp] = fitbbb(X2,beta);