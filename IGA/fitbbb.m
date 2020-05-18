function [R,RMSEP,y] = fitbbb(x,X,L,beta)     %适应度函数


    %X=xlsread('C:\Users\c430\Desktop\GA\HY1.xlsx');%导入数据
    j = 1;
    %f = 7;
    for i=1:L
        if x(i) == 1
            XX(:,j:j+19) = X(:,20*i-18:20*i+1);
            j = j+20;
        end
    end
    YY = X(:,1);
    RMSEP = sqrt(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2)/28);
    y = XX*beta(2:end,:)+beta(1,:);
    R = sqrt(1-(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2))/(sum((YY-mean(YY)).^2)));

    z = R/(1+RMSEP);
    %disp(z)
end