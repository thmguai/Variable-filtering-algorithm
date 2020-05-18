function [RMSEC,Rc,RMSEP,Rp] = fitness_1(X1,X2)     %适应度函数
    f = 7;
    XX = X1(:,2:end);
    YY = X1(:,1);
    [xl,yl,xs,ys,beta,pctvar,mse]=plsregress(XX,YY,f);%对xr和Y进行pls回归
    RMSEC = sqrt(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2)/50);
    Rc = sqrt(1-(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2))/(sum((YY-mean(YY)).^2)));
    
    XX = X2(:,2:end);
    YY = X2(:,1);
    RMSEP = sqrt(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2)/28);
    Rp = sqrt(1-(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2))/(sum((YY-mean(YY)).^2)));
end