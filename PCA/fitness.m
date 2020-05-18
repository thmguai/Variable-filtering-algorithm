function z = fitness(i,X)     %适应度函数
    f = 7;
    j = 20*i-18;
    %z = 20*i+1;
    XX = X(:,j:j+19);
    YY = X(:,1);
    [xl,yl,xs,ys,beta,pctvar,mse]=plsregress(XX,YY,f);%对xr和Y进行pls回归
    RMSEC = sqrt(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2)/50);

    R = sqrt(1-(sum((YY-(XX*beta(2:end,:)+beta(1,:))).^2))/(sum((YY-mean(YY)).^2)));

    z = R/(1+RMSEC);
   
end