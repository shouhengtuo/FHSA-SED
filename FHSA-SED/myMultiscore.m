function [score1,score2] = myMultiscore(snp_com,state)
[xrow,~] = size(snp_com);
%% K2 score 

subs = snp_com+1;
sample = accumarray(subs,ones(xrow,1)); %% 统计各个基因型组合 对应 的样本总数
disease = accumarray(subs,state) ;   %% 统计基因型组合 对应 的疾病样本数
control = sample-disease;    %% 计算基因型组合 对应 的control 样本数

sample(4,4) = 0;
disease(4,4) = 0;
control(4,4) = 0;
z=0;


for i = 1:3
    for j = 1:3
        y=My_factorial(sample(i,j)+1);
        r=My_factorial(disease(i,j))+My_factorial(control(i,j));
        z=z+(r-y);        
    end
end
score1=abs(z);

%% GINI score

P=sample(1:3,1:3)./xrow;  %% 每种基因型组合在样本中的频率。
Pcase=zeros(3,3);
Pcontrol=zeros(3,3);
for i=1:3
    for j=1:3
       if sample(i,j)~=0
           Pcase(i,j)=disease(i,j)./sample(i,j);
           Pcontrol(i,j)=control(i,j)./sample(i,j);
       end
    end
end

score2=sum(sum(P.*(1-(Pcase.^2+Pcontrol.^2))));


%score2=P1*(1-sum((disease/N1).^2)) + P2*(1-sum((control/N2).^2))

%% f is function used to calculate log form factorial
    function f=My_factorial(e)
        f=0;
        if e>0
            for o=1:e
                f=f+log(o);
            end
        end
    end
end


