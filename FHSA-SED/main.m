
%clc;
clear;
samplesize = 2000;
data_num = 100;
dim = 100;  %% number of SNPs in the mode data
Bp = 0.01;
dim_epi = 2;  %% two-way
HMS = 10;     %% Harmony memory size

max_iter =4980;  %% In HS algorithm, only one new harmony (solution) will be evaluated at each iteration.
comb = nchoosek(dim,dim_epi);
pvalue = max(min(Bp/comb,1e-10),1e-16);  %% P-value Bonferroni correction 

datafolder='DATA\';%%store results data   数据存放
ResultName='Result';
PowerName='Power';



warning('off');

Runtimes=2;  %% run times : 

 %% known functional marker (99,100)
Cx=99;
Cy=100;

slash = '\';

 power=[{'K2'},{'GINI'},{'K2-GINI'},{'G-test'},{'MOFHSST'},{'Evaluation Times'},{'Mean time'}];
           result=[{'True positive rate (TPR,Sensitivity, Recall)'},{'True negative rate(TNR, Specificity, SPC) '},{'Positive predictive value (PPV, Precision) '},{'False discovery rate(FDR) '},{'Accuracy (ACC) '},{'F1-Score'},{'Negative predictive value (NPV) '},{'FPR'},{'FNR'}];
           dataFile=strcat(datafolder,PowerName);
            dataFile=strcat(dataFile,'.xls'); 
           xlswrite(dataFile,power,'b1:h1');
           
           dataFile=strcat(datafolder,ResultName);
             dataFile=strcat(dataFile,'.xls');
             xlswrite(dataFile,result,'b1:J1');

for MODEL=1:1
        switch(MODEL)
            
            case 1
                model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.2';
                 root='.\Mode data\h2=0.02,Pd=0.1,MAF=0.2\1000CASE_EDM-1\1000CASE_EDM-1';
                
           %% add you model data that can be generated using GAMETES_2.1
            case 2
                %% %% add you model data that can be generated using GAMETES_2.1
            case 3 % ...   
        end

%          STIME=zeros(1,100*Runtimes);  %% 记录运行时间，绘制合图
%          SITER=zeros(1,100*Runtimes);
        for T = 1:Runtimes  %% 

            TP=0;
            FP=0;
            TN=0;
            FN=0;
            Hc=0;
            Kc=0;
            Gc=0;
            Lc=0;
            ac = 0;
            ac2=0;
            ac_inter = 0;
            fprintf(' model:%s, parameter:%s\n',  model, parameter);
            for d = 0:data_num-1

                data_iteration = d;
                fprintf('data-%2d: ',d+1); 
      %% get the file path of modedata 
               flag1=0; flag2=0;
               if d+1<10
                   STR=['_','00',num2str(d+1),'.txt'];
               elseif d+1<100
                   STR=['_','0',num2str(d+1),'.txt'];
               else
                   STR=['_',num2str(d+1),'.txt'];
               end
               filename=[root,STR];

       %%read data      %     a = dlmread(filename,'\t',1,0);
                a = dlmread(filename,'\t',1,0);
                data = a(:,[1:dim+1]);
      %% swap the functional markers randomly for avoiding preference to the position of SNP markers
                 cx=ceil(rand*dim);
                 cy=ceil(rand*dim);
                 while(cx>=cy)
                     cy=ceil(rand*dim); cx=ceil(rand*dim);
                 end                     
        
                temp=data(:,[Cx,Cy]);
                data(:,[Cx,Cy])=data(:,[cx cy]);
                data(:,[cx cy])=temp;

       %%   Harmony search the functional markers

                    [filter_snps, filter_snps2, iter(d+1),Time(d+1)] = HS_Score3(data,dim_epi,HMS,max_iter,cx,cy);
                    
        %             [score1,score2] = myMultiscore(data(:,[cx,cy]),data(:,dim+1))
        %% computer the power of K2SCORE AND GINI-SCORE       
                if ismember([cx,cy],filter_snps,'rows')
                    Kc=Kc+1;flag1=1;
                end
                if ismember([cx,cy],filter_snps2,'rows')
                    Gc=Gc+1;flag2=1;
                end

                if flag1==1||flag2==1
                    Hc=Hc+1;
                end
                k = d+1;

         %% UNIQUE THE DATA %%%%%%%%%%%%%%%%%%%%%%%%
                filter_snps2=[filter_snps;filter_snps2];
                filter_snps2=unique(filter_snps2,'rows');
                LS=length(filter_snps2(:,1));
                Score=zeros(LS,1);
                state=data(:,end);
                xrow=samplesize;      
   
         %% G-TEST
             for i=1:LS

                [score(i),Score(i)] = Gtest_score3(data(:,filter_snps2(i,:)),state);
             end
             
             [Score,sind]=sort(Score);
             filter_snps2=filter_snps2(sind,:);
             
             for i=1:LS
               snp_com=filter_snps2(i,:);
%                mdl=fitglm(data(:,snp_com),state,modelSpec,'Distribution','binomial','Link','logit'); 
%                   Cqtest=devianceTest(mdl);
                    
                   if Score(i)<pvalue
                       if  cx==filter_snps2(i,1) && cy==filter_snps2(i,2)
                             Lc=Lc+1;
                             TP=TP+1;
                             fprintf('^_^ TP=true,1st stage rank:%3d ,  pvalue=%-8.2e;  MEs=%-6d,  Mean time=%-10.2f\n',i,Score(i),iter(d+1),Time(d+1));
                       else
                           FP=FP+1;
                          %% fprintf('xxx FP=true,(%d,%d)    pvalue=%-13.5e; \n',filter_snps2(i,1),filter_snps2(i,2),Score(i));
                       end
                   else
                       if  cx==filter_snps2(i,1) && cy==filter_snps2(i,2)                            
                           FN=FN+1;
                           fprintf('!!! FN=true, 1st stage rank:%3d ,  pvalue=%-8.2e;  MEs=%-6d,  Mean time=%-10.2f\n',i,Score(i),iter(d+1),Time(d+1));
                       else
                           TN=TN+1;
                       end
                   end
             end

        %        True positive(TP): Sick people correctly diagnosed as sick
        %         False positive (FP): Healthy people incorrectly identified as sick
        %         True negative(TN): Healthy people correctly identified as healthy
        %         False negative(FN): Sick people incorrectly identified as healthy


            end
            
            %% BEGINING THE POWER AND OTHER METRICS
            Power_K2=Kc/100;
            Power_Gini=Gc/100;
            Power_K2andGini=Hc/100;
            GTEST_POWER=Lc/100;
            FHSA_POWER =TP/100;
            Result(T,:)=[Power_K2, Power_Gini,Power_K2andGini,GTEST_POWER,FHSA_POWER, mean(iter),mean(Time)];
            
            STIME(100*(T-1)+1:T*100)=Time;
            SITER(100*(T-1)+1:T*100)=iter;
           
              TPR=TP/(TP+FN); %% sensitivity or ture positive rate (TPR)
          
               SPC=TN/(FP+TN); %% specificity(SPC) or true negative rate(TNR)
                    
             
            FPR=FP/(FP+TN);            %% * False positive rate (FPR)
            FNR=FN/(FN+TP);            %% * False Negative Rate
            TNR=TN/(TN+FP);
            
            %%PPV=TP/(TP+FP);% precision or positive predictive value (PPV)
            PPV=TPR/(TPR+FPR); 
            
            %% FDR=FP/(FP+TP);  %% 
            FDR=FPR/(FPR+TPR);
            
            ACC=(TPR+TNR)/(TPR+FNR+FPR+TNR); %% accuracy;
            
            F1=2*TPR/(2*TPR+FPR+FNR);       % F1 score is the harmonic mean of precision and sensitivity
           
            
            %% NPV =TN/(TN+FN);           %%  Negative predictive value (NPV)
              NPV=TNR/(TNR+FNR); %%NPV=(TN/(TN+FP))/( (TN/(TN+FP))+ (FN/(TP+FN)));  
            
            RESULT2(T,:)=[TPR,SPC,PPV,FDR,ACC,F1,NPV,FPR,FNR];
        end
          
         
        
        %% save Power data into dataFile  [[Power_K2, Power_Gini,Power_K2andGini,CHi_Square,Power_MOFHSST, mean(iter),mean(Time)];]
            dataFile=strcat(datafolder,PowerName);
            dataFile=strcat(dataFile,'.xls');
           
            ModelParameter= {strcat(strcat(model,':  '),parameter)}; 
              cellPosition=strcat('A',num2str(MODEL+1));
              cellPosition=strcat(cellPosition,strcat(':',cellPosition));
              
            xlswrite(dataFile,ModelParameter,cellPosition);
            
            cellPosition=strcat('B',num2str(MODEL+1));
            cellPosition=strcat(cellPosition,':H');
             cellPosition=strcat(cellPosition,num2str(MODEL+1));      
          
               xlswrite(dataFile,mean(Result,1),cellPosition); 
          
         %%  save Power data into dataFile  [TPR,SPC,PPV,FDR,ACC,F1,NPV,FPR,FNR]
             dataFile=strcat(datafolder,ResultName);
             dataFile=strcat(dataFile,'.xls');
             
              ModelParameter= {strcat(strcat(model,':  '),parameter)}; 
              cellPosition=strcat('A',num2str(MODEL+1));
              cellPosition=strcat(cellPosition,strcat(':',cellPosition));
              
            xlswrite(dataFile,ModelParameter,cellPosition);
            
            cellPosition=strcat('B',num2str(MODEL+1));
            cellPosition=strcat(cellPosition,':J');
            cellPosition=strcat(cellPosition,num2str(MODEL+1));          
           
            xlswrite(dataFile,mean(RESULT2,1),cellPosition);


       
end