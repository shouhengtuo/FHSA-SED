function [X,X2,NC,totaltime] = HS_Score3(data,dim_epi,HMS,max_iter,cx,cy)
%input--------------------------------------------------------------------
% data-----------------input dataset
% epi_dim--------------the epistasis dimension
% HMS--------------the size of harmony memory(HM)
% max_iter--------- the maximum evaluation times of disease model.
% (cx,cy) --------- the given interaction markers that have been swapped
% randomly in main.

%%-------------------------------------------------------------------------
% initial arguments
HMCR=0.95;        %% harmony memory consideration rate.
PAR=0.35;         %% pitch-adjustment rate
n=size(data,2);  
State=data(:,n); %% disease state

%% ---------------------------------------------------------------

SNPs=n-1;    %% the number of all SNP markers

%% ��ʼ�����ʼ�¼�� Initilize the table for trace the visited position


TD=(SNPs+1)*(SNPs-2)/2+1;
TB=zeros(5,TD) ;
head=1;
tail=TD;

FoundNumber=0;

for i=1:SNPs-1
    for j=i+1:SNPs
        sno=SNPs*(i-1)-i*(i-1)/2+j-i;
        TB(2,sno)=i;
        TB(3,sno)=j;
        if sno==1
          TB(4,sno)=-1; %% the first position which has not precursor ǰ��
        else
            TB(4,sno)=sno-1;%%  The precursor of this position in linked listǰ��
        end
        if sno==TD
            TB(5,sno)=-1;  %% the last position which has not successor ���һ��Ԫ�صĺ��λ��
        else
            TB(5,sno)=sno+1; %%  The successor of this position in linked list ���
        end
    end
end


X=zeros(HMS,dim_epi);

for i=1:HMS
    snp(1)=ceil(rand*SNPs);
    snp(2)=ceil(rand*SNPs);  
    snp=sort(snp);   
          while snp(2)==snp(1) || ismember(snp,X,'rows')
             snp(2)=ceil(rand*SNPs);
             snp=sort(snp);
          end
          
       sno= SNPs*(snp(1)-1)-snp(1)*(snp(1)-1)/2+snp(2)-snp(1);  %%���ɼ�¼���λ�� 
         
       TB(1,sno)=1; %% ����Ϊ�ѷ���
       %% ά�� ǰ�� ������ϵ
       %% ��¼ǰ��
        if TB(4,sno)==-1;  %% ��ʾ�Ǳ�ͷԪ��
            %% �鿴����Ԫ��
            next=TB(5,sno);            
            TB(4,next)=-1;
            head=next;
        elseif TB(5,sno)==-1;  %% ��ʾ�Ǳ�βԪ��
            pre=TB(4,sno);
            TB(5,pre)=-1;
            tail=pre;
        else  %% ά��ǰ�������ӹ�ϵ
            pre=TB(4,sno);  %% �ҵ�ǰ��
            next=TB(5,sno);
            TB(5,pre)=next;  %% ǰ���ĺ����sno��nextԪ��
            TB(4,next)=pre;  %% ��̵�ǰ����sno��preԪ��
        end        
        
    X(i,:)=snp;

    [Fit(i,1),Fit2(i,1)] = myMultiscore(data(:,X(i,:)),State);    
end

if sno<tail
   currPos=TB(5,sno);  %% ��Ԫ�ص���һ��Ԫ��λ��
else
    currPos=TB(4,sno); %% ��Ԫ�ص���һ��Ԫ��λ��
end

FoundNumber=FoundNumber+HMS; %% �Ѿ�����Ԫ�ظ���

X2=X;

LT=0;
NC=0;

%%-------------------------------------------------------------------------
tic;
while NC <= max_iter      
    [best1,idbest1]=min(Fit);
    [best2,idbest2]=min(Fit2);
    Xbest1=X(idbest1,:);
    Xbest2=X2(idbest2,:);
    BW=ceil((TD-FoundNumber)*rand/50);

     for i=1: 2  
        if rand<HMCR
            a=ceil(rand*HMS*2);
            if a<=HMS
                Xnew(i)=X(a,i); 
                if rand<PAR
                    %% Xnew(i)=Xnew(i)+round(2*(rand-0.5)*BW);
                   Xnew(i)=Xnew(i)+((-1)^Xnew(i))*abs(Xbest1(i)-X(ceil(rand*HMS),i));
                    Xnew(i)=max(min(Xnew(i),SNPs),1);
                end   
            else
%                 Xnew(i)=TB(3,currPos);
%                 Xnew(i)=ceil(rand*SNPs);
                 Xnew(i)=X2(a-HMS,i);
                 if rand<PAR
                     %Xnew(i)=Xnew(i)+round(2*(rand-0.5)*BW);
                    %%Xnew(i)=Xnew(i)+round(normrnd(0,10));
                    Xnew(i)=Xnew(i)+((-1)^Xnew(i))*abs(Xbest2(i)-X2(ceil(rand*HMS),i));
                    Xnew(i)=max(min(Xnew(i),SNPs),1);
                end 
            end
                
        else
            Xnew(i)=ceil(rand*SNPs);
        end
     end
    while Xnew(1)==Xnew(2)
        Xnew(2)=ceil(rand*SNPs);
    end
    Xnew=sort(Xnew);
       
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sno= SNPs*(Xnew(1)-1)-Xnew(1)*(Xnew(1)-1)/2+Xnew(2)-Xnew(1);  %%���ɼ�¼���λ�� 
         
     if  TB(1,sno)==1 %% �ѷ���
         %% ��head �� tail֮�������һ��
%              Rn=ceil(rand*(TD-FoundNumber));      
%              q=head;
%              p=TB(5,head);       
%              rn=1;
%              while rn<Rn && p<tail
%                  q=p;                
%                  p=TB(5,p);                
%                 rn=rn+1;
%              end
         %%  ---------�ֲ�����------------------------------------------------   
             Rn=BW;
             
             q=currPos;
            
             if rand<0.5  %% ���̽��
                 p=q;
                 rn=1;
                 while rn<Rn && p<tail
                     q=p;
                     p=TB(5,p);
                     rn=rn+1;
                 end
             else
                 p=q; %% ��ǰ̽��
                 rn=1;
                 while rn<Rn && p>head
                     q=p;
                     p=TB(4,p);
                     rn=rn+1;
                 end                 
             end
             TB(1,p)=1;
        %% -------------------------------------------------------------------
     else
         % Xnew(1)=TB(2,sno); Xnew(2)=TB(3,sno);
          TB(1,sno)=1;
          p=sno;
     end
     
     Xnew(1)=TB(2,p); Xnew(2)=TB(3,p);   
     FoundNumber=FoundNumber+1;
       %% ά�� P�� ǰ�� ������ϵ
       %% ��¼ǰ��
        if TB(4,p)==-1;  %% ��ʾ�Ǳ�ͷԪ��
            %% �鿴����Ԫ��
            next=TB(5,p);  
            
                TB(4,next)=-1;
                head=next;
            
        elseif TB(5,p)==-1;  %% ��ʾ�Ǳ�βԪ��
            pre=TB(4,p);
            TB(5,pre)=-1;
            tail=pre;
        else  %% ά��ǰ�������ӹ�ϵ
            pre=TB(4,p);  %% �ҵ�ǰ��
            next=TB(5,p);
            TB(5,pre)=next;  %% ǰ���ĺ����sno��nextԪ��
            TB(4,next)=pre;  %% ��̵�ǰ����sno��preԪ��
        end
        
        if p<tail
           currPos=TB(5,p);  %% ��Ԫ�ص���һ��Ԫ��λ��
        else
            currPos=TB(4,p); %% ��Ԫ�ص���һ��Ԫ��λ��
        end
%         p
%         currPos
        
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    Xnew=sort(Xnew);

  [score,score2] = myMultiscore(data(:,Xnew),State);
  
%   if Xnew==[cx,cy]
%         fprintf('--- input ---\n');
%     end
  
   [fworst,idworst]=max(Fit(:,1));
    [fworst2,idworst2]=max(Fit2(:,1));
    
        if score<=fworst || score2<=fworst2
            if  score<=fworst
                Fit(idworst,1)=score;           
                X(idworst,:)=Xnew; 
            end
            if score2<=fworst2
               Fit2(idworst2,1)=score2;           
               X2(idworst2,:)=Xnew; 
            end
        else
            LT=LT+1;
            if LT>SNPs
                Fit(idworst,1)=score;             
                X(idworst,:)=Xnew; 
                
                 Fit2(idworst2,1)=score2;             
                X2(idworst2,:)=Xnew; 
                
                LT=0;
            end
        end
    
    
    NC=NC+1;

 
   if ismember([cx,cy],X,'rows') ||  ismember([cx,cy],X2,'rows')

        break;
    end
end




totaltime=toc;
% 





