%������ͬ�ĳ�����ֵ������array_calculate����������λ��
%��Բ����������
%��Ⱥ�Ż����ó���Ϊ����Ķ��ٱ����ɡ�
clc
clear
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamda=1;
circle_num=6;
H=circle_num;
L=9.4;
MULTI=1.2;  %����/����
element_sum=142;
element_space=0.5*lamda;
%-----�Ŵ��㷨����
genetic_num=3;%�Ŵ�����
group_num=10;%��Ⱥ��
Population_Init=zeros(circle_num+1,group_num);
%---------------------------------------------------------------------------------------

basic_distance=element_space:element_space:circle_num*element_space;
for group_i=1:group_num    %��ʼ����Ⱥ�뾶���ܹ���ʼ����10����Ⱥ�뾶������10���뾶���У�������죬ѡ�����õġ�
    exprnd_num=exprnd(0.5,1,circle_num);%����ָ������
    redundance=L/2-circle_num*element_space;    %����Բ��֮����Ԫ��������£�ʣ�µ�����
    normalization_num=sort(exprnd_num/max(exprnd_num)*redundance);  %��һ��*1.7
    group_i_radius=normalization_num+basic_distance;         %����ԭ���İ뾶,���ɵİ뾶��Ϊa����Ϊa����С���룬
    for n=1:circle_num+1
        if n==1
           Population_Init(n,group_i)=0;
        else
            Population_Init(n,group_i)=group_i_radius(n-1); %R(1)=0
        end
    end
end

%------------------------------------------------------------------------------------------------
Population_Init=sort(Population_Init);%%%%���д�С��������
Population_InitA=Population_Init*MULTI;  %������Ⱥ
Population_InitB=Population_Init;  %������Ⱥ

Tem_rsll=zeros(genetic_num,1);%��ʱrsll��Сֵ
Array=[];

for genetic_i=1:genetic_num
    genetic_i
    Rsll_it=zeros(group_num,1);%ÿһ���ĸ��������ֵ�԰��ƽ�洢�ռ�
    for m=1:group_num
        [Array]=ArrayGroup(Population_InitA,Population_InitB,circle_num,m,element_space);
         Rsll_it(m)=ovalRSLLofCircle(Array);   %�����з���ͼ�Լ���ֵ�԰��ƽ
    end
     Tem_rsll(genetic_i)=min(Rsll_it);%�ҳ���С��ֵ�԰��ƽ
     %----������һ����Ⱥ�������
     [Population_next]=nextgroup_bak(Population_InitB,Rsll_it,group_num,circle_num,L,element_space);
     Population_InitB=sort(Population_next);%%%%���д�С��������
     Population_InitA=Population_InitB*MULTI;
     rsll=min(Rsll_it);
     %%%�ӽ�����
end
% plot(Array(:,1)',Array(:,2)','*');
%-----������������
% [Array]=ArrayGroup(Population_Init,circle_num,1,element_space);
% Bestallarray=Array;
Tem_rsll
rsll=min(Tem_rsll)


%-----����ͼ
figure

plot(Tem_rsll);

%------------------------------------------------------------------------------------------
% 
% for i=1:length(data) %����
%     a=data(i,1);
%     b=data(i,2);
%     [position] = array_calculate(a,b,minSpace);  %����������a,b��ֵ��ÿ����Բ�ϵ���С��Ԫ��࣬�����Բ����Ԫ�ļ������ֱ�������λ��
%     for index=1:length(position)
%         allPosition(len+index,1)=position(index,1);
%         allPosition(len+index,2)=position(index,2);
%     end
%     len=len+length(position);
%     %allPosition{i}=position;
% end
% 
% %-----------------------------------------------------------------------------
% allPosition;  %���һ�������������Բ�ļ�����ֵ�ۼ�
% rsll=ovalRSLLofCircle(allPosition);
% 
% plot(allelementPointX1,allelementPointY1,'*');
% 
% for i=1:length(eachOvalPosition)
%     len=length(eachOvalPosition{i})/2;
%     for j=1:len
%         if(j<len)
%             line([eachOvalPosition{i}(j),eachOvalPosition{i}(j+1)],[eachOvalPosition{i}(len+j),eachOvalPosition{i}(len+j+1)]);
%         else
%             line([eachOvalPosition{i}(j),eachOvalPosition{i}(1)],[eachOvalPosition{i}(len+j),eachOvalPosition{i}(len+1)]);
%         end
%     end
% end
% 
% maxValue=max(data);
% maxValue=maxValue(1,1);
%axis([-maxValue maxValue -maxValue maxValue]);