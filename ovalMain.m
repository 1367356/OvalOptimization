%给定不同的长短轴值，调用array_calculate，生成坐标位置
%椭圆阵列主函数
%种群优化设置长轴为短轴的多少倍即可。
clc
clear
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamda=1;
circle_num=6;
H=circle_num;
L=9.4;
MULTI=1.2;  %长轴/短轴
element_sum=142;
element_space=0.5*lamda;
%-----遗传算法参数
genetic_num=3;%遗传代数
group_num=10;%种群数
Population_Init=zeros(circle_num+1,group_num);
%---------------------------------------------------------------------------------------

basic_distance=element_space:element_space:circle_num*element_space;
for group_i=1:group_num    %初始化种群半径，总共初始化了10个种群半径，从这10个半径组中，交叉变异，选择出最好的。
    exprnd_num=exprnd(0.5,1,circle_num);%产生指数序列
    redundance=L/2-circle_num*element_space;    %满足圆环之间阵元间距的情况下，剩下的冗余
    normalization_num=sort(exprnd_num/max(exprnd_num)*redundance);  %归一化*1.7
    group_i_radius=normalization_num+basic_distance;         %加上原来的半径,生成的半径作为a，因为a有最小距离，
    for n=1:circle_num+1
        if n==1
           Population_Init(n,group_i)=0;
        else
            Population_Init(n,group_i)=group_i_radius(n-1); %R(1)=0
        end
    end
end

%------------------------------------------------------------------------------------------------
Population_Init=sort(Population_Init);%%%%按列从小到大排列
Population_InitA=Population_Init*MULTI;  %长轴种群
Population_InitB=Population_Init;  %短轴种群

Tem_rsll=zeros(genetic_num,1);%临时rsll最小值
Array=[];

for genetic_i=1:genetic_num
    genetic_i
    Rsll_it=zeros(group_num,1);%每一代的各个个体峰值旁瓣电平存储空间
    for m=1:group_num
        [Array]=ArrayGroup(Population_InitA,Population_InitB,circle_num,m,element_space);
         Rsll_it(m)=ovalRSLLofCircle(Array);   %求阵列方向图以及峰值旁瓣电平
    end
     Tem_rsll(genetic_i)=min(Rsll_it);%找出最小峰值旁瓣电平
     %----生成下一代种群的随机数
     [Population_next]=nextgroup_bak(Population_InitB,Rsll_it,group_num,circle_num,L,element_space);
     Population_InitB=sort(Population_next);%%%%按列从小到大排列
     Population_InitA=Population_InitB*MULTI;
     rsll=min(Rsll_it);
     %%%加进度条
end
% plot(Array(:,1)',Array(:,2)','*');
%-----导出最优阵列
% [Array]=ArrayGroup(Population_Init,circle_num,1,element_space);
% Bestallarray=Array;
Tem_rsll
rsll=min(Tem_rsll)


%-----收敛图
figure

plot(Tem_rsll);

%------------------------------------------------------------------------------------------
% 
% for i=1:length(data) %长轴
%     a=data(i,1);
%     b=data(i,2);
%     [position] = array_calculate(a,b,minSpace);  %给定长短轴a,b的值，每个椭圆上的最小阵元间距，求出椭圆中阵元的极坐标和直角坐标的位置
%     for index=1:length(position)
%         allPosition(len+index,1)=position(index,1);
%         allPosition(len+index,2)=position(index,2);
%     end
%     len=len+length(position);
%     %allPosition{i}=position;
% end
% 
% %-----------------------------------------------------------------------------
% allPosition;  %求得一个个体的所有椭圆的极坐标值累加
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