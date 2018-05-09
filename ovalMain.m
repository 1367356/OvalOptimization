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
element_sum=120;
element_space=0.4*lamda;
%-----遗传算法参数
genetic_num=200;%遗传代数  
group_num=20;%种群数
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

%  画一个半径均分满阵圆,计算结果，作为对比。
% for group_i=1:group_num    %初始化种群半径，总共初始化了10个种群半径，从这10个半径组中，交叉变异，选择出最好的。
%     for n=1:circle_num+1
%         if n==1
%             Population_Init(n,group_i)=0;
%         else
%             Population_Init(n,group_i)=(L/12)*(n-1); %R(1)=0
%         end
%     end
% end

%------------------------------------------------------------------------------------------------
Population_Init=sort(Population_Init);%%%%按列从小到大排列
Population_InitA=Population_Init*MULTI;  %长轴种群
Population_InitB=Population_Init; %短轴种群

Tem_rsll=zeros(genetic_num,1); %临时rsll最小值
Rsll_it=zeros(group_num,1);
Array=[];
optimal_path={};
for genetic_i=1:genetic_num
    genetic_i
    Rsll_it=zeros(group_num,1);%每一代的各个个体峰值旁瓣电平存储空间
    for m=1:group_num
        [Array,allElementPoint]=ArrayGroup(Population_InitA,Population_InitB,circle_num,m,element_space,element_sum);
         Rsll_it(m)=ovalRSLLofCircle(Array);   %求阵列方向图以及峰值旁瓣电平
    end
    min(Rsll_it)
    Tem_rsll(genetic_i)=min(Rsll_it);%找出最小峰值旁瓣电平
    [aftersort_Rsll_it,Index]=sort(Rsll_it);
%     optimal_path{genetic_i}=Population_InitB(:,Index(1))';
     %----生成下一代种群的随机数
     if genetic_i<genetic_num
         [Population_next]=nextgroup_bak(Population_InitB,Rsll_it,group_num,circle_num,L,element_space);
         Population_InitB=Population_next;%%%%按列从小到大排列
         Population_InitA=Population_InitB*MULTI;
     end
     rsll=min(Tem_rsll);
end

sort(Tem_rsll)
% [aftersort_Rsll_it,Index]=sort(Tem_rsll);
[aftersort_Rsll_it,Index_Rsll]=sort(Rsll_it);
% global_optimal_path=optimal_path{Index(1)};
% allElementPointGenetici{Index(1)}  %最优阵元位置

% sort(Tem_rsll);
% [aftersort_Rsll_it,Index]=sort(Tem_rsll);
%global_optimal_path=optimal_path{Index(1)}
%-----导出最优阵列
n=Index_Rsll(1)
m=1;
[Array,allElementPoint]=ArrayGroup(Population_InitA,Population_InitB,circle_num,m,element_space,element_sum);  %Array是阵列中所有阵元组成的极坐标值。

Population_InitB(:,m)   %打印出最优圆环的半径

Bestallarray=Array;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%优化结果展示%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----导出方向图
N_sam=1024;  %采样点数
u=-1:2/N_sam:1-2/N_sam;  %u=sin(sita)*cos(fai)，U轴
v=-1:2/N_sam:1-2/N_sam;  %v=sin(sita)*sin(fai)，V轴
[X,Y]=meshgrid(u,v);
FF=zeros(N_sam,N_sam);%采样矩阵
lamda=1;%波长

fun_x=find(Bestallarray~=0);     
%得到的结果为一列向量，有阵元的位置索引值，维数一般比quartered_matrix低，就是把
%没有阵元的位置去掉
q_position=[real(Bestallarray(fun_x)) imag((Bestallarray(fun_x)))]; 
%得到的结果为二维列向量，第一列幅值，第二列相位，行数和fun_x一样    
%-----下面这段代码就是方向图的计算 公式（2）和（3）只是转换为在直角坐标下面计算
for n=1:N_sam
    for m=1:N_sam
        if abs(v(m))<=sqrt(1-u(n)^2)% 数学推导合理性保证，限定在单位圆内部
           temp=0;
           for a=1:length(fun_x)%把所有阵元的用上的，所以应该是整个圆平面  
               temp=temp+exp(j*2*pi*q_position(a,1)*(cos(q_position(a,2))...
                       *u(n)+sin(q_position(a,2))*v(m)));
           end
           FF(n,m)=temp+1;%加上圆心        
       else 
       end
   end
end

FF(find(FF==0))=eps; %eps是非常小的一个数
ff=20*log10(abs(FF)/max(max(abs(FF))));%%归一化
bottom=-40;%设置底平台电平
ff(find(ff<=-40))=bottom;%最低电平设置为-80db
%when fai=0 ---ff(m,:) ---u axis   
%when fai=90 ---ff(:,m) ---v axis  when fai=45----u=v
m=ceil(find(ff==max(max(ff)))/N_sam);%在一个100*100的矩阵中找到最大值的位置
fai0=ff(m,:); fai90=ff(:,m);  

%-----uv图
figure
plot(u,fai0,'-b');
hold on
plot(v,fai90,'--r');
legend('u=0','v=0');
%ylabel('F(u,v)/dB');
%xlabel('μ,ν');
figure_FontSize=14;
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
xlabel('u,v');
ylabel('Radiation pattern (dB)');
%-----方向图
figure
mesh(X,Y,ff);%画出曲面图
shading interp;
colormap(gray);
% plot(Array(:,1)',Array(:,2)','*');
%-----导出最优阵列
% [Array]=ArrayGroup(Population_Init,circle_num,1,element_space);
% Bestallarray=Array;
xlabel('u=sin\theta cos φ');
ylabel('u=sin\theta sin φ');
zlabel('Radiation pattern(dB)');

%-----收敛图
figure
plot(Tem_rsll);
allElementPoint
xlabel('generation');
ylabel('PSLL(dB)');

allElementPoint{1}(element_sum)=0;
allElementPoint{2}(element_sum)=0;
%阵元分布图
figure
plot(allElementPoint{1},allElementPoint{2},'bo');
xlabel('x(λ)');
ylabel('y(λ)');
%title('XLamda YLamda title');

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