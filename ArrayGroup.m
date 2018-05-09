function [Array,allElementPoint]=ArrayGroup(Population_InitA,Population_InitB,circle_num,group_i,element_space,element_sum)
data=zeros(circle_num+1,2);                                                                          %
allPosition=[];
%将Population_Init第group_i个种群移到data中
Perimeter=zeros(circle_num+1,1);
Circle_element_num=zeros(circle_num+1,1);
Circle_element_num(1,1)=1;

%Element_alpha=zeros(circle_num+1,1);
% Element_space_on_every_oval=zeros(circle_num+1,1); %每个椭圆上的阵元间距
%allElementPoint={};
for i=1:circle_num+1
     data(i,1)=Population_InitA(i,group_i);
     data(i,2)=Population_InitB(i,group_i);
     a=data(i,1);
     b=data(i,2);
     Perimeter(i)=2*pi*b+4*(a-b);%椭圆周长
end

a =     -0.1538 ;% (-0.2078, -0.09973)
b =    0.004088  ;%(0.002347, 0.005828)
c =      0.5607  ;%(0.4597, 0.6618)

for circle_i=2:circle_num
         p=save_rate_of_radius(data(i,2),a,b,c);
           Circle_element_num(circle_i)=fix(p*(Perimeter(circle_i)/(element_space)));  %每个圆环放置的个数
end
Circle_element_num(circle_num+1)=element_sum-sum(Circle_element_num(1:circle_num));  %不能为负值，如果为负值，下面将会出现数组角标越界。

len=0;
% num_of_each_oval=zeros(circle_num+1,1);
% num_of_each_oval(1)=1;
allElementPointX=[];
allElementPointY=[];
for i=2:circle_num+1 %长轴,中心一点没算
    a=data(i,1);
    b=data(i,2);
%     minSpace=Element_space_on_every_oval(i);
    num_of_element=Circle_element_num(i)
    [position,elementPointX,elementPointY] = array_calculate_bak(a,b,num_of_element);  %给定长短轴a,b的值和椭圆上阵元的个数，求出椭圆中阵元的极坐标和直角坐标的位置
    
    for index=1:length(position)
        allPosition(len+index,1)=position(index,1);
        allPosition(len+index,2)=position(index,2);
        allElementPointX(len+index)=elementPointX(index);
        allElementPointY(len+index)=elementPointY(index);
    end
    len=len+length(position);
end

Array=allPosition;  %求得一个个体的所有椭圆的极坐标值累加
%plot(allElementPointX,allElementPointY,'*');
allElementPoint={allElementPointX,allElementPointY};
