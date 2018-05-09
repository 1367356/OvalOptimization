function [Array,allElementPoint]=ArrayGroup(Population_InitA,Population_InitB,circle_num,group_i,element_space,element_sum)
data=zeros(circle_num+1,2);                                                                          %
allPosition=[];
%��Population_Init��group_i����Ⱥ�Ƶ�data��
Perimeter=zeros(circle_num+1,1);
Circle_element_num=zeros(circle_num+1,1);
Circle_element_num(1,1)=1;

%Element_alpha=zeros(circle_num+1,1);
% Element_space_on_every_oval=zeros(circle_num+1,1); %ÿ����Բ�ϵ���Ԫ���
%allElementPoint={};
for i=1:circle_num+1
     data(i,1)=Population_InitA(i,group_i);
     data(i,2)=Population_InitB(i,group_i);
     a=data(i,1);
     b=data(i,2);
     Perimeter(i)=2*pi*b+4*(a-b);%��Բ�ܳ�
end

a =     -0.1538 ;% (-0.2078, -0.09973)
b =    0.004088  ;%(0.002347, 0.005828)
c =      0.5607  ;%(0.4597, 0.6618)

for circle_i=2:circle_num
         p=save_rate_of_radius(data(i,2),a,b,c);
           Circle_element_num(circle_i)=fix(p*(Perimeter(circle_i)/(element_space)));  %ÿ��Բ�����õĸ���
end
Circle_element_num(circle_num+1)=element_sum-sum(Circle_element_num(1:circle_num));  %����Ϊ��ֵ�����Ϊ��ֵ�����潫���������Ǳ�Խ�硣

len=0;
% num_of_each_oval=zeros(circle_num+1,1);
% num_of_each_oval(1)=1;
allElementPointX=[];
allElementPointY=[];
for i=2:circle_num+1 %����,����һ��û��
    a=data(i,1);
    b=data(i,2);
%     minSpace=Element_space_on_every_oval(i);
    num_of_element=Circle_element_num(i)
    [position,elementPointX,elementPointY] = array_calculate_bak(a,b,num_of_element);  %����������a,b��ֵ����Բ����Ԫ�ĸ����������Բ����Ԫ�ļ������ֱ�������λ��
    
    for index=1:length(position)
        allPosition(len+index,1)=position(index,1);
        allPosition(len+index,2)=position(index,2);
        allElementPointX(len+index)=elementPointX(index);
        allElementPointY(len+index)=elementPointY(index);
    end
    len=len+length(position);
end

Array=allPosition;  %���һ�������������Բ�ļ�����ֵ�ۼ�
%plot(allElementPointX,allElementPointY,'*');
allElementPoint={allElementPointX,allElementPointY};
