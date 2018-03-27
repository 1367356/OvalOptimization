function [Array]=ArrayGroup(Population_InitA,Population_InitB,circle_num,group_i,element_space)
data=zeros(circle_num+1,2);                                                                          %
allPosition=[];
%��Population_Init��group_i����Ⱥ�Ƶ�data��
Perimeter=zeros(circle_num+1,1);
Circle_element_num=zeros(circle_num+1,1);
Circle_element_num(1,1)=1;

%Element_alpha=zeros(circle_num+1,1);
% Element_space_on_every_oval=zeros(circle_num+1,1); %ÿ����Բ�ϵ���Ԫ���

for i=1:circle_num+1
     data(i,1)=Population_InitA(i,group_i);
     data(i,2)=Population_InitB(i,group_i);
     a=data(i,1);
     b=data(i,2);
     Perimeter(i)=2*pi*b+4*(a-b);%��Բ�ܳ�
end

for circle_i=2:circle_num
         Circle_element_num(circle_i)=fix(Perimeter(circle_i)/(element_space*(1+circle_i*0.2)));  %ÿ��Բ�����õĸ���
end
Circle_element_num(circle_num+1)=100-sum(Circle_element_num(1:circle_num));

len=0;
% num_of_each_oval=zeros(circle_num+1,1);
% num_of_each_oval(1)=1;
allElementPointX=[];
allElementPointY=[];
for i=2:circle_num+1 %����,����һ��û��
    a=data(i,1);
    b=data(i,2);
%     minSpace=Element_space_on_every_oval(i);
    num_of_element=Circle_element_num(i);
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
plot(allElementPointX,allElementPointY,'*');
