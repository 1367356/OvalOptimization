function [position,elementPointX,elementPointY] = array_calculate_bak(a,b,num_of_element)  
%��������
%a,����뾶
%b,����뾶
% num_of_element  ��Բ��Ҫ���õ���Ԫ����
%position,��Բ��������Ԫ�ļ�����λ��
%k ����Ϊa������Ϊb����Բ����Ԫ�ĸ���Ϊ k
%elementPointX,��Ԫ��ֱ������ϵ�µ�X�������λ��
%elementPointY����Ԫ��ֱ������ϵ�µ�Y�������λ��
% a=4,b=3,num_of_element=10;
minSpace=(2*pi*b+4*(a-b))/num_of_element;
h=0.000001; %����
data=0:h:a;  
x1=fliplr(data);  %��data���з�ת
y1=sqrt(b^2-b^2*x1.^2/a^2); %��Բ��ʽ
x2=-data;
y2=sqrt(b^2-b^2*x2.^2/a^2);
x3=-fliplr(data);
y3=-sqrt(b^2-b^2*x3.^2/a^2);
x4=data;
y4=-sqrt(b^2-b^2*x4.^2/a^2);
%y4=-fliplr(y1);
x=[x1,x2,x3,x4];
y=[y1,y2,y3,y4];
n=length(x);  %����
lamda=minSpace;   %lamda����
i=1;  %��i����
k=1;   %û�����ʼ�㣬���յ���k=num_of_element+1,
% e=0.0005;   %���,Ҫ�ǲ�����n��������ᵼ��ĳ��������֮�󣬾�һֱû�к��ʵĵ�
sumlength=0;
elementPointX=[];  
elementPointY=[];
%while (sum<2*pi*b+4*(a-b))
elementPointX(k)=a;
elementPointY(k)=0;
while(1<2)
    while (i<n)  %n��ĸ���
        m=sqrt( (x(i+1)-x(i))^2+(y(i+1)-y(i))^2);  %����֮��ľ���,mͻȻһ�±�ĺܴ󣬲���lamda��Χ����
        sumlength=sumlength+m;
        if((lamda-lamda*0.00001)<sumlength && sumlength<=(lamda+lamda*0.00001))  %��i����Χ��ʱ������
               k=k+1;
    %            v(k)=i;        
               sumlength=0;
               elementPointX(k)=x(i+1);
               elementPointY(k)=y(i+1);
               i = i + 1;
         else
                i = i + 1;
         end
    end 
    if(lamda-lamda*0.05<sumlength)
        break;
    else
        lamda=lamda-lamda*0.001;
        i=1;  
        k=1;   
        sumlength=0;
    end
end
%plot(elementPointX,elementPointY);  //ֱ������

angle=0;%�Ƕ�
radius=0;%�뾶
len=int8(length(elementPointX));
position=zeros(len,2); %�����Բ����ĵ�İ뾶�ͽǶ�

for index=1:len
    radius=sqrt(elementPointX(index)^2+elementPointY(index)^2);  %�뾶
    if (elementPointX(index)>0 && elementPointY(index)>0)
        angle=atan(elementPointY(index)/elementPointX(index));
    else if (elementPointX(index)<0 && elementPointY(index)>0)
        angle=pi/2+abs(atan(elementPointX(index)/elementPointY(index)));
    else if (elementPointX(index)<0 && elementPointY(index)<0)
        angle=pi+abs(atan(elementPointY(index)/elementPointX(index)));
    else if (elementPointX(index)>0 && elementPointY(index)<0) 
        angle=pi*3/2+abs(atan(elementPointX(index)/elementPointY(index)));
        end
        end
        end
    end
    position(index,1)=radius;  %�뾶
    position(index,2)=angle;   %�Ƕ�
end
% position


