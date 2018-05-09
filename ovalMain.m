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
element_sum=120;
element_space=0.4*lamda;
%-----�Ŵ��㷨����
genetic_num=200;%�Ŵ�����  
group_num=20;%��Ⱥ��
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

%  ��һ���뾶��������Բ,����������Ϊ�Աȡ�
% for group_i=1:group_num    %��ʼ����Ⱥ�뾶���ܹ���ʼ����10����Ⱥ�뾶������10���뾶���У�������죬ѡ�����õġ�
%     for n=1:circle_num+1
%         if n==1
%             Population_Init(n,group_i)=0;
%         else
%             Population_Init(n,group_i)=(L/12)*(n-1); %R(1)=0
%         end
%     end
% end

%------------------------------------------------------------------------------------------------
Population_Init=sort(Population_Init);%%%%���д�С��������
Population_InitA=Population_Init*MULTI;  %������Ⱥ
Population_InitB=Population_Init; %������Ⱥ

Tem_rsll=zeros(genetic_num,1); %��ʱrsll��Сֵ
Rsll_it=zeros(group_num,1);
Array=[];
optimal_path={};
for genetic_i=1:genetic_num
    genetic_i
    Rsll_it=zeros(group_num,1);%ÿһ���ĸ��������ֵ�԰��ƽ�洢�ռ�
    for m=1:group_num
        [Array,allElementPoint]=ArrayGroup(Population_InitA,Population_InitB,circle_num,m,element_space,element_sum);
         Rsll_it(m)=ovalRSLLofCircle(Array);   %�����з���ͼ�Լ���ֵ�԰��ƽ
    end
    min(Rsll_it)
    Tem_rsll(genetic_i)=min(Rsll_it);%�ҳ���С��ֵ�԰��ƽ
    [aftersort_Rsll_it,Index]=sort(Rsll_it);
%     optimal_path{genetic_i}=Population_InitB(:,Index(1))';
     %----������һ����Ⱥ�������
     if genetic_i<genetic_num
         [Population_next]=nextgroup_bak(Population_InitB,Rsll_it,group_num,circle_num,L,element_space);
         Population_InitB=Population_next;%%%%���д�С��������
         Population_InitA=Population_InitB*MULTI;
     end
     rsll=min(Tem_rsll);
end

sort(Tem_rsll)
% [aftersort_Rsll_it,Index]=sort(Tem_rsll);
[aftersort_Rsll_it,Index_Rsll]=sort(Rsll_it);
% global_optimal_path=optimal_path{Index(1)};
% allElementPointGenetici{Index(1)}  %������Ԫλ��

% sort(Tem_rsll);
% [aftersort_Rsll_it,Index]=sort(Tem_rsll);
%global_optimal_path=optimal_path{Index(1)}
%-----������������
n=Index_Rsll(1)
m=1;
[Array,allElementPoint]=ArrayGroup(Population_InitA,Population_InitB,circle_num,m,element_space,element_sum);  %Array��������������Ԫ��ɵļ�����ֵ��

Population_InitB(:,m)   %��ӡ������Բ���İ뾶

Bestallarray=Array;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�Ż����չʾ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----��������ͼ
N_sam=1024;  %��������
u=-1:2/N_sam:1-2/N_sam;  %u=sin(sita)*cos(fai)��U��
v=-1:2/N_sam:1-2/N_sam;  %v=sin(sita)*sin(fai)��V��
[X,Y]=meshgrid(u,v);
FF=zeros(N_sam,N_sam);%��������
lamda=1;%����

fun_x=find(Bestallarray~=0);     
%�õ��Ľ��Ϊһ������������Ԫ��λ������ֵ��ά��һ���quartered_matrix�ͣ����ǰ�
%û����Ԫ��λ��ȥ��
q_position=[real(Bestallarray(fun_x)) imag((Bestallarray(fun_x)))]; 
%�õ��Ľ��Ϊ��ά����������һ�з�ֵ���ڶ�����λ��������fun_xһ��    
%-----������δ�����Ƿ���ͼ�ļ��� ��ʽ��2���ͣ�3��ֻ��ת��Ϊ��ֱ�������������
for n=1:N_sam
    for m=1:N_sam
        if abs(v(m))<=sqrt(1-u(n)^2)% ��ѧ�Ƶ������Ա�֤���޶��ڵ�λԲ�ڲ�
           temp=0;
           for a=1:length(fun_x)%��������Ԫ�����ϵģ�����Ӧ��������Բƽ��  
               temp=temp+exp(j*2*pi*q_position(a,1)*(cos(q_position(a,2))...
                       *u(n)+sin(q_position(a,2))*v(m)));
           end
           FF(n,m)=temp+1;%����Բ��        
       else 
       end
   end
end

FF(find(FF==0))=eps; %eps�Ƿǳ�С��һ����
ff=20*log10(abs(FF)/max(max(abs(FF))));%%��һ��
bottom=-40;%���õ�ƽ̨��ƽ
ff(find(ff<=-40))=bottom;%��͵�ƽ����Ϊ-80db
%when fai=0 ---ff(m,:) ---u axis   
%when fai=90 ---ff(:,m) ---v axis  when fai=45----u=v
m=ceil(find(ff==max(max(ff)))/N_sam);%��һ��100*100�ľ������ҵ����ֵ��λ��
fai0=ff(m,:); fai90=ff(:,m);  

%-----uvͼ
figure
plot(u,fai0,'-b');
hold on
plot(v,fai90,'--r');
legend('u=0','v=0');
%ylabel('F(u,v)/dB');
%xlabel('��,��');
figure_FontSize=14;
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
xlabel('u,v');
ylabel('Radiation pattern (dB)');
%-----����ͼ
figure
mesh(X,Y,ff);%��������ͼ
shading interp;
colormap(gray);
% plot(Array(:,1)',Array(:,2)','*');
%-----������������
% [Array]=ArrayGroup(Population_Init,circle_num,1,element_space);
% Bestallarray=Array;
xlabel('u=sin\theta cos ��');
ylabel('u=sin\theta sin ��');
zlabel('Radiation pattern(dB)');

%-----����ͼ
figure
plot(Tem_rsll);
allElementPoint
xlabel('generation');
ylabel('PSLL(dB)');

allElementPoint{1}(element_sum)=0;
allElementPoint{2}(element_sum)=0;
%��Ԫ�ֲ�ͼ
figure
plot(allElementPoint{1},allElementPoint{2},'bo');
xlabel('x(��)');
ylabel('y(��)');
%title('XLamda YLamda title');

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