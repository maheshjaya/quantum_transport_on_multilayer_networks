clear;
clc;
clear all;
close all;

i=sqrt(-1);
%Multilayer with 4 nodes and 2 layers
%In this test we are calclauting the probability of finding the quantum walker at each vertex vs time
%The probability of finding the quantum walker on layer L vs time
M =[1 1 1 1 1 0 0 0; 
    1 1 1 1 0 1 0 0; 
    1 1 1 1 0 0 1 0; 
    1 1 1 1 0 0 0 1;
    1 0 0 0 0 1 0 1;
    0 1 0 0 1 0 1 0;
    0 0 1 0 0 1 0 1;
    0 0 0 1 1 0 1 0]; 

V = size(M,1);
d = sum(M,2);

f = cell(V,1);

 for r=1:V
 f{r}=find(M(r,:));
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step_list=0:100;

T=5:5:length(step_list); %Time gaps for time average probability

Temp_array=zeros(1,V);
Time_avarage_probability=zeros(V,length(T));

ww=1;

for qq=1:length(step_list)

step=step_list(qq);

% Probability amplitute matrix
A=zeros(V,V);
Atemp=zeros(V,V);


%initial state1
% A(1,2) = 1/sqrt(4);
% A(1,3) = 1/sqrt(4);
% A(1,4) = 1/sqrt(4);
% A(1,5) = 1/sqrt(4);

A(1,2) = 1/sqrt(2);
A(5,6) = 1/sqrt(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%............................Fourier Coin.........................
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:step

for x=1:V

        for r=1:d(x)
            for s=1:d(x)
                Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*(1/sqrt(d(x)))*exp((2*pi*i*(r-1)*(s-1))/d(x)); 
            end
        end
    
end

A = transpose(Atemp);
Atemp = zeros(V,V);

end

%Total Proability calculation
VertexProb = zeros(1,V);

for k=1:V
   VertexProb(k) = A(k,:)*conj(transpose(A(k,:)));
end


%Time average Proability calculation 
Temp_array=Temp_array+VertexProb;

if step==T(ww)
    for aa=1:V
    Time_avarage_probability(aa,ww)=(1/(T(ww)+1))*Temp_array(1,aa);
    end
    ww=ww+1;
end



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time average of Layer 1 and Layer 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time_average_Layer1=zeros(1,length(T));
Time_average_Layer2=zeros(1,length(T));

for pp=1:V/2
Time_average_Layer1=Time_average_Layer1+Time_avarage_probability(pp,:);
end

for pp=(V/2+1):V
Time_average_Layer2=Time_average_Layer2+Time_avarage_probability(pp,:);
end


figure(1)
plot(T,Time_average_Layer1,'-*','LineWidth',2,'DisplayName','Layer 1');
hold on
plot(T,Time_average_Layer2,'-.','LineWidth',2,'DisplayName','Layer 2');
title('$$Time \ average \ probability \ of \ finding \ the \ quantum \ walker \ at \ different \ layers $$','Interpreter','latex')
xlim([5,100])
%ylim([0,1])
lablel_x=xlabel('$$T$$','Interpreter','latex');
label_y =ylabel('$$Probability$$','Interpreter','latex');
label_y.Rotation		= 90;

lgd = legend;
lgd.NumColumns = 2;


