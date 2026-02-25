clear;
clc;
clear all;
close all;

i=sqrt(-1);

%Adjacency Matrix of the multilayer network
M =[0 1 1 1 0 0 0 0; 
    1 0 1 0 0 0 0 0; 
    1 1 0 0 0 0 0 0; 
    1 0 0 0 0 0 0 1;
    0 0 0 0 0 1 1 0;
    0 0 0 0 1 0 1 0;
    0 0 0 0 1 1 0 1;
    0 0 0 1 0 0 1 0]; 

step=100;

V = size(M,1);
d = sum(M,2);

f = cell(V,1);

 for r=1:V
 f{r}=find(M(r,:));
 end
 
% Probability amplitute matrix
A=zeros(V,V);
Atemp=zeros(V,V);


theta=45;

C=zeros(2,2);

C(1,1)=cosd(theta);
C(1,2)=sind(theta);
C(2,1)=sind(theta);
C(2,2)=-cosd(theta);

% phi=38;
% 
% C=zeros(2,2);
% C(1,1)=1/sqrt(2);
% C(1,2)=(1/sqrt(2))*exp(i*phi*(pi/180));
% C(2,1)=(1/sqrt(2))*exp(-i*phi*(pi/180));
% C(2,2)=-1/sqrt(2);

%initial state1
A(1,2) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grover Coin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:step

for x=1:V
    for r=1:d(x)
        for s=1:d(x)
            if (r==s) 
                Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*(1/d(x))*(2-d(x));
            else
                Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*(2/d(x));
            end
        end
    end
end

A = transpose(Atemp);
Atemp = zeros(V,V);

end


%Total Proability calculation
P = 0;
VertexProb = zeros(V,1);

for k=1:V
   P = P + A(k,:)*transpose(conj(A(k,:)));
   VertexProb(k) = A(k,:)*conj(transpose(A(k,:)));
   x(k)=k;
end



% Total_probability 
disp('Probabilities related to QW ');
fprintf('Total_probability_QW: %f\n',P);
fprintf('Probability in first layer: %f\n',VertexProb(1)+VertexProb(2)+VertexProb(3)+VertexProb(4));
fprintf('Probability in second layer: %f\n',VertexProb(5)+VertexProb(6)+VertexProb(7)+VertexProb(8));


subplot(1,2,1);
%bar(x,VertexProb,0.5)
bar(x,VertexProb)
title("QW")
xlabel("Nodes")

hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Classical Walk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Probability amplitute matrix
B=zeros(V,1);
Btemp=zeros(V,1);


%initial state1
B(1) = 1;

for t=1:step
   for x=1:V % Select the vertex
        for r=1:d(x)
            Btemp(x) = Btemp(x)+ B(f{x}(r))*(1/d(f{x}(r)));
        end
   end    

B = Btemp;
Btemp = zeros(V,1);
end

%Total Proability calculation
PCal = 0;
VertexProb_Cal = zeros(V,1);

for k=1:V
   PCal = PCal + B(k);
   x(k)=k;
end

% Total_probability 
disp('Probabilities related to CRW ');
fprintf('Total_probability_CRW: %f\n',PCal);
fprintf('Probability in first layer: %f\n',B(1)+B(2)+B(3)+B(4));
fprintf('Probability in second layer: %f\n',B(5)+B(6)+B(7)+B(8));

subplot(1,2,2);
%bar(x,VertexProb,0.5)
bar(x,B)
title("CRW")
xlabel("Nodes")



