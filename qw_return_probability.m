clear;
clc;
i=sqrt(-1);


M =[0 1 1 1 1 0 0 0; 
    1 0 1 1 0 1 0 0; 
    1 1 0 1 0 0 1 0; 
    1 1 1 0 0 0 0 1;
    1 0 0 0 0 1 0 1;
    0 1 0 0 1 0 1 0;
    0 0 1 0 0 1 0 1;
    0 0 0 1 1 0 1 0]; 

step=100;

T=1:5:step; %Polya number parameter

initial_node=4; % We initialize the walker from this node 

V = size(M,1);
d = sum(M,2);

f = cell(V,1);

 for r=1:V
 f{r}=find(M(r,:));
 end
 
% Probability amplitute matrix
A=zeros(V,V);
Atemp=zeros(V,V);


%initial state
for h=1:d(initial_node)
A(initial_node,f{initial_node}(h)) = 1/sqrt(d(initial_node));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grover Coin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P=ones(1,length(T)); %Polya number calculation
X=1;
Y=1;

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

        %Polya number calculation
        X=X*(1-A(initial_node,:)*conj(transpose(A(initial_node,:))));
        if sum(T==t)
            P(1,Y) = 1-X;
            Y=Y+1;
        end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fourier Coin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P=ones(1,length(T)); %Polya number calculation
% X=1;
% Y=1;
% 
% for t=1:step
%    for x=1:V % Select the vertex
%     for r=1:d(x)
%         for s=1:d(x)
%             Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*(1/sqrt(d(x)))*exp((2*pi*i*(r-1)*(s-1))/d(x));
%         end
%     end
%    end    
%  %end
% A = transpose(Atemp);
% Atemp = zeros(V,V);
% 
%         %Polya number calculation
%         X=X*(1-A(5,:)*conj(transpose(A(5,:))));
%         if sum(T==t)
%             P(1,Y) = 1-X;
%             Y=Y+1;
%         end
% 
% end


plot(T,P)




