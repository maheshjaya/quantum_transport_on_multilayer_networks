clear;
clc;
i=sqrt(-1);

%Adjacency Matrix

M =[0 1 1 1 1 0 0 0; 
    1 0 1 1 0 1 0 0; 
    1 1 0 1 0 0 1 0; 
    1 1 1 0 0 0 0 1;
    1 0 0 0 0 1 0 1;
    0 1 0 0 1 0 1 0;
    0 0 1 0 0 1 0 1;
    0 0 0 1 1 0 1 0]; 

%step size
step=10000;


V = size(M,1);
d = sum(M,2);

%function f_x(r)
f = cell(V,1);
for r=1:V
 f{r}=find(M(r,:));
end
 

% Probability amplitute matrix
A=zeros(V,V);
Atemp=zeros(V,V);


%initial state
A(1,2) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fourier Coin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:step
   for x=1:V % Select the vertex
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
P = 0;
VertexProb = zeros(V,1);

for k=1:V
   P = P + A(k,:)*transpose(conj(A(k,:)));
   VertexProb(k) = A(k,:)*conj(transpose(A(k,:)));
   x(k)=k;
end

% Total_probability 
fprintf('Total_probability: %f\n',P);
%bar(x,VertexProb,0.5)
bar(x,VertexProb)
xlabel("Nodes")
ylabel("Probability")


