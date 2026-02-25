clear
clc

%***********************Building Matrix M ****************************

input_path1 = 'scalefree_set1.txt';
input_path2 = 'scalefree_set2.txt';

% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = 50;
N = 2*n1;
M = zeros(N,N);

j = n1+1;
for i=1:n1
   M(i,j) = 1; 
   M(j,i) = 1;
   j = j + 1; 
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%---LAYER-1----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % %------------------------SF--------------------------
[A1,n1] = read_large_sparse_network(input_path1);
A1 = full(A1);

% %-------------------------ER------------------------
% k1 = 3; p1 = k1/n1;
% [A1,~] = ER_Random(n1, p1); 

% %::::::::::::::: Copy A1 to M ::::::::::::::::::::::::::
M(1:n1,1:n1) = A1;
clear A1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---LAYER-2----%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %------------------------SF--------------------------
[A2,n1] = read_large_sparse_network(input_path2);
A2 = full(A2);

% %-------------------------ER------------------------
% k1 = 3; p1 = k1/n1;
% [A2,~] = ER_Random(n1, p1); 


% %::::::::::::::: Copy A2 to M ::::::::::::::::::::::::::
M(n1+1:N,n1+1:N) = A2;
clear A2

%***********************End of Building Matrix M **************************

V = size(M,1);
d = sum(M,2);
f = cell(V,1);

for r=1:V
f{r}=find(M(r,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cd1= cell(7,1);
Cd2= cell(7,1);

%Broken Links
L1=[1,10,30,40,50,60,83]; 
L2=[f{1}(1),f{10}(1),f{30}(1),f{40}(1),f{50}(1),f{60}(1),f{83}(1)];

% L1=[1,10,30,40,50]; 
% L2=[f{1}(1),f{10}(1),f{30}(1),f{40}(1),f{50}(1)];

for h=1:length(L1)

    Cd_temp1=zeros(d(L1(h))-1,d(L1(h))-1);
    Cd_temp2=zeros(d(L2(h))-1,d(L2(h))-1);

    for r=1:d(L1(h))-1
        for s=1:d(L1(h))-1
            Cd_temp1(r,s)=(1/sqrt(d(L1(h))-1))*exp((2*pi*sqrt(-1)*(r-1)*(s-1))/(d(L1(h))-1));
        end
    end

     for r=1:d(L2(h))-1
        for s=1:d(L2(h))-1
            Cd_temp2(r,s)=(1/sqrt(d(L2(h))-1))*exp((2*pi*sqrt(-1)*(r-1)*(s-1))/(d(L2(h))-1));
        end
    end
    
    Cd_temp1(end+1, :) = 0; 
    Cd_temp1([find(f{L1(h)}==L2(h)) end], :) = Cd_temp1([end find(f{L1(h)}==L2(h))], :);
    Cd_temp1(:, end+1) = 0;
    Cd_temp1(:, [find(f{L1(h)}==L2(h)) end]) = Cd_temp1(:, [end find(f{L1(h)}==L2(h))]);
    Cd_temp1(find(f{L1(h)}==L2(h)),find(f{L1(h)}==L2(h)))=1;

    Cd1{h}=Cd_temp1;

    
    Cd_temp2(end+1, :) = 0; 
    Cd_temp2([find(f{L2(h)}==L1(h)) end], :) = Cd_temp2([end find(f{L2(h)}==L1(h))], :); 
    Cd_temp2(:, end+1) = 0;
    Cd_temp2(:, [find(f{L2(h)}==L1(h)) end]) = Cd_temp2(:, [end find(f{L2(h)}==L1(h))]);
    Cd_temp2(find(f{L2(h)}==L1(h)),find(f{L2(h)}==L1(h)))=1;

    Cd2{h}=Cd_temp2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

W=zeros(1,V);

trial=1000;

for u=1:trial

% Probability amplitute matrix
  A = zeros(V,V);
  Atemp = zeros(V,V);

%Initial state1
A(1,f{1}(1))=1;
  
step = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fouriour Coin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:step

    qq=rand;

    for x=1:V % Select the vertex

            if qq<=0.5 && x==L1(1)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd1{1}(r,s);
                end
            end

            elseif qq<=0.5 && x==L1(2)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd1{2}(r,s);
                end
            end

            elseif qq<=0.5 && x==L1(3)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd1{3}(r,s);
                end
            end

            elseif qq<=0.5 && x==L1(4)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd1{4}(r,s);
                end
            end

            elseif qq<=0.5 && x==L1(5)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd1{5}(r,s);
                end
            end

            elseif qq<=0.5 && x==L1(6)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd1{6}(r,s);
                end
            end
            elseif qq<=0.5 && x==L1(7)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd1{7}(r,s);
                end
            end

            elseif qq<=0.5 && x==L2(1)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd2{1}(r,s);
                end
            end

            elseif qq<=0.5 && x==L2(2)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd2{2}(r,s);
                end
            end

            elseif qq<=0.5 && x==L2(3)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd2{3}(r,s);
                end
            end

            elseif qq<=0.5 && x==L2(4)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd2{4}(r,s);
                end
            end

            elseif qq<=0.5 && x==L2(5)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd2{5}(r,s);
                end
            end

            elseif qq<=0.5 && x==L2(6)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd2{6}(r,s);
                end
            end

            elseif qq<=0.5 && x==L2(7)
            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*Cd2{7}(r,s);
                end
            end

            else

            for r=1:d(x)
                for s=1:d(x)
                    Atemp(x,(f{x}(r))) = Atemp(x,(f{x}(r)))+ A(x,(f{x}(s)))*(1/sqrt(d(x)))*exp((2*pi*sqrt(-1)*(r-1)*(s-1))/d(x));
                end
            end

            end
        
    end  

    A = transpose(Atemp);
    Atemp = zeros(V,V);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Total Proability calculation
VertexProb = zeros(1,V);

for k=1:V
   VertexProb(1,k) = A(k,:)*conj(transpose(A(k,:)));
end

W=W+VertexProb;

end

for k=1:V
   x(k)=k;
end


bar(x,W/trial)

sum(W/trial)


% %xlim([0,1])
% ylim([0,1])
lablel_x=xlabel({'$$ Nodes $$','$$ (d)$$'},'Interpreter','latex');
label_y =ylabel('$$Probability$$','Interpreter','latex');
label_y.Rotation		= 90;
% 
% 
set(gca,'FontSize',40)
set(gca,'Linewidth',8)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
fig = gcf;
exportgraphics(fig,'Fouriour_SF_SF_decoherence_7link.png','Resolution',300)
% exportgraphics(fig,'Heatmap_Fouriour_Coin_second_initial_state.png','Resolution',1000)






