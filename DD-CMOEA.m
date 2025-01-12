classdef DD-CMOEA < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>


    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population1 = Problem.Initialization(); % Main task
%             Fitness1    = CalFitness(Population1.objs,Population1.cons);
            Population2 = Problem.Initialization(); % Global auxiliary task
            Population = [Population1,Population2];
    
           
           %% Generate the weight vectors
            [RW,WVs] = UniformPoint(Problem.N,Problem.M);
            T = ceil(WVs/10);
            nr = ceil(WVs/100);
            % 创建一个新矩阵，用于存储反转后的行
            W = zeros(size(RW));
            % 逐行反转
            for i = 1:WVs
                W(i, :) = RW(WVs - i + 1, :);
            end
            
            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);
 
            
           %% Angle
            angle=acos(1-pdist2(W,W,'cosine'));
            temp_angle=angle;
            temp_angle(logical(eye(size(temp_angle))))=inf;
            theta_min=min(temp_angle');
            theta_min=theta_min';
            theta=theta_min.*0.5;
            
            
            arch = ArchiveUpdate(Population,Problem.N);
            Population = arch;
            gen  = 1;	% index of generation
            flag = 0;
            VAR0=0;

            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                 
                Z = min(Population2.objs,[],1);
                if flag == 0
                    std_obj(gen,:) = std(Population2.objs,[],1);
                    if gen>100
                        if  sum(std(std_obj(gen-100:gen,:),[],1)<0.1) == Problem.M
                            flag = 1;
                            Fes = Problem.FE;
                            cons = Population1.cons;
                            cons(cons<0) = 0;
                            cons =sum(cons,2);%（矩阵行和）
                            index =find(cons>0);
                            if isempty(index)
                                VAR0 = 0;
                            else
                                VAR0 =  mean(cons(index));%O1中所有不可行解均值
                            end
                        end
                    end
                end
           
                %%
                for i = 1:Problem.N
                    ii = i;
                    if i > WVs
                        ii = randi([1,WVs]);
                    end                   
                    
                    if rand < 0.2
                        PP{i} = randperm(WVs);
                        P = PP{i};
                        E(i,:)=P(1:2);
                    else
                        PP{i} = B(ii,randperm(size(B,2)));%在矩阵 B 中选择第 i 行的元素，并对这些元素进行随机排列
                        P = PP{i};
                        E(i,:)=P(1:2);
                    end
                end

                Offspring1 = OperatorGA(Problem,Population1(randi(Problem.N,1,Problem.N)));  
                if flag == 0
                    Offspring2 = OperatorGAhalf(Problem,[Population2(E(:,1)),Population2(E(:,2))]); 
                else
                    Offspring2 = OperatorDE(Problem,Population2,Population2(E(:,1)),Population2(E(:,2)));
                end
                
               %%                
               
               Offspring = Offspring2;
               % For each solution
               for i = 1 : Problem.N
                   P = PP{i};
                   Z = min(Z,Offspring(i).obj);
                   PopObj2=Population2(P).objs-repmat(Z,length(P),1);
                   Angle0 = acos(1 - pdist2(real(PopObj2),W(P,:),'cosine'));
                   Angle2 = diag(Angle0);
                   NewAngle2=acos(1-pdist2(real(Offspring(i).objs-Z),W(P,:),'cosine'));
                   NewAngle2=NewAngle2';
                   
                   g_old2 = max(abs(Population2(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                   g_new2 = max(repmat(abs(Offspring(i).obj-Z),length(P),1).*W(P,:),[],2);
                   
                   if (flag == 0)
                       case1=NewAngle2<= theta(P) & Angle2<= theta(P) & g_old2>=g_new2;%both in niche
                       case2=NewAngle2> theta(P) & Angle2> theta(P) & g_old2>=g_new2;%both not in niche
                       case3=NewAngle2<= theta(P) & Angle2> theta(P)& g_old2>=g_new2;
                       case4=NewAngle2> theta(P) & Angle2<= theta(P) & g_old2>=g_new2;
                       Population2(P(find(case1 | case2 | case3 | case4 ,nr))) = Offspring(i);
                   else
                       CVO2 = sum(max(0,Offspring(i).con));
                       CVP2 = sum(max(0,Population2(P).cons),2);
                       case1=NewAngle2<= theta(P) & Angle2<= theta(P) & (CVP2>CVO2 | (CVP2==CVO2 & g_old2>=g_new2));%both in niche
                       case2=NewAngle2<= theta(P) & Angle2> theta(P) & (CVP2>CVO2 | (CVP2==CVO2 & g_old2>=g_new2));
                       case3=NewAngle2<= theta(P) & Angle2> theta(P);
                       indices = mod(find([case1 ; case2 ; case3],nr),length(P));
                       if indices > 0
                           Population2(P(indices)) = Offspring(i);
                       end
                   end
               end
               
               
               if flag == 0
                   [Population1,~] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N,true);
               else
                   % Adaptive constraint search boundary
                   [Population1,~] = EnvironmentalSelection_VAR([Population1,Offspring1,Offspring2,arch],Problem.N,VAR0);
                   cons = Offspring1.cons;
                   cons(cons<0) = 0;
                   cons = sum(cons,2);
                   index = find(cons>0);
                   if isempty(index)
                       VAR0 = 0;
                   else
                       tanhFE = 8 * ((Problem.FE-Fes) / (Problem.maxFE-Fes)) - 4;
                       K = 0.5-tanh(tanhFE)/2;
                       VAR0 = K*mean(cons(index));
                   end
               end
               
               %%
               
               Population = [Population1,Population2];
               % Output the non-dominated and feasible solutions.
               arch = ArchiveUpdate([arch,Population],Problem.N);
               
               gen =gen +1;
               
               if Problem.FE >=Problem.maxFE
                   Population = arch;
               end
               
            end
        end
    end
end
