% 自定义目标函数fun
%PIO鸽群算法
clc,clear,close all
warning off
format longG
T1=90;           % Global search algebra，迭代次数
T2=15;           % Local search algebra，迭代次数
pigeonnum=30;    % 种群数量
nvar = 1;             % 未知量个数
R = 0.3;           % 地磁场参数parameters of magnetic field
bound=[-1,1];    % 搜索范围
%% 初始化种群
for i=1:pigeonnum
    pop(i,1) = bound(1) + (bound(2)-bound(1))*rand;
    fitness(i) = fun( pop(i,1)  );   % 适应度函数
    v(i,1) = rand;                   % 飞行速度
end
% 记录一组最优值
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   %全局最佳
gbest = pop;                % 个体最佳
fitnessgbest=fitness;     %个体最佳适应度值
fitnesszbest=bestfitness; %全局最佳适应度值
%% 地图和指南针算子 magnetic compass and solar operator
for t=1:T1  % 迭代次数
    for i=1:pigeonnum
        v(i,:)=v(i,:)*(1-exp(-R*t))+rand*(gbest(i,:)-pop(i,:));
        pop(i,:)=pop(i,:)+v(i,:);   %check whether beyond the searching space
        for j=1:nvar
            if abs(i-1)<=eps
                if pop(i,j)<bound(1)||pop(i,j)>bound(2)
                    pop(i,j)=bound(1)+rand*(bound(2)-bound(1));
                    pop(i,j)=rand;
                end
            else
                if pop(i,j)<bound(1)||pop(i,j)>bound(2)
                    pop(i,j)=pop(i-1,j);
                    v(i,j)=v(i-1,j);
                end   
            end
        end
        fitness(i) = fun( pop(i,:)  );   % 适应度函数
        % 比较  个体间比较
        if fitness(i)<fitnessgbest(i)
            fitnessgbest(i) = fitness(i);
            gbest(i,:) = pop(i,:);
        end
        if fitness(i)<bestfitness
            bestfitness = fitness(i);
            zbest =  pop(i,:);
        end
    end
    fitness_iter(t) = bestfitness;
end
%% 地标算子 landmark operator
pop = gbest;                % 个体最佳
fitness = fitnessgbest;     %个体最佳适应度值
for t=1:T2
    % sort the pigeons
    [a0,b0] = sort( fitness, 'ascend' );
    pop = pop(b0, :);
    fitness = a0;
    % 取前一半的最优解进行分析
    pigeonnum1=ceil(pigeonnum/2);               % remove half of the pigeons according to the landmark
    % 鸽子的中心值
    addpigeonnum = sum( pop(1:pigeonnum1, :) )  ;                    
    pigeoncenter=ceil(addpigeonnum./pigeonnum); % calculate central position
    for i=1:pigeonnum
        v(i,:)=v(i,:)*(1-exp(-R*t))+rand*(gbest(i,:)-pop(i,:));
        pop(i,:)=pop(i,:)+v(i,:);   %check whether beyond the searching space
        for j=1:nvar
            if abs(i-1)<=eps
                if pop(i,j)<bound(1)||pop(i,j)>bound(2)
                    pop(i,j)=bound(1)+rand*(bound(2)-bound(1));
                    pop(i,j)=rand;
                end
            else
                if pop(i,j)<bound(1)||pop(i,j)>bound(2)
                    pop(i,j)=pop(i-1,j);
                    v(i,j)=v(i-1,j);
                end   
            end
        end
        fitness(i) = fun( pop(i,:)  );   % 适应度函数
        % 比较  个体间比较
        if fitness(i)<fitnessgbest(i)
            fitnessgbest(i) = fitness(i);
            gbest(i,:) = pop(i,:);
        end
        if fitness(i)<bestfitness
            bestfitness = fitness(i);
            zbest =  pop(i,:);
        end
    end
    fitness_iter(T1+t) = bestfitness;
end
disp('最优解')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
