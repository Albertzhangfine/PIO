% �Զ���Ŀ�꺯��fun
%PIO��Ⱥ�㷨
clc,clear,close all
warning off
format longG
T1=90;           % Global search algebra����������
T2=15;           % Local search algebra����������
pigeonnum=30;    % ��Ⱥ����
nvar = 1;             % δ֪������
R = 0.3;           % �شų�����parameters of magnetic field
bound=[-1,1];    % ������Χ
%% ��ʼ����Ⱥ
for i=1:pigeonnum
    pop(i,1) = bound(1) + (bound(2)-bound(1))*rand;
    fitness(i) = fun( pop(i,1)  );   % ��Ӧ�Ⱥ���
    v(i,1) = rand;                   % �����ٶ�
end
% ��¼һ������ֵ
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   %ȫ�����
gbest = pop;                % �������
fitnessgbest=fitness;     %���������Ӧ��ֵ
fitnesszbest=bestfitness; %ȫ�������Ӧ��ֵ
%% ��ͼ��ָ�������� magnetic compass and solar operator
for t=1:T1  % ��������
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
        fitness(i) = fun( pop(i,:)  );   % ��Ӧ�Ⱥ���
        % �Ƚ�  �����Ƚ�
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
%% �ر����� landmark operator
pop = gbest;                % �������
fitness = fitnessgbest;     %���������Ӧ��ֵ
for t=1:T2
    % sort the pigeons
    [a0,b0] = sort( fitness, 'ascend' );
    pop = pop(b0, :);
    fitness = a0;
    % ȡǰһ������Ž���з���
    pigeonnum1=ceil(pigeonnum/2);               % remove half of the pigeons according to the landmark
    % ���ӵ�����ֵ
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
        fitness(i) = fun( pop(i,:)  );   % ��Ӧ�Ⱥ���
        % �Ƚ�  �����Ƚ�
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
disp('���Ž�')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
