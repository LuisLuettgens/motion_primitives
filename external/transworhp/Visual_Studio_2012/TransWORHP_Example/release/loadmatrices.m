%% WORHP Matrizen

subplot(3,1,1)
plotworhpmatrix('DF_0',1);

subplot(3,1,2)
plotworhpmatrix('DG_0',1);

subplot(3,1,3)
plotworhpmatrix('HM_0',1);

%% WORHP Zen Matrizen

subplot(3,2,1)
plotworhpmatrix('ZenDX_0',0);

subplot(3,2,2)
plotworhpmatrix('ZenDMu_0',0);

subplot(3,2,3)
plotworhpmatrix('ZenDF_0',0);

subplot(3,2,4)
plotworhpmatrix('ZenDF2_0',0);

subplot(3,2,5)
plotworhpmatrix('ZenDG_0',0);


%% Referenzwert (0) und 3 St?rungen laden, um es mit SensAbl. zu vergleichen
clear

X0 = load('X_0.m', '-ASCII');
X1 = load('X_1.m', '-ASCII');
X2 = load('X_2.m', '-ASCII');
X3 = load('X_3.m', '-ASCII');
ZenDX0 = load('ZenDX_0.m', '-ASCII');

r = 4;
step = 0.01
%%

for i=1:4
    subplot(r,4,i)

    plot( X0(i:4:end) , 'k*')
    hold on
    plot( X1(i:4:end) , 'r' )
    plot( X2(i:4:end) , 'g' )
    plot( X3(i:4:end) , 'b' )
    if i==1
        legend('Ref','pert 1','pert 2','pert 3')
    end
    if i==1
        title('State X1')
    end
    if i==2
        title('State X2')
    end
    if i==3
        title('State X3')
    end
    if i==4
        title('Control U1')
    end
end


% Sens 1
for i=1:4
    subplot(r,4,4+i)
    plot( (X1(i:4:end)- X0(i:4:end)) / step )
    hold on
    plot( ZenDX0(i:4:end,1), 'r')
    %plot( ZenDX0(i:4:end,48), 'g')
    plot( ZenDX0(i:4:end,78), 'g')
    if i==1
        legend('Diffquot','Zen p', 'Zen q')
        title('Sens 1: in Anfwert X1')
    end
end
% Sens 2
for i=1:4
    subplot(r,4,2*4+i)
    plot( (X2(i:4:end)- X0(i:4:end)) / step )
    hold on
    plot( ZenDX0(i:4:end,2), 'r')
    plot( ZenDX0(i:4:end,79), 'g')
 
    if i==1
        legend('Diffquot','Zen p','Zen q')
        title('Sens 2: in Endwert X2')
    end
end
% Sens 3
for i=1:4
    subplot(r,4,3*4+i)
    plot( (X3(i:4:end)- X0(i:4:end)) / step )
    hold on
    plot( ZenDX0(i:4:end,3), 'r')
 
    if i==1
        legend('Diffquot','Zen')
        title('Sens 3: in Dynamik')
    end
end
