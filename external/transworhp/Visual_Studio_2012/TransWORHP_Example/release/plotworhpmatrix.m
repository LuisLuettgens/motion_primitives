function plotworhpmatrix(filename, style)

A = load( strcat(filename, '.m'),'-ASCII')

if (style==0)
    spy(A)
else
    bar3(A,'r.')
    axis([ 0 size(A,1) 0 size(A,2) ])
end

title(filename)

