clear

load boylat.mat

% save file
fileID = fopen('buoies.txt','w');

for i=1:length(boylat)
    
    fprintf(fileID,'%12.8f ',boylat{i,1});
    fprintf(fileID,'%12.8f ',boylat{i,2});
    fprintf(fileID,'%i ',boylat{i,3});
    fprintf(fileID,'%s\n',boylat{i,4});
    
end
fclose(fileID);

