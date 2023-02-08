function SaveMatrix(M,filename);
% Purpose: To save a matrix in a text file that is Fortran readable where the first line is the size of the matrix
% Inputs:   M           -  Matrix to be saved
%           filename    -  String of name of the file (without file extension)    

% Author: Sam MacPherson

%writematrix(size(M),strcat(filename,".dat"),"delimiter","tab");
% writematrix(M,strcat(filename,".dat"),"delimiter","tab","WriteMode","append",'-double');

fid = fopen(strcat(filename,".dat"),'wt');

[a,b]= size(M);
fprintf(fid,'%i',a);
fprintf(fid,'%c',' ');
fprintf(fid,'%i',b);
fprintf(fid,'\n');
for ii = 1:size(M,1)
    fprintf(fid,'%58.50e',M(ii,:));
    fprintf(fid,'\n');
end


fclose(fid);
return
    