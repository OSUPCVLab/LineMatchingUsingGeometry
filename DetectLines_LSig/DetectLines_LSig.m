function [data] = DetectLines_LSig(inname,outname,minLength)

%     fpath = fileparts(which(mfilename));
%     exec_str = ['"' fpath '\LineDetect.exe"'];
% 
%     system([exec_str  '  ' inname '  ' outname ]);
%     
% %     fid = fopen(outname);
% %     if fid==-1
% %     error('Cannot load results from binary executable');
% %     end
%     %file_header = textscan(fid, '%s', 1, 'delimiter', '|');
% 
%     % read numeric data
%     data = dlmread(outname, '%f %f %f %f');
%     %data=([data{[1:length(data)]}])'; 
%     data=data';
%     ind=(sqrt(sum((data(1:2,:)-data(3:4,:)).^2,1)))>=minLength;
%     data=data(:,ind);
% %     fclose(fid);

%%
    fpath = fileparts(which(mfilename));
    exec_str = ['"' fpath '\LineDetect.exe"'];

    system([exec_str  '  ' inname '  ' outname ]);
    
    fid = fopen(outname);
    if fid==-1
    error('Cannot load results from binary executable');
    end
    file_header = textscan(fid, '%s', 1, 'delimiter', '|');

    % read numeric data
    data = textscan(fid, '%f %f %f %f');
    data=([data{[1:length(data)]}])'; 
    ind=((sqrt(sum((data(1:2,:)-data(3:4,:)).^2,1)))>=minLength);
    data=data(:,ind);
    fclose(fid);

end