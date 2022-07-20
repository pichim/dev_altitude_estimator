function [data]=read_bin_data(fn,nb_in)
%% Read Binary Data from SD-card
% typically there is a binary file .bin and a description file .txt of the
% same name. The .txt file describes columns like:
% gyr, 1:3
% acc, 4:6
% ...           please take care of correct syntax (comma right after
% desc., line break etc.)

% R. Altenburger 5.7.2019
if nargin<2
    nb_in=[];
end
if strcmp(fn(end-3),'.')        % if extension is given it is cut away
    fn=fn(1:end-4);
end
% ------------- open file
fid=fopen(strcat(fn,'.txt'),'r');       % open textfile 
if fid < 0          % old version with no text description
    vers = 0;
else
    vers = 1;
    c=fscanf(fid,'%c');
    fclose(fid);                
end
% - read binary data in one ...
fid=fopen(strcat(fn,'.bin'),'r');
if fid>0
    nb = fread(fid,1,'char'); % number data columns + time Vector at very beginning
    if ~isempty(nb_in)
        nb=nb_in;
    end
    a=fread(fid,inf,'float32');
    fclose(fid);
    N = floor(length(a)/nb)*nb;
    b=reshape(a(1:N),nb,N/nb)';
    data.ti=b(:,1);
    data.alldata=b(:,2:end);
else
    disp('cannot read file')
    data=[];
end
if vers ==1         % new version with file description
    data = convert_data(data,c);
end

% % hack pmic, 18.03.2021, maybe altb has a better idea for this
% if ( size(data.Lidar_z, 2) == 5 && size(data.Baro_alt, 2) == 3 ) && ( size(data.OF, 2) == 6 && size(data.RS, 2) == 11 )
%     time_offset = data.ti( find( diff(data.ti(1:10)) == max(diff(data.ti(1:10)))) + 1 );
%     data.ti              = data.ti - time_offset;
%     data.Lidar_z(:,end)  = data.Lidar_z(:,end) - time_offset;
%     data.Baro_alt(:,end) = data.Baro_alt(:,end) - time_offset;
%     data.OF(:,end)       = data.OF(:,end) - time_offset;
%     data.RS(:,end)       = data.RS(:,end) - time_offset;
% end

function data = convert_data(data,c)
while contains(c,',')
    k1=strfind(c,13);   % find line endings
    k1=k1(1)-1;         % select line
    li=c(1:k1);         % 
    c=c(k1+3:end);      % delete first line
    cm=strfind(li,','); % find comma
    eval(sprintf('data.%s = data.alldata(:,%s);',li(1:cm-1),li(cm+1:end)));
end