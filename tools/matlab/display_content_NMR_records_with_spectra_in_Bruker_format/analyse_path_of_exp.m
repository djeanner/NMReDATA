function        [path_exp, acqno, procno]=analyse_path_of_exp(input_path)

if nargin ==0
    input_path='/AN-menthol/12/pdata/1/';
    % input_path='toto/AN-menthol/12/pdata/1';
    % input_path='toto/AN-menthol/1';
    %input_path='AN-menthol/12';
    %input_path='./AN-menthol/12';
end
%% in case still there remove "file=" of "file:"
tmp=strfind(input_path,'file');
if size(tmp,2)>0
    p=tmp(1,1);
    if tmp(1,1)>0
        %  disp(input_path)
        input_path=input_path(1,p+5:end);
        %         disp(input_path)
        
    end
end

%% in case remove "backslash"
tmp=strfind(input_path,'\');
if size(tmp,2)>0
    p=tmp(1,1);
    if tmp(1,1)>0
        %   disp(input_path)
        input_path=input_path(1,1:p-1);
        %   disp(input_path)
        
    end
end

%% if no "/pdata/1/"... will add it
tmp=strfind(input_path,'pdata');
if size(tmp,2)==0
    %disp(input_path)
    if input_path(1,end)=='/'
        input_path=[input_path 'pdata/1/'];
    else
        input_path=[input_path '/pdata/1/'];
        
    end
    %      disp(input_path)
end

%% if no "/pdata/1/"... will add it
tmp=strfind(input_path,'pdata');
if size(tmp,2)==0
    %disp(input_path)
    if input_path(1,end)=='/'
        input_path=[input_path 'pdata/1/'];
    else
        input_path=[input_path '/pdata/1/'];
        
    end
    %      disp(input_path)
    
end

%% add "/" if not finishing with it...
if input_path(1,end)~='/'
    input_path=[input_path '/'];
end

%disp(input_path)
tmp=strfind(input_path,'/');
si=size(tmp,2);
if si<4
    error('problem with number of elements in path')
else
    %     if si==4
    %     tmp=[0 tmp(1,si-3:end)];
    %     else
    tmp=[ tmp(1,si-3:end)];
    %     end
    % tmp
end
path_exp=input_path(1:tmp(1,1));
% disp(path_exp)
acqno=str2num(input_path(1,tmp(1,1)+1:tmp(1,2)-1));
procno=str2num(input_path(1,tmp(1,3)+1:tmp(1,4)-1));
% disp(num2str(acqno))
% disp(num2str(procno))
end
