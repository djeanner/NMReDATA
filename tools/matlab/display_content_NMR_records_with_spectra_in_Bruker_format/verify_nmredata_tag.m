function [obj]=verify_nmredata_tag(tag_identyfyer,content_of_tag,header_tag)
verbose=1;

obj.tag_name=tag_identyfyer;
if verbose==2
    disp(['Tag:  <' tag_identyfyer '>'])
    %   disp(content_of_tag)
end

field_separator=',';
comment_separator=';';
%content_of_tag=strrep(content_of_tag,10,47)
%ascii code end-of-line:10
%ascii code of '\' :92
line_separator=[92];% inversion >1.0
list_of_line_break=strfind(content_of_tag,line_separator);

if size(list_of_line_break,2)==0 %check that there is at least one end of line...
    line_separator=[10];% inversion 1.0 only
    list_of_line_break=strfind(content_of_tag,line_separator);
    warning('no backslash ... using char 10 as seperator ... this is probably a verions 1 NMReDATA file')
else
    content_of_tag=strrep(content_of_tag,10,[]);
    list_of_line_break=strfind(content_of_tag,line_separator);
end


%list_of_line_break=findstr(content_of_tag,10);
from=1;
pieces_of_information=1;%will be incremented while reading the tag
if contains(tag_identyfyer,[header_tag '1D_'])
    obj.nb_dim=1;
    disp(['This a 1D spectrum' tag_identyfyer ' ' header_tag])
end
if contains(tag_identyfyer,[header_tag '2D_'])
    obj.nb_dim=2;
    disp(['This a 2D spectrum' tag_identyfyer ' ' header_tag])
end
for loop_over_content=list_of_line_break
    % extract each line from the tag paragraph
    current_line=content_of_tag(from:loop_over_content-size(line_separator,2));
    from=loop_over_content+size(line_separator,2);%set start position of next line
    comment='';
    % separate the comment form the rest
    list_of_line_comment_separator=strfind(current_line,comment_separator);
    if size(list_of_line_comment_separator,1)>0
        comment=current_line(list_of_line_comment_separator(1,1)+1:end);
        current_line=current_line(1:list_of_line_comment_separator(1,1)-1);%shorten the line
        if verbose==0
            disp(['comment : ' comment])
        end
        
    end
    if verbose==2
        disp(['<' current_line '>']);
    end
    if length(current_line)>1%there is something usefull in this line...
        %         if verbose==2
        %             disp(['Working on<' current_line '>']);
        %         end
        pos_equal=strfind(current_line,'=');
        found_variable=0;
        if size(pos_equal,1)>0
            if contains(current_line,'Larmor')
                got=sscanf(current_line(pos_equal(1,1)+1:end),'%f');
                if size(got,1)>0
                    obj.larmor=got(1,1);found_variable=1;
                end
            elseif contains(current_line,'Decoupled')
                obj.decoupled=current_line(pos_equal(1,1)+1:end);found_variable=1;
            elseif contains(current_line,'Interchangeable')
                obj.interchangeable=current_line(pos_equal(1,1)+1:end);found_variable=1;
            elseif contains(current_line,'CorrType')
                obj.cortype=current_line(pos_equal(1,1)+1:end);found_variable=1;
            elseif contains(current_line,'Pulseprogram')
                obj.pulseprogram=current_line(pos_equal(1,1)+1:end);found_variable=1;
            elseif contains(current_line,'Sequence')
                obj.sequence=current_line(pos_equal(1,1)+1:end);found_variable=1;
            elseif contains(current_line,'Spectrum_ID')
                obj.spectrumid=current_line(pos_equal(1,1)+1:end);found_variable=1;
            elseif contains(current_line,'Spectrum_Location')
                obj.spectrum_location=current_line(pos_equal(1,1)+1:end);found_variable=1;
            elseif contains(current_line,'zip_file_Location')
                obj.zip_file_Location=current_line(pos_equal(1,1)+1:end);found_variable=1;
            end
        end
        if ~found_variable%normal line of data
            obj.main_data{pieces_of_information}=current_line;
            obj.main_comment{pieces_of_information}=comment;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if contains(tag_identyfyer,[header_tag '1D'])
                obj=extract_data_1d(obj,current_line,pieces_of_information,field_separator);
            end
            if contains(tag_identyfyer,[header_tag '2D'])
                obj=extract_data_2d(obj,current_line,pieces_of_information,field_separator);
            end
            if contains(tag_identyfyer,[header_tag 'J'])
                warning('not reading J yet... write function first....')
             %   obj=extract_data_J(obj,current_line,pieces_of_information,field_separator);

            end
            if contains(tag_identyfyer,[header_tag 'ASSIGNMENT'])
                obj=extract_data_signals(obj,current_line,pieces_of_information,field_separator);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pieces_of_information=pieces_of_information+1;
            
        else
            if isfield(obj,'comment')
                obj.comment=[obj.comment  comment 10];
            else
                obj.comment=[  comment 10];
            end
        end
    else
        if isfield(obj,'comment')
            obj.comment=[obj.comment  comment 10];
        else
            obj.comment=[  comment 10];
        end
    end
end

if verbose==2
    obj
end

end

