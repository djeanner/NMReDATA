function obj=extract_data_J(obj,current_line,pieces_of_information,field_separator,verbose)
if nargin<5
    verbose=0;
end
%example of current_line : H1(C2), 8.2433, H1
list_separator=strfind(current_line,field_separator);

if size(list_separator,2)<2%for 3
    error_message=['signals lines expect three fields : '   'inside'  current_line];
    disp(['Error1 :  ' error_message]);
else
    jp1 =current_line(1:list_separator(1,1)-1);
    jp2=current_line(list_separator(1,1)+1:list_separator(1,2)-1);
    
    j_coupling=current_line(list_separator(1,2)+1:end);
    
    obj.label_signal1{pieces_of_information}=jp1;
    obj.label_signal2{pieces_of_information}=jp2;
    value=sscanf(j_coupling,'%f');
    obj.J{pieces_of_information}=value;
    if verbose
        disp(['J(' jp1 ',' jp2 ')=' num2str(value)  ' Hz'])
    end
end

end

%td-alliance.github.io/metat^ata-directory
