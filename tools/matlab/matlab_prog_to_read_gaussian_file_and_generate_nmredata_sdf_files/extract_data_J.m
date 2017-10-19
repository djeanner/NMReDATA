function obj=extract_data_J(obj,current_line,pieces_of_information,field_separator)

%example of current_line : (2), H1(C2), 149.61;1J(C,H)
list_separator=strfind(current_line,field_separator);

if size(list_separator,2)<2
    error_message=['J expect three fiels : '   'inside'  current_line];
    disp(['Error :  ' error_message]);
else
    first_label_or_chemshift =current_line(1:list_separator(1,1)-1);
    second_label_or_chemshift=current_line(list_separator(1,1)+1:list_separator(1,2)-+1);
    coupling                 =current_line(list_separator(1,2)+1:end);
    %    disp(['correlation between ' first_label_or_chemshift ' and ' second_label_or_chemshift])
    obj.signalf1{pieces_of_information}=first_label_or_chemshift;
    obj.signalf2{pieces_of_information}=second_label_or_chemshift;
    value=sscanf(coupling,'%f');
    obj.J{pieces_of_information}=value;
end


end