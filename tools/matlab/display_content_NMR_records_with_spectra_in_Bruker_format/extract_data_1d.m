function obj=extract_data_1d(obj,current_line,pieces_of_information,field_separator)

%example of current_line : 9.06, L=H31(C25), S=d, J=8.44(H28(C27)), E=1
%extract first field, the chemical shift can be a float or a two float
%separated by "-" (eg. 7.5-7.2)
verbose=0;
list_separator=strfind(current_line,field_separator);
if size(list_separator,1)>0
    chem_shift_string=current_line(1:list_separator(1,1)-1);current_line=current_line(list_separator(1,1)+1:end);
else
    chem_shift_string=current_line;current_line='';%no other field
end
chem_shifts=sscanf(chem_shift_string,'%f-%f');%this returns either one or two values.
if size(chem_shifts,1)<1
    error_message=['1 expect chemical shift as first field of line : '   current_line ' ' chem_shift_string];
    disp(['Error :  ' error_message]);
else
    obj.chemical_shift{pieces_of_information}=chem_shifts;
    if verbose
        disp(['Found chemical shift :' num2str(obj.chemical_shift{pieces_of_information}) ])
    end
end
%at this point current_line is shorter and only contains one-letter fields:
%example of current_line : L=H31(C25), S=d, J=8.44(H28(C27)), E=1

list_equalsign=strfind(current_line,'=');% delimater determined by position of "," and "="
% check real separator... if not set to zero
for loop_over_fields=1:size(list_equalsign,2)
    %One character code of the field:
    ok=0;
    if loop_over_fields>1%if first
        if strcmp(current_line(list_equalsign(1,loop_over_fields)-2),field_separator(1,1)) ok=1;end% a separator should be found 2 or three caracters earlyer
        if strcmp(current_line(list_equalsign(1,loop_over_fields)-3),field_separator(1,1)) ok=1;end
        if ~ok
            list_equalsign(1,loop_over_fields)=0;
            disp(['not ok ' current_line]);
            
        end
    end
end
list_equalsign=list_equalsign(find(list_equalsign~=0));%removes the zero from the list
%extract fields
for loop_over_fields=1:size(list_equalsign,2)
    code=current_line(list_equalsign(1,loop_over_fields)-1);
    start_pos=list_equalsign(1,loop_over_fields)+1;%the content of the field starts here...
    if loop_over_fields==size(list_equalsign,2)%last
        end_pos=length(current_line);
    else
        for test_backwards=(list_equalsign(1,loop_over_fields+1)-2):-1:1%
            if strcmp(current_line(test_backwards),field_separator(1,1)) end_pos=test_backwards-1;break;
            end
        end
    end
    content=current_line(start_pos:end_pos);
    % now can process the content.
    if code=='E'%integral
        values=sscanf(content,'%f');
        if size(values,1)<1
            error_message=['expect float integral value following "E=" '   current_line];
            disp(['Error :  ' error_message]);
        else
            obj.integral{pieces_of_information}=values;
            if verbose
                disp(['Found integral :' num2str(obj.integral{pieces_of_information}) ])
            end
        end
    end
     if code=='N'%nb atoms...
        values=sscanf(content,'%d');
        if size(values,1)<1
            error_message=['expect int integral value following "E=" '   current_line];
            disp(['Error :  ' error_message]);
        else
            obj.integral{pieces_of_information}=values;
            if verbose
                disp(['Found number of atom N :' num2str(obj.N{pieces_of_information}) ])
            end
        end
    end
    if code=='I'%intensity
        values=sscanf(content,'%f');
        if size(values,1)<1
            error_message=['expect float intensity value following "I=" '   current_line];
            disp(['Error :  ' error_message]);
        else
            obj.intensity{pieces_of_information}=values;
            if verbose
                disp(['Found intensity :' num2str(obj.intensity{pieces_of_information}) ])
            end
        end
    end
    if code=='S'%signal multiplicity
        
        if length(content)<1
            error_message=['expect multiplicity after "S=" '   current_line];
            disp(['Error :  ' error_message]);
        else
            obj.multiplicity{pieces_of_information}=content;
            if verbose
                disp(['Found multiplicity :' num2str(obj.multiplicity{pieces_of_information}) ])
            end
        end
    end
     if code=='L'%signal label

        if length(content)<1
            error_message=['expect label after "L=" '   current_line];
            disp(['Error :  ' error_message]);
        else
            obj.label{pieces_of_information}=content;
            if verbose
                disp(['Found label :' num2str(obj.label{pieces_of_information}) ])
            end
        end
    end
    if code=='J'%couplings
        list_elements_J= strsplit(content,',');
        for loop_over_couplings=1:size(list_elements_J,2)
            curent_J_txt= list_elements_J{1,loop_over_couplings};
            %    disp([num2str(loop_over_couplings) ':'  curent_J_txt]);
            [ coupling_value   out_of_scan ]=sscanf(curent_J_txt,'%f(%*s)');
            
            if size(coupling_value,2)<1
                error_message=['expect coupling after "J=" '   current_line ' content : ' content];
                disp(['Error :  ' error_message]);
            else
                
                obj.J{pieces_of_information,loop_over_couplings}=coupling_value;
                if verbose
                    disp(['Found coupling  J ' num2str(loop_over_couplings) ':' num2str(obj.J{pieces_of_information,loop_over_couplings}) ])
                end
                label_pos_start=strfind(curent_J_txt,'(');
                label_pos_end=strfind(curent_J_txt,')');
                if size(label_pos_end,2)>0
                    label_pos_end=label_pos_end(1,end);
                else
                    label_pos_end=0;
                end
                label='';
                
                if size(label_pos_start,2)>0
                    if size(label_pos_end,2)>0
                        if (label_pos_start(1,1)+1)<=(label_pos_end(1,1)-1)
                            label=curent_J_txt(label_pos_start+1:label_pos_end-1);
                            obj.J_label{pieces_of_information,loop_over_couplings}=label;
                            %           disp(['For coupling number ' num2str(loop_over_couplings) ' found J=' num2str(coupling_value) 'with label:' label])
                            if verbose
                                disp(['Found label for ' num2str(loop_over_couplings) ':' obj.J_label{pieces_of_information,loop_over_couplings} ])
                            end
                        end
                    end
                end
                
            end
        end
        
    end
    
end

end