function obj=extract_data_2d(obj,current_line,pieces_of_information,field_separator)

%example of current_line :H28(C27)/H30(C29), Ja=6.13, J1=7.85(H31(C25)), J2=7.40(H32(C26))
%extract first field, the chemical shift or labels separated by "/"
%separated by "-" (eg. 7.5-7.2)
list_separator=strfind(current_line,field_separator);
if size(list_separator,1)>0
    pair_of_signals=current_line(1:list_separator(1,1)-1);current_line=current_line(list_separator(1,1)+1:end);
else
    pair_of_signals=current_line;current_line='';%no other field
end
position_of_separator_of_the_two=strfind(pair_of_signals,'/');
% first_label_or_chemshift=pair_of_signals(1:
% =sscanf(pair_of_signals,'%f/%f');%this returns either one or two values.
if size(position_of_separator_of_the_two,1)<1
    error_message=['2 expect separator slash in : '  pair_of_signals ' in line: '  current_line];
    disp(['Error :  ' error_message]);
    error(['Error :  ' error_message]);
else
    first_label_or_chemshift =pair_of_signals(1:position_of_separator_of_the_two(1,1)-1);
    second_label_or_chemshift=pair_of_signals(position_of_separator_of_the_two(1,1)+1:end);
  %  disp(['correlation between ' first_label_or_chemshift ' and ' second_label_or_chemshift])
        obj.signalf1{pieces_of_information}=first_label_or_chemshift;
        obj.signalf2{pieces_of_information}=second_label_or_chemshift;

end
%at this point current_line is shorter and only contains one-letter fields:
%example of current_line : Ja=6.13, J1=7.85(H31(C25)), J2=7.40(H32(C26))

list_equalsign=strfind(current_line,'=');% delimater determined by position of "," and "="
% check real separator... if not set to zero
for loop_over_fields=1:size(list_equalsign,2)
    %One character code of the field:
    ok=0;
    if loop_over_fields>1%if first
        if strcmp(current_line(list_equalsign(1,loop_over_fields)-2),field_separator(1,1)) ok=1;end% a separator should be found 2 or three or four caracters earlyer
        if strcmp(current_line(list_equalsign(1,loop_over_fields)-3),field_separator(1,1)) ok=1;end
        if strcmp(current_line(list_equalsign(1,loop_over_fields)-4),field_separator(1,1)) ok=1;end
        if ~ok
            list_equalsign(1,loop_over_fields)=0;
            disp(['not ok ' current_line]);
            asdf
        end
    end
end
list_equalsign=list_equalsign(find(list_equalsign~=0));%removes the zero from the list
%extract fields
for loop_over_fields=1:size(list_equalsign,2)
    code=current_line(list_equalsign(1,loop_over_fields)-1);
    code2=current_line(list_equalsign(1,loop_over_fields)-2);
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
        end
    end
     if code=='I'%intensity
        values=sscanf(content,'%f');
        if size(values,1)<1
            error_message=['expect float intensity value following "E=" '   current_line];
            disp(['Error :  ' error_message]);
        else
            obj.intensity{pieces_of_information}=values;
        end
    end
    
    if code2=='J'%couplings threee kinds active (Ja) passive F1; (J1) passive F2 (J2)
        if code=='a'
                [ coupling_value   out_of_scan ]=sscanf(content,'%f');
                if size(coupling_value,2)<1
                    error_message=['expect coupling after "Ja=" '   current_line];
                    disp(['Error :  ' error_message]);
                else
                    obj.Ja{1,pieces_of_information}=coupling_value;
       %              disp(['For active coupling  found J=' num2str(coupling_value) ])
                end
        end
        if code=='1'
            list_elements_J= strsplit(content,field_separator);
            for loop_over_couplings=1:size(list_elements_J,2)
                curent_J_txt= list_elements_J{1,loop_over_couplings};
                %    disp([num2str(loop_over_couplings) ':'  curent_J_txt]);
                [ coupling_value   out_of_scan ]=sscanf(curent_J_txt,'%f(%*s)');
                if size(coupling_value,2)<1
                    error_message=['expect coupling after "J1=" '   current_line];
                    disp(['Error :  ' error_message]);
                else
                  obj.J1{pieces_of_information,loop_over_couplings}=coupling_value;
                    label_pos_start=strfind(curent_J_txt,'(');
                    label_pos_end=strfind(curent_J_txt,')');
                    label_pos_end=label_pos_end(1,end);
                    label='';
                    if size(label_pos_start,2)>0
                        if size(label_pos_end,2)>0
                            if (label_pos_start(1,1)+1)<=(label_pos_end(1,1)-1)
                                label=curent_J_txt(label_pos_start+1:label_pos_end-1);
                                obj.J1_label{pieces_of_information,loop_over_couplings}=label;
                     %                      disp(['For coupling number ' num2str(loop_over_couplings) ' found J1=' num2str(coupling_value) 'with label:' label])
                            end
                        end
                    end
                end
            end
        end
         if code=='2'
            list_elements_J= strsplit(content,field_separator);
            for loop_over_couplings=1:size(list_elements_J,2)
                curent_J_txt= list_elements_J{1,loop_over_couplings};
                %    disp([num2str(loop_over_couplings) ':'  curent_J_txt]);
                [ coupling_value   out_of_scan ]=sscanf(curent_J_txt,'%f(%*s)');
                if size(coupling_value,2)<1
                    error_message=['expect coupling after "J2=" '   current_line];
                    disp(['Error :  ' error_message]);
                else
                  obj.J2{pieces_of_information,loop_over_couplings}=coupling_value;
                    label_pos_start=strfind(curent_J_txt,'(');
                    label_pos_end=strfind(curent_J_txt,')');
                    label_pos_end=label_pos_end(1,end);
                    label='';
                    if size(label_pos_start,2)>0
                        if size(label_pos_end,2)>0
                            if (label_pos_start(1,1)+1)<=(label_pos_end(1,1)-1)
                                label=curent_J_txt(label_pos_start+1:label_pos_end-1);
                                obj.J2_label{pieces_of_information,loop_over_couplings}=label;
                         %                  disp(['For coupling number ' num2str(loop_over_couplings) ' found J2=' num2str(coupling_value) 'with label:' label])
                            end
                        end
                    end
                end
            end
        end
    end
end

end
