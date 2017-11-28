function obj=extract_data_signals(obj,current_line,pieces_of_information,field_separator)

%example of current_line : H1(C2), 8.2433, H1
list_separator=strfind(current_line,field_separator);

if size(list_separator,2)<2%for 3
    error_message=['signals lines expect three fields : '   'inside'  current_line];
    disp(['Error1 :  ' error_message]);
else
    label_or_chemshift =current_line(1:list_separator(1,1)-1);
    chemical_shift=current_line(list_separator(1,1)+1:list_separator(1,2)-1);
    
    obj.label_signal{pieces_of_information}=label_or_chemshift;
    value=sscanf(chemical_shift,'%f');
    obj.chemical_shift{pieces_of_information}=value;
    
  %  disp(['signal: ' label_or_chemshift ' at ' num2str(value) 'ppm corresponds to atoms : '])
    
    for loop_over_atoms=2:size(list_separator,2)
        if loop_over_atoms==size(list_separator,2)
            current_string=current_line(list_separator(1,loop_over_atoms)+1:end);
        else
            current_string=current_line(list_separator(1,loop_over_atoms)+1:list_separator(1,loop_over_atoms+1)-1);
        end
        %    disp(['Atom' num2str(loop_over_atoms-1) ':' current_string ])
        OK=0;
        if strcmp(current_string(1),' ')
            current_string=current_string(2:end);
        end
        if strcmp(current_string(1),' ')
            current_string=current_string(2:end);
        end
        
      %  disp(['here is  <' current_string(1,1) '>'])
        if strcmp(current_string(1),'H')
            current_string(2:end);
            atomnum    =sscanf(current_string(2:end),'%f');
            if size(atomnum)<1
                error_message=['expect a number after H of implicit atom  : ' current_string  ' inside '  current_line];
                disp(['Error2 :  ' error_message]);
                sadf
                
            else
            %    disp(['ImplicitAtom' num2str(loop_over_atoms-1) ': bound to atom ' num2str(atomnum) ])
                ok=1;
            end
            obj.atom_num{pieces_of_information,loop_over_atoms-1}=-atomnum;
            
        else
            atomnum    =sscanf(current_string(1:end),'%f');
            if size(atomnum)<1
                error_message=['expect a number for explicit atom  : ' current_string  ' inside '  current_line];
                disp(['Error3 :  ' error_message]);
                asdf
            else
           %     disp(['Explicit' num2str(loop_over_atoms-1) ': bound to atom ' num2str(atomnum) ])
                ok=1;
            end
            obj.atom_num{pieces_of_information,loop_over_atoms-1}=atomnum;%negative are for implicit hydrogens (directly bound heavy atom)
            
        end
%         if loop_over_atoms==4
%             obj
%             obj.atom_num
%         end
    end
end

end
