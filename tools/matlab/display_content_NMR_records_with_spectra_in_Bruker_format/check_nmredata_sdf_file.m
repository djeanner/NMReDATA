function [super_obj, returned_value, text_of_the_problem]=check_nmredata_sdf_file(file_name,opt)
% reads and verifies the integrity of .nmredata.sdf files
% for format version V 1.0
if isfield(opt,'reading_sdf_file_verbose')
    verbose=opt.reading_sdf_file_verbose;
else
    verbose=0;
end
opt.dim1=1;
opt.dim2=2;
opt.dim3=3;
fig_plot=1;

forceexit=0;
super_obj=[];
tag_identyfyer='NMREDATA_';

returned_value=1;%will be set to zero in case of problem
text_of_the_problem='';
content_of_tag='';

% if nargin==0
%         file_name='./nmredata_sdf_demo_files/benzoapyrene.sdf';
%         file_name='./nmredata_sdf_demo_files/benzocchrysene.sdf';
%         file_name='./nmredata_sdf_demo_files/androsten.sdf';
% end
if ~exist(file_name,'file')
    text_of_the_problem=[text_of_the_problem 'file ' file_name ' not found'];
    disp(text_of_the_problem);
    returned_value=0;
    return
end
if  verbose>0
    disp(['Testing file : ' file_name]);
end
fid = fopen(file_name);

tline = fgetl(fid);
state=0;%=0 when reading the mol part. 1 between tag 2 inside tag
inc=1;
object_counter=1;
mol_block=[];lo_bl=1;
while ischar(tline)
   
    if state==0
        mol_block{lo_bl}= tline;
        lo_bl=lo_bl+1;
    end
    
    if contains(tline,'$$$$') && (state==1)%end of record
        if  verbose>0
            disp(['End of record found... stop reading file ' ]);
        end
        break
    end
    
    if contains(tline,'>  <')% && (state==1)%start of tag
        if (state==2)%when reading tag...
            warning('should insert empty line betweeb tags')
            forceexit=1;
            tag_name=next_tag_name;
        end
        
        next_tag_name=tline(5:end);% skip first 4 char '>  <'
        post=findstr(next_tag_name,'>');%look for the end of the field
        if size(post,1)<1
            returned_value=0;text_of_the_problem=[text_of_the_problem 'No closing caracter in tag label (>)'];
        end
        next_tag_name=next_tag_name(1:post(1,1)-1);%removes the closing ">"
        if length(next_tag_name)<1
            returned_value=0;text_of_the_problem=[text_of_the_problem 'Invalid tag name'];
        end
        if  verbose>0
            disp(['Tag found: ' next_tag_name]);
        end
        state=2;
        first=1;% so that will not store this line
        
    end
    if ((length(tline)==1)||(isempty(tline))) && (state==1)%empty line: end of tag (when reading tag)
        warning(['There should be no emply lines between tags in .sdf files']);
    end
    if (((length(tline)==1)||(isempty(tline))) && (state==2)) || forceexit%empty line: end of tag (when reading tag)
        if ~forceexit
            tag_name=next_tag_name;
        end
        
        forceexit=0;
        %disp(       content_of_tag)
        if  verbose>0
            disp(['End of tag ' tag_name ]);
        end
        if contains(tag_name,tag_identyfyer)%check NMREDATAT_ is in the name
            if  verbose>0
                disp(['Analysing ' tag_name ' tag']);
            end
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            obj=verify_nmredata_tag(tag_name,content_of_tag,tag_identyfyer);
            if contains(obj.tag_name,[tag_identyfyer '1D_'])
                super_obj{object_counter+2}=obj;
                object_counter=object_counter+1;
            end
            if contains(obj.tag_name,[tag_identyfyer '2D_'])
                super_obj{object_counter+2}=obj;
                object_counter=object_counter+1;
                
            end
            if strcmp(obj.tag_name,[tag_identyfyer 'J'])
                super_obj{2}=obj;
                
            end
            if strcmp(obj.tag_name,[tag_identyfyer 'ASSIGNMENT'])
                super_obj{1}=obj;
                super_obj{1}.structure=structure;
            end
            inc=inc+1;
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            
        else
            if  verbose>0
                disp(['Ignoring ' tag_identyfyer ' tag']);
            end
        end
        tag_name=next_tag_name;
        
        state=1;
        content_of_tag='';
        
    end
    if contains(tline,'M  END') && (state==0)%end of the mol part
        if  verbose>0
            disp('End of mol. section');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  verbose>0
            disp('Start analysis of molblock');
        end
        structure=read_mol_block(mol_block,opt);
        if  verbose>0
            disp('Finished analysis of molblock');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if verbose
            disp(['Starts generation of table of 1J, 2J, 3J between atoms ' ])
        end
        structure.nb_bond_between=create_nb_bonds_between_atoms(structure,opt);
        if verbose
            disp(['Ends generates table of 1J, 2J, 3J between atoms ' ])
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        state=1;
        
    end
    drawnow
    
    if verbose==2
        if state==0
            disp(['Mol. part : ' tline])
        end
        if state==2
            disp(['TAG: ' next_tag_name ':'  tline ])
        end
        if state==1
            disp(['Between tags '   tline])
        end
    end
    if state==2
        if ~first
            content_of_tag=[content_of_tag  tline 10];%'\n'
        end
        first=0;
    end
    tline = fgetl(fid);
    
end

if state==2
    returned_value=0;text_of_the_problem=[text_of_the_problem 'End of file before end of TAG'];
end
if state==0
    returned_value=0;text_of_the_problem=[text_of_the_problem 'End of file before end of mol part'];
end
fclose(fid);
if ~returned_value
    % disp(text_of_the_problem)
end

end