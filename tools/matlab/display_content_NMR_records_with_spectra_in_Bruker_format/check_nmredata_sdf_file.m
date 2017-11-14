function [super_obj, returned_value, text_of_the_problem]=check_nmredata_sdf_file(file_name,verbose)
% reads and verifies the integrity of .nmredata.sdf files
% for format version V 0.98
super_obj=[];
tag_identyfyer='NMREDATA_';
if nargin<2
verbose=1;
end
returned_value=1;%will be set to zero in case of problem
text_of_the_problem='';
% if nargin==0    
%         file_name='./nmredata_sdf_demo_files/benzoapyrene.sdf';
%         file_name='./nmredata_sdf_demo_files/benzocchrysene.sdf';
%         file_name='./nmredata_sdf_demo_files/androsten.sdf';
% end
if ~exist(file_name,'file')
    text_of_the_problem=[text_of_the_problem 'file ' file_name ' not found'];
    disp(text_of_the_problem)
    
    returned_value=0;man warning
    
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
    if contains(tline,'M  END') && (state==0)%end of the mol part
        if  verbose>0
            disp('End of mol. section');
        end
        line_num_before_list_atoms=4;
        ret_val= sscanf(mol_block{line_num_before_list_atoms},'%d %d');
        if size(ret_val,1)<2
            disp(['Could not find the number of atoms and bouds in the line ' num2str(line_num_before_list_atoms)  'of the molblock!'])
        else
            structure=[];
            nb_atoms=ret_val(1,1);
            nb_bounds=ret_val(2,1);
            structure.atom.XYZ=zeros(3,nb_atoms);
            
            structure.bond.a1=zeros(1,nb_bounds);
            structure.bond.a2=zeros(1,nb_bounds);
            structure.bond.type=zeros(1,nb_bounds);
            for loop_over_atom_list=1:nb_atoms
                disp(['Found atom line : ' mol_block{line_num_before_list_atoms+1+loop_over_atom_list-1}])
                % work on atoms
                pieces=textscan(mol_block{line_num_before_list_atoms+1+loop_over_atom_list-1},'%f %f %f %s ');
                tt= pieces{1,4};
                tta= tt{1,1};
                structure.atom.XYZ(:,loop_over_atom_list)=[pieces{1,1}(1,1) pieces{1,2}(1,1) pieces{1,3}(1,1)];
                structure.atom.n{1,loop_over_atom_list}=tta;
            end
             for loop_over_bond_list=1:nb_bounds
                disp(['Found bond line : ' mol_block{line_num_before_list_atoms+1+loop_over_bond_list+nb_atoms-1}])
                % work on atoms
                pieces=textscan(mol_block{line_num_before_list_atoms+1+loop_over_bond_list+nb_atoms-1},'%d %d %d ');
                structure.bond.a1(1,loop_over_bond_list)=pieces{1,1}(1,1) ;
                structure.bond.a2(1,loop_over_bond_list)=pieces{1,2}(1,1) ;
                structure.bond.type(1,loop_over_bond_list)=pieces{1,3}(1,1) ;
             end
          %caluclates the number of bond betweeb pairs of atoms
          
% % %           function nb_bond_between=create_nb_bonds_between_atoms(bond_list)
% % %           nb_bond_between=bond_list;
% % %           nb_bond_between=bond_list*0+(nb_bond_between>0);
% % %           
% % %           
% % %           %symmetrize in case it is not
% % %           for l1=1:size(nb_bond_between,1)
% % %               for l2=1:size(nb_bond_between,1)
% % %                   if nb_bond_between(l1,l2)
% % %                       nb_bond_between(l2,l1)=nb_bond_between(l1,l2);
% % %                   end
% % %                   if nb_bond_between(l2,l1)
% % %                       nb_bond_between(l1,l2)=nb_bond_between(l2,l1);
% % %                   end
% % %               end
% % %           end
% % %           for dist_in_bounds=1:100000%number of maximal number of bonds
% % %               count=0;
% % %               for l1=1:size(nb_bond_between,1)
% % %                   for l2=1:size(nb_bond_between,1)
% % %                       if l1~=l2
% % %                           if nb_bond_between(l1,l2)==dist_in_bounds
% % %                               for other=1:size(nb_bond_between,1)
% % %                                   if (other~=l1) && (other~=l2)
% % %                                       if nb_bond_between(l1,other)==1
% % %                                           if nb_bond_between(l2,other)==0
% % %                                               nb_bond_between(l2,other)=dist_in_bounds+1;
% % %                                               nb_bond_between(other,l2)=dist_in_bounds+1;
% % %                                               count=count+1;
% % %                                           end
% % %                                       end
% % %                                       if nb_bond_between(l2,other)==1
% % %                                           if nb_bond_between(l1,other)==0
% % %                                               nb_bond_between(l1,other)=dist_in_bounds+1;
% % %                                               nb_bond_between(other,l1)=dist_in_bounds+1;
% % %                                               count=count+1;
% % %                                           end
% % %                                       end
% % %                                   end
% % %                               end
% % %                               
% % %                           end
% % %                       end
% % %                   end
% % %               end
% % %               if count==0 break;end
% % %           end
% % %           
% % %           
% % %           end
% % %           
         %   super_obj.structure=structure;
        end
        
        
        state=1;
    end
    if contains(tline,'$$$$') && (state==1)%end of record
        if  verbose>0
            disp(['End of record found... stop reading file ' ]);
        end
        break
    end
    if contains(tline,'>  <') && (state==1)%start of tag
        tag_name=tline(5:end);% ignores '>  <'
        post=findstr(tag_name,'>');%look for the end of the field
        if size(post,1)<1
            returned_value=0;text_of_the_problem=[text_of_the_problem 'No closing caracter in tag label (>)'];
        end
        tag_name=tag_name(1:post(1,1)-1);%removes the closing ">"
        if length(tag_name)<1
            returned_value=0;text_of_the_problem=[text_of_the_problem 'Invalid tag name'];
        end
        if  verbose>0
            disp(['Tag found: ' tag_name]);
        end
        state=2;
        content_of_tag='';
        first=1;
    end
    if ((length(tline)==1)||(length(tline)==0)) && (state==2)%empty line: end of tag (when reading tag)
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
        state=1;
    end
    if verbose==2
        if state==0
            disp(['Mol. part : ' tline])
        end
        if state==2
            disp(['TAG: ' tag_name ':'  tline ])
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
    disp(text_of_the_problem)
end

end