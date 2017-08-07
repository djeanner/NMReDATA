function [super_obj, returned_value, text_of_the_problem]=check_nmredata_sdf_file(file_name)
super_obj=[];
tag_identyfyer='NMREDATA_';
verbose=0;
returned_value=1;%will be set to zero in case of problem
text_of_the_problem='';
if nargin==0
    file_name='../../../benzoapyrene/benzoapyrene.sdf';
end
if ~exist(file_name,'file')
    text_of_the_problem=[text_of_the_problem 'file ' file_name ' not found'];
    return
end
if  verbose>0
    disp(['Testing file : ' file_name]);
end
fid = fopen(file_name);

tline = fgetl(fid);
state=0;%=0 when reading the mol part. 1 between tag 2 inside tag
inc=1;
while ischar(tline)
    if contains(tline,'M  END') && (state==0)%end of the mol part
        if  verbose>0
            disp('End of mol. section');
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
        tag_name=tline(5:end);
        post=findstr(tag_name,'>');%look for the end of the field
        if size(post,1)<1
            returned_value=0;text_of_the_problem=[text_of_the_problem 'No closing caracter in tag label (>)'];
        end
        tag_name=tag_name(1:post(1,1)-1);
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
        if contains(tag_name,tag_identyfyer)
            if  verbose>0
                disp(['Analysing ' tag_name ' tag']);
            end
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
            obj=verify_nmredata_tag(tag_name,content_of_tag,tag_identyfyer);
            if strcmp(obj.tag_name,[tag_identyfyer '1D_1H'])
                super_obj{1}=obj;
            end
               if strcmp(obj.tag_name,[tag_identyfyer '1D_13C'])
                super_obj{2}=obj;
               end
            if strcmp(obj.tag_name,[tag_identyfyer '2D_1H_NJ_1H'])
                super_obj{3}=obj;
            end
            if strcmp(obj.tag_name,[tag_identyfyer '2D_13C_1J_1H'])
                super_obj{4}=obj;
            end
            if strcmp(obj.tag_name,[tag_identyfyer '2D_13C_NJ_1H'])
                super_obj{5}=obj;
            end
             if strcmp(obj.tag_name,[tag_identyfyer 'J'])
                super_obj{6}=obj;
             end
             if strcmp(obj.tag_name,[tag_identyfyer 'SIGNALS'])
                super_obj{7}=obj;
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