function [table_1h_cs, table_coup_for_return, table_1h_J, label_atoms, ref_to_atom, androsten_better_nmrJ_log, cosy, hsqc, hmbc]=read_dft_nmr_data_fn(dft_folder,min_J)
%only for C-H (no OH !)
table_1h_cs=zeros(0,0);
table_coup_for_return=zeros(0,0);
disp_3d=1;
if disp_3d
    figure(888);clf;hold on;axis('equal');%pbaspect([1 1 1])
end
dist_max_ch=1.6; %maximal distance of CH bond
dist_max_cc=3; %maximal distance of CH bond
folder=dft_folder;

file=[ dft_folder 'log.txtj'];

load(file,'androsten_better_nmrJ_log');
file=[ dft_folder 'log.txt'];

fid=fopen(file);
%    C = textscan(fid, '%s%s%f32%d8%u%f%f%s%f');
C = textscan(fid, '%d%s%f');
fclose(fid);
num=C{1};
atom=C{2};
cs=C{3};

%strcmp(atom(5),'C')
file=[ dft_folder 'log.min.XYZt'];

fid=fopen(file);
%1    C = textscan(fid, '%s%s%f32%d8%u%f%f%s%f');
C = textscan(fid, '%d%f%f%f');
fclose(fid);
xyz=[ C{2} C{3} C{4}];
atom_number=C{1} ;

for looop_of_c=1:size(atom_number,1),%for each C ...
    alias_table{looop_of_c}=num2str(looop_of_c);
    label_atoms{looop_of_c}=alias_table{looop_of_c};
    
end

%% identify methyl
%% also make HSQC table and cosy 2J

hsqc=zeros(size(atom_number,1),size(atom_number,1));
cosy=hsqc*0;

number_of_h=zeros(1,size(atom_number,1));
for looop_of_c=1:size(atom_number,1)%for each C ...
    
    if atom_number(looop_of_c,1)==6%if carbon
        label_atoms{looop_of_c}=['(' alias_table{looop_of_c} ')'];
        for looop_of_h=1:size(atom_number,1)
            if atom_number(looop_of_h,1)==1% if H...
                tmp=xyz(looop_of_c,:)-xyz(looop_of_h,:);
                dist=sum(sum(tmp.*tmp));
                if dist<dist_max_ch
                    number_of_h(  looop_of_c)=  number_of_h(  looop_of_c)+1;
                    label_atoms{looop_of_h}=[alias_table{looop_of_h} '(' alias_table{looop_of_c} ')'];
                    if disp_3d
                        plot3([xyz(looop_of_c,1) xyz(looop_of_h,1)],[xyz(looop_of_c,2) xyz(looop_of_h,2)],[xyz(looop_of_c,3) xyz(looop_of_h,3)],'b-')
                    end
                    hsqc(looop_of_h,looop_of_c)=1;
                    hsqc(looop_of_c,looop_of_h)=1;
                    %% cosy 2J
                    for loop_over_h2=1:size(atom_number,1)
                        if (atom_number(loop_over_h2,1)==1) && (loop_over_h2~=looop_of_h) && hsqc(looop_of_c,loop_over_h2)
                            cosy(loop_over_h2,looop_of_h)=2;
                            cosy(looop_of_h,loop_over_h2)=2;
                            if disp_3d
                                plot3([xyz(looop_of_h,1) xyz(loop_over_h2,1)],[xyz(looop_of_h,2) xyz(loop_over_h2,2)],[xyz(looop_of_h,3) xyz(loop_over_h2,3)],'m:');%2JHH
                            end
                        end
                    end
                end
            end
        end
        
    end
end

%% list carbons connected to carbons (for list of HMBC correlations)
hmbc=hsqc*0;

for looop_of_c=1:size(atom_number,1)%for each C ...
    if atom_number(looop_of_c,1)==6%if carbon
        list_c{looop_of_c}=[];
        label_atoms{looop_of_c}=['(' alias_table{looop_of_c} ')'];
        count=1;
        for looop_of_c2=1:size(atom_number,1)
            if (atom_number(looop_of_c2,1)==6) && (looop_of_c2~=looop_of_c)% if C...
                tmp=xyz(looop_of_c,:)-xyz(looop_of_c2,:);
                dist=sum(sum(tmp.*tmp));
                if dist<dist_max_cc
                    list_c{looop_of_c}=[list_c{looop_of_c} looop_of_c2];
                    if disp_3d
                        plot3([xyz(looop_of_c,1) xyz(looop_of_c2,1)],[xyz(looop_of_c,2) xyz(looop_of_c2,2)],[xyz(looop_of_c,3) xyz(looop_of_c2,3)],'k-')% CC
                    end
                    for loop_over_h1=1:size(atom_number,1)%
                        if hsqc(loop_over_h1,looop_of_c)
                            hmbc(loop_over_h1,looop_of_c2)=2;
                            hmbc(looop_of_c2,loop_over_h1)=2;
                            if disp_3d
                                plot3([xyz(loop_over_h1,1) xyz(looop_of_c2,1)],[xyz(loop_over_h1,2) xyz(looop_of_c2,2)],[xyz(loop_over_h1,3) xyz(looop_of_c2,3)],'m--')% 2JCH
                            end
                            for loop_over_h2=1:size(atom_number,1)%
                                if hsqc(loop_over_h2,looop_of_c2)
                                    cosy(loop_over_h2,loop_over_h1)=3;
                                    cosy(loop_over_h1,loop_over_h2)=3;
                                    if disp_3d
                                        plot3([xyz(loop_over_h1,1) xyz(loop_over_h2,1)],[xyz(loop_over_h1,2) xyz(loop_over_h2,2)],[xyz(loop_over_h1,3) xyz(loop_over_h2,3)],'c:')%3JHH
                                    end
                                end
                            end
                            
                        end
                    end
                end
            end
        end
    end
end
% similar as above...
for looop_of_c=1:size(atom_number,1)%for each C ...
    if atom_number(looop_of_c,1)==6%if carbon
        list_c{looop_of_c}=[];
        label_atoms{looop_of_c}=['(' alias_table{looop_of_c} ')'];
        count=1;
        for looop_of_c2=1:size(atom_number,1)
            if (atom_number(looop_of_c2,1)==6) && (looop_of_c2~=looop_of_c)% if C...
                tmp=xyz(looop_of_c,:)-xyz(looop_of_c2,:);
                dist=sum(sum(tmp.*tmp));
                if dist<dist_max_cc
                    list_c{looop_of_c}=[list_c{looop_of_c} looop_of_c2];
                    if disp_3d
                        plot3([xyz(looop_of_c,1) xyz(looop_of_c2,1)],[xyz(looop_of_c,2) xyz(looop_of_c2,2)],[xyz(looop_of_c,3) xyz(looop_of_c2,3)],'k-')% CC
                    end
                    for loop_over_h1=1:size(atom_number,1)%
                        if hsqc(loop_over_h1,looop_of_c)
                            
                            %%for 3JCH hmbc
                            if size(list_c,2)>=looop_of_c2
                                if size(list_c{looop_of_c2},2)>0
                                    list_c{looop_of_c2}
                                    for loop_over_next_c=list_c{looop_of_c2}
                                        if loop_over_next_c~=looop_of_c
                                            hmbc(loop_over_h1,loop_over_next_c)=3;
                                            hmbc(loop_over_next_c,loop_over_h1)=3;
                                            if disp_3d
                                                plot3([xyz(loop_over_h1,1) xyz(loop_over_next_c,1)],[xyz(loop_over_h1,2) xyz(loop_over_next_c,2)],[xyz(loop_over_h1,3) xyz(loop_over_next_c,3)],'c--')%3JCH
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%
% list_c
% disp('-------------')
% for looop_of_c=1:size(list_c,2)%for each C ...
%     disp(num2str(looop_of_c))
%     list_c{looop_of_c}
% end
%
% disp('-------------')
%% fill hmbc table with number of bond between any pair of C and H
for looop_of_c=1:size(atom_number,1)%for each C ...
    if atom_number(looop_of_c,1)==6%if carbon
        list_c{looop_of_c}=[];
        label_atoms{looop_of_c}=['(' alias_table{looop_of_c} ')'];
        count=1;
        for looop_of_c2=1:size(atom_number,1)
            if (atom_number(looop_of_c2,1)==6) && (looop_of_c2~=looop_of_c)% if C...
                tmp=xyz(looop_of_c,:)-xyz(looop_of_c2,:);
                dist=sum(sum(tmp.*tmp));
                if dist<dist_max_cc
                    list_c{looop_of_c}=[list_c{looop_of_c} looop_of_c2];
                    
                end
                
            end
        end
    end
end


%% change lablels if file exists
file=[ dft_folder 'log.alias'];
if exist(file),
    fid=fopen(file);
    %    C = textscan(fid, '%s%s%f32%d8%u%f%f%s%f');
    C = textscan(fid, '%d%s');
    fclose(fid);
    numin=C{1};
    aliin=C{2};
    for loodpf=1:size(numin,1),
        tttpl=aliin(loodpf);
        label_atoms{numin(loodpf)}=tttpl{1};
    end
end


for looop_of_c=1:size(atom_number,1)%for each C ...
    if disp_3d
        text(xyz(looop_of_c,1) ,xyz(looop_of_c,2) ,xyz(looop_of_c,3) ,label_atoms(looop_of_c))
    end
end

number_of_h_toth=(sum(number_of_h));

number_of_different_h=sum(1*(number_of_h==1))+sum(2*(number_of_h==2))+sum(0001*(number_of_h==3));
table_of_shift_ch3=zeros(1,size(atom_number,1));
table_1h_cs=zeros(1,number_of_different_h);
table_1h_J=zeros(number_of_different_h,size(atom_number,1));
loop_over_1H=1;
loop_over_all_1H=1;
text2='';
count_of_c=0;
other_count_of_h=0;
table_ok=[];

for looop_of_c=1:size(atom_number,1),%for each C ...
    if atom_number(looop_of_c,1)==1,
        
        other_count_of_h=other_count_of_h+1;
    end
    if atom_number(looop_of_c,1)==6,
        count_table_of_shift_ch3=1;%table_of_shift_ch3
        textb='';
        
        count_of_c=count_of_c+1;
        textb=[textb ' C'];
        here=number_of_h(looop_of_c);%1 for ch 2 for CH" 3 for CH3
        if here==0,
            textb=[textb 'q '];
        end
        if here==1,
            textb=[textb 'H '];
        end
        if here==2,
            textb=[textb 'H2'];
        end
        if here==3,
            textb=[textb 'H3'];
        end
        textb=[textb '(' num2str(count_of_c,'%3d') ')'];
        
        text2=[text2 textb ' '];
        
        for looop_of_h=1:size(atom_number,1),%for each C ...
            if atom_number(looop_of_h,1)==1,
                
                tmp=xyz(looop_of_c,:)-xyz(looop_of_h,:);
                dist=sum(sum(tmp.*tmp));
                if dist<dist_max_ch,
                    
                    if here==3%for methyl sum (to average)
                        if count_table_of_shift_ch3==1
                            %  local=loop_over_all_1H;
                            % local=other_count_of_h;
                            local=looop_of_h;
                            table_ok=[table_ok looop_of_h];
                            
                        else
                            %table_of_shift_ch3(count_table_of_shift_ch3)=local;
                            %   table_of_shift_ch3(loop_over_all_1H)=local;
                            table_of_shift_ch3(looop_of_h)=local;
                            
                        end
                        count_table_of_shift_ch3=count_table_of_shift_ch3+1;%table_of_shift_ch3
                        
                        table_1h_cs(1,loop_over_1H)=table_1h_cs(1,loop_over_1H)+cs(looop_of_h,1);
                        text2=[text2 num2str(cs(looop_of_h,1),'%.2f') '-'];
                        num_h=1;
                        for looop_of_h2=1:size(atom_number,1),%for each C ...
                            %  if atom_number(looop_of_h2,1)==1,
                            table_1h_J(loop_over_1H , num_h)=table_1h_J(loop_over_1H , num_h)+androsten_better_nmrJ_log(looop_of_h,looop_of_h2);
                            num_h=num_h+1;
                            %  end
                        end
                        loop_over_all_1H=loop_over_all_1H+1;
                    else
                        table_ok=[table_ok looop_of_h];
                        
                        table_1h_cs(1,loop_over_1H)=cs(looop_of_h);
                        text2=[text2 num2str(table_1h_cs(1,loop_over_1H),': %.2f') ','];
                        laelsh{loop_over_1H}=[textb num2str(cs(looop_of_c,1),'%.3f') ' '];
                        
                        num_h=1;
                        for looop_of_h2=1:size(atom_number,1),%for each C ...
                            % if atom_number(looop_of_h2,1)==1,
                            table_1h_J(loop_over_1H , num_h)=table_1h_J(loop_over_1H , num_h)+androsten_better_nmrJ_log(looop_of_h,looop_of_h2);
                            num_h=num_h+1;
                            % end
                        end
                        ref_to_atom(loop_over_1H,1)=looop_of_h;
                        
                        loop_over_1H=loop_over_1H+1;
                        loop_over_all_1H=loop_over_all_1H+1;
                        
                    end
                    
                end
            end
        end
        if here==3
            table_1h_cs(loop_over_1H)=table_1h_cs(loop_over_1H)/3;%average for 3 of the mythel
            text2=[text2 num2str(table_1h_cs(loop_over_1H),': %.2f') ','];
            laelsh{loop_over_1H}=[textb num2str(cs(looop_of_c,1),'%.3f') ' '];
            
            num_h=1;
            for looop_of_h2=1:size(atom_number,1),%for each H ...
                % if atom_number(looop_of_h2,1)==1,
                table_1h_J(loop_over_1H , num_h)=table_1h_J(loop_over_1H , num_h)/3;
                num_h=num_h+1;
                %  end
            end
            ref_to_atom(loop_over_1H,1)=looop_of_c;
            
            loop_over_1H=loop_over_1H+1;
            
        end
        %dist=dist+100*(dist==0);
        %  min(min(dist))
    end
    
end




%table_of_shift_ch3;
table_ok2=[];
for loo1=1:size(atom_number,1),%  number_of_different_h
    if table_of_shift_ch3(loo1)==0,
        if atom_number(loo1,1)==1,
            table_ok2=[table_ok2 loo1];
        end
    end
end
for loo1=1:number_of_different_h,%  number_of_different_h
    
    %  for loo2=1:size(atom_number,1),%  number_of_different_h
    %  if   abs(table_1h_J(loo1,loo2))>min_J,
    % laelsh{loo1}= [ laelsh{loo1} ' ' num2str(table_1h_J(loo1,loo2),'%.1f') ' '];
    %    if atom_number(loo2,1)==1,
    %   laelsh{loo1}= [ laelsh{loo1}  '(' num2str(loo1) ',' num2str(loo2) ')' num2str(table_1h_J(loo1,loo2),'%.1f') ' '];
    
    %     text2=[text2 '(' num2str(loo1) ',' num2str(loo2) ')' num2str(table_1h_J(loo1,loo2),'%.1f') ' '];
    % end
    % end
    % end
    %for loo1=1:size(atom_number,1),%
    if table_of_shift_ch3(loo1)~=0,
        work_on=table_of_shift_ch3(loo1);
        for loo2=1:size(atom_number,1),%
            if work_on==table_of_shift_ch3(loo2),
                for loo3=1:number_of_different_h,
                    table_1h_J(loo3,work_on)=table_1h_J(loo3,work_on)+table_1h_J(loo3,loo2);
                end
                table_of_shift_ch3(loo2)=0;
            end
            
        end
        for loo3=1:number_of_different_h,
            table_1h_J(loo3,work_on)=table_1h_J(loo3,work_on)/3;
            %  disppp=table_1h_J(loo3,work_on)
        end
    end
end

table_1h_J=table_1h_J(:,table_ok);




%text
text2='';
table_coup_for_return=zeros(number_of_different_h,number_of_different_h);
for loo1=1:number_of_different_h,
    inc2=1;
    text2=[text2 laelsh{loo1}] ;
    text2=[text2 num2str(table_1h_cs(1,loo1),'%.3f') ' J= '];
    
    for loo2=1:number_of_different_h,
        if abs(table_1h_J(loo1,loo2))>min_J,
            % text2=[text2 num2str(table_1h_J(loo1,loo2),'%.1f') ' '];
            text2=[text2 '(' num2str(loo1) ',' num2str(loo2) ') ' num2str(table_1h_J(loo1,loo2),'%.1f') ' '];
            %    text2=[text2 '(' num2str(loo2) ',' num2str(loo1) ')' num2str(table_1h_J(loo2,loo1),'%.1f') ' '];
            table_coup_for_return(loo1,inc2)=table_1h_J(loo1,loo2);
            inc2=inc2+1;
        end
    end
    
    text2=[text2 '\n'];
    
end
%
% end
read_dft_nmr_data_fn_disp= sprintf(text2)


table_coup_for_return=abs(table_coup_for_return);

if disp_3d
    savefig('Structure_with_COSY_HSQC_HMBC_correlations.fig')
end
end
