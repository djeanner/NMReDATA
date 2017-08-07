function               red_system=reduce_system(system,list_in)

number_of_items=size(system.label_atoms,2);
first=list_in(1,1);
list=[];
for lo=1:number_of_items
    if size(find(list_in(1,2:end)== lo),2)==0
        list=[list lo];
    end
end

%% copy most of it...
red_system.xyz=system.xyz(list,:);
red_system.atom_number=system.atom_number(list,1);
red_system.J=system.J(list,list);
red_system.cs=system.cs(list,1);
red_system.bond_list=system.bond_list(list,list);

%red_system.label_atoms{1,1:size(red_system.atom_number,2)}=system.label_atoms{1,list};
inc=1;
for loo=list
    red_system.label_atoms{1,inc}=system.label_atoms{1,loo};
    red_system.lookup{1,inc}=system.lookup{1,loo};
    inc=inc+1;
end
for lo2=list_in(1,2:end)
        red_system.lookup{first}=[red_system.lookup{1,first}  system.lookup{lo2}];

    for lo33=lo2+1
   % red_system.lookup{lo33}=system.lookup{lo33}-1;
    end
end

%% copy the left over averaging
%shoudl be 1 anyways not needed

%red_system.atom_number(first,1)=system.atom_number(first,1);
red_system.xyz(first,:)=sum((system.xyz(list_in,:)))/size(list_in,2);
red_system.cs(first,1)=sum(sum(system.cs(list_in,1)))/size(list_in,2);
red_system.J(:,first)=sum((system.J(list,list_in)'))/   (   size(list_in,2) );
red_system.J(first,:)=sum((system.J(list_in,list)))/   (   size(list_in,2) );
red_system.J(first,first)=sum(sum(system.J(list_in,list_in)))/   (   (size(list_in,2)*size(list_in,2))   -size(list_in,2) );%diag elements =0 so not counted for average

% fasdf
% red_system.bond_list=system.bond_list(list_in,list_in);
% red_system.J=system.J(list_in,list_in);
% red_system.label_atoms=system.label_atoms{1,list_in};
%
%
%
%
% xyz=zeroes(number_of_items,3);
% atom_number=zeroes(number_of_items,1);
% cs=zeroes(number_of_items,1);
% J=zeros(number_of_items,number_of_items);
% bond_list=zeros(number_of_items,number_of_items);
% %     atom_number: [36?1 int32]
% %               J: [36?36 double]
% %              cs: [36?1 double]
% %       bond_list: [36?36 double]
% %     label_atoms: {1?36 cell}
% inc=1;
% found=0;
% for loop=1:size(system.label_atoms,2)%loo over all
%     if size(find(list==loop),2)
%         if found==0
%
%             inc=inc+1;
%         end
%     else
%         inc=inc+1;
%     end
% end
end