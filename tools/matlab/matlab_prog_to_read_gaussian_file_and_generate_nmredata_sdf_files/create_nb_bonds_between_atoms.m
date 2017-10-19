function nb_bond_between=create_nb_bonds_between_atoms(bond_list)
nb_bond_between=bond_list;
nb_bond_between=bond_list*0+(nb_bond_between>0);


%symmetrize in case it is not
for l1=1:size(nb_bond_between,1)
    for l2=1:size(nb_bond_between,1)
        if nb_bond_between(l1,l2)
            nb_bond_between(l2,l1)=nb_bond_between(l1,l2);
        end
        if nb_bond_between(l2,l1)
            nb_bond_between(l1,l2)=nb_bond_between(l2,l1);
        end
    end
end
for dist_in_bounds=1:100000%number of maximal number of bonds
    count=0;
    for l1=1:size(nb_bond_between,1)
        for l2=1:size(nb_bond_between,1)
            if l1~=l2
                if nb_bond_between(l1,l2)==dist_in_bounds
                    for other=1:size(nb_bond_between,1)
                        if (other~=l1) && (other~=l2)
                            if nb_bond_between(l1,other)==1
                                if nb_bond_between(l2,other)==0
                                    nb_bond_between(l2,other)=dist_in_bounds+1;
                                    nb_bond_between(other,l2)=dist_in_bounds+1;
                                  count=count+1;
                                end
                            end
                            if nb_bond_between(l2,other)==1
                                if nb_bond_between(l1,other)==0
                                    nb_bond_between(l1,other)=dist_in_bounds+1;
                                    nb_bond_between(other,l1)=dist_in_bounds+1;
                                  count=count+1;
                                end
                            end
                        end
                    end
                    
                end
            end
        end
    end
    if count==0 break;end
end


end