function         structure= read_mol_block(mol_block,opt)
%% create a structure with atom and bond as fields
% atom will list coordinates in the xyz field and n as the element number
%
% structure.atom.xyz (3,n) x, y, z, coordinates of the atoms
% structure.atom.n   (1,n) x, y, z, the atomic number of the eleemnt
% structure.bond.a1  (1,m) the first atom number of the bond (1..n)
% structure.bond.a2  (1,m)the second atom number of the bond (1..n)
% structure.bond.type(1,m) the type of bond 1:single bond, 2: double bond, etc.
% structure.nb_bonds_between  (a,n,n) if non-zero, the path of a bonds

%dim1,dim2,dim3 inditace which dimension is used for 3D plot.
 structure=[];
    structure.mol_block=mol_block;

if isfield(opt,'dim1')
    dim1=opt.dim1;
else
    dim1=1;
end
if isfield(opt,'dim2')
    dim2=opt.dim2;
else
    dim2=2;
end
if isfield(opt,'dim3')
    dim3=opt.dim3;
else
    dim3=3;
end

if isfield(opt,'draw_verbose')
    verbose=opt.draw_verbose;
else
    verbose=0;
end

line_num_before_list_atoms=4;
ret_val= sscanf(mol_block{line_num_before_list_atoms},'%d %d');
if size(ret_val,1)<2
    error(['Could not find the number of atoms and bouds in the line ' num2str(line_num_before_list_atoms)  'of the molblock!'])
else
   

    nb_atoms=ret_val(1,1);
    nb_bounds=ret_val(2,1);
    structure.type_bond=zeros(nb_atoms,nb_atoms);
    
    structure.atom.XYZ=zeros(3,nb_atoms);
    
    structure.bond.a1=zeros(1,nb_bounds);
    structure.bond.a2=zeros(1,nb_bounds);
    structure.bond.type=zeros(1,nb_bounds);
    for loop_over_atom_list=1:nb_atoms
        if verbose
            disp(['Found atom line : ' mol_block{line_num_before_list_atoms+1+loop_over_atom_list-1}])
        end
        % work on atoms
        pieces=textscan(mol_block{line_num_before_list_atoms+1+loop_over_atom_list-1},'%f %f %f %s ');
        tt= pieces{1,4};
        tta= tt{1,1};
        structure.atom.XYZ(:,loop_over_atom_list)=[pieces{1,1}(1,1) pieces{1,2}(1,1) pieces{1,3}(1,1)];
        structure.atom.n{1,loop_over_atom_list}=tta;
    end
    for loop_over_bond_list=1:nb_bounds
        if verbose
            disp(['Found bond line : ' mol_block{line_num_before_list_atoms+1+loop_over_bond_list+nb_atoms-1}])
        end
        % work on atoms
        pieces=textscan(mol_block{line_num_before_list_atoms+1+loop_over_bond_list+nb_atoms-1},'%d %d %d ');
        structure.bond.a1(1,loop_over_bond_list)=pieces{1,1}(1,1) ;
        structure.bond.a2(1,loop_over_bond_list)=pieces{1,2}(1,1) ;
        structure.bond.type(1,loop_over_bond_list)=pieces{1,3}(1,1) ;
        structure.type_bond(pieces{1,1}(1,1),pieces{1,2}(1,1))=pieces{1,3}(1,1);%0 no bond, 1: single bond, 2: double bond...
        structure.type_bond(pieces{1,2}(1,1),pieces{1,1}(1,1))=pieces{1,3}(1,1);%0 no bond, 1: single bond, 2: double bond...
    end
    
    %caluclates the number of bond betweeb pairs of atoms
    %%%%%
    %%%%%
    %%%%%
   
end

end