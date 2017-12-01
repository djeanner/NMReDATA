function draw_structure(obj,opt)
structure=obj.structure;
% Draws the structure
dim1=opt.dim1;
dim2=opt.dim2;
dim3=opt.dim3;


figure(1134);clf;hold on;axis('equal')
%    figure(343);clf;hold on;axis('equal')
for lop=1:size(structure.atom.XYZ,2)
    % plot(structure.atom.XYZ(1,lop),structure.atom.XYZ(2,lop),'k.')
    lab='';
    for lolab=1:size(obj.atom_num,1)
        if size(find(obj.atom_num{lolab}==lop),2)>0
            lab=['' obj.label_signal{1,lolab} ''];
            break;
        end
    end
    atom_txt=structure.atom.n{1,lop};
    if strcmp(atom_txt,'C')
        atom_txt='';
    else
        atom_txt=[atom_txt ' '];
    end
    text_atom=[ atom_txt  lab '['  num2str(lop) ']'];
    text(structure.atom.XYZ(dim1,lop),structure.atom.XYZ(dim2,lop),structure.atom.XYZ(dim3,lop),text_atom,'FontSize',16)
end
for lop=1:size(structure.bond.a1,2)
    % plot(structure.atom.XYZ(1,lop),structure.atom.XYZ(2,lop),'k.')
    l1=structure.bond.a1(1,lop);
    l2=structure.bond.a2(1,lop);
    line_width_used=1;
    if     structure.bond.type(1,lop)>1
        line_width_used=2;
    end
    if     structure.bond.type(1,lop)>3
        line_width_used=3;
    end
    plot3([structure.atom.XYZ(dim1,l1) structure.atom.XYZ(dim1,l2)],[structure.atom.XYZ(dim2,l1) structure.atom.XYZ(dim2,l2)],[structure.atom.XYZ(dim3,l1) structure.atom.XYZ(dim3,l2)],'k-','LineWidth',line_width_used)
end
axis off

drawnow

end