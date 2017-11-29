function [nb_bond_between]=create_nb_bonds_between_atoms(structure,opt)
% returns a 3D array:
%
% structure.nb_bonds_between  (a,n,n) 
%
% if non-zero, the path of "a" bonds
% links the two atoms of the second and third arguments
% Important because sometime there is a 1J and a 3J between the same pair of
% atoms (eg. in cyclobutane). Need to know the presence of both the explai
% HSQC & HMBC correlations
%
dim1=opt.dim1;
dim2=opt.dim2;
dim3=opt.dim3;

if isfield(opt,'draw_verbose')
    verbose=opt.draw_verbose;
else
    verbose=0;
end

nb_bond_between=zeros(3,size(structure.type_bond,1),size(structure.type_bond,2));
nb_bond_between(1,:,:)=(structure.type_bond>0)*2;%*2  is to avoid table of logical... but integer...
nb_bond_between=nb_bond_between/2;%/2 is to avoid table of logical... but integer...
nn=size(nb_bond_between,2);

if verbose>1
    figure(1133);clf;hold on;axis('equal')
    %    figure(343);clf;hold on;axis('equal')
    for lop=1:size(structure.atom.XYZ,2)
        % plot(structure.atom.XYZ(1,lop),structure.atom.XYZ(2,lop),'k.')
        text(structure.atom.XYZ(dim1,lop),structure.atom.XYZ(dim2,lop),structure.atom.XYZ(dim3,lop),structure.atom.n{1,lop})
    end
    for lop=1:size(structure.bond.a1,2)
        % plot(structure.atom.XYZ(1,lop),structure.atom.XYZ(2,lop),'k.')
        l1=structure.bond.a1(1,lop);
        l2=structure.bond.a2(1,lop);
        plot3([structure.atom.XYZ(dim1,l1) structure.atom.XYZ(dim1,l2)],[structure.atom.XYZ(dim2,l1) structure.atom.XYZ(dim2,l2)],[structure.atom.XYZ(dim3,l1) structure.atom.XYZ(dim3,l2)],'k-','LineWidth',2)
    end
end
drawnow
sto=0;
%stage 1/2 (main stage) % note that this algorythm is not perfect. it is
%general of any number of bounds but we don't need this. We only need 1-4
%bonds. It may be better to  better to adapt the one of stage 2
for dist_in_bounds=1:2%number of maximal number of bonds
    current=0;
    if verbose>2
        if dist_in_bounds==1 colo='g-';colob='go'; colo2='g--';current=1;liwi=0.5;end
        if dist_in_bounds==2 colo='b-';colob='bo'; colo2='b--';current=1;liwi=0.5;end
        if dist_in_bounds==3 colo='m-';colob='mo'; colo2='m:';current=1;liwi=0.5;end
        if dist_in_bounds==4 colo='c-';colob='co'; colo2='c:';current=1;liwi=0.5;end
    end
    count=0;
    for l1=1:nn
        for l2=1:nn
            %  if l1~=l2
            if l1~=l2
                if nb_bond_between(dist_in_bounds,l1,l2)>0
                    %   for other=1:nn
                    for other=1:nn
                        if (other~=l1) && (other~=l2)
                            for curl=[0 1]
                                if curl
                                    ll1=l1;ll2=l2;
                                else
                                    ll2=l1;ll1=l2;
                                end
                                if nb_bond_between(1,ll1,other)>0
                                    
                                    if sum(nb_bond_between(1:dist_in_bounds,ll2,other))==0
                                        nb_bond_between(dist_in_bounds+1,ll2,other)=ll1;
                                        nb_bond_between(dist_in_bounds+1,other,ll2)=ll1;
                                        
                                        count=count+1;
                                        p1=ll2;
                                        p2=other;
                                        if current
                                            plot3([structure.atom.XYZ(dim1,p1) structure.atom.XYZ(dim1,p2)],[structure.atom.XYZ(dim2,p1) structure.atom.XYZ(dim2,p2)],[structure.atom.XYZ(dim3,p1) structure.atom.XYZ(dim3,p2)],colo,'LineWidth',liwi)
                                            plot3([structure.atom.XYZ(dim1,p1)*0.9+0.1* structure.atom.XYZ(dim1,p2)],[structure.atom.XYZ(dim2,p1) *0.9+0.1*structure.atom.XYZ(dim2,p2)],[structure.atom.XYZ(dim3,p1) *0.9+0.1*structure.atom.XYZ(dim3,p2)],colob,'LineWidth',liwi)
                                        end
                                        
                                        
                                    else
                                        if dist_in_bounds==1% make two bonds
                                            %   nb_bond_between(:,ll2,other)
                                            nb_bond_between(dist_in_bounds+1,ll2,other)=ll1;
                                            nb_bond_between(dist_in_bounds+1,other,ll2)=ll1;
                                            
                                            count=count+1;
                                            % if dist_in_bounds==stored_nb_bond_between{dist_in_bounds-1}(ll2,other)+1
                                            %    disp('close loop');
                                            %  nb_bond_between(l2,other)=99;
                                            %   nb_bond_between(other,l2)=99;
                                            %      disp(num2str(nb_bond_between(dist_in_bounds,ll2,other)))
                                            p1=ll2;
                                            p2=other;
                                            if current
                                                plot3([structure.atom.XYZ(dim1,p1) structure.atom.XYZ(dim1,p2)],[structure.atom.XYZ(dim2,p1) structure.atom.XYZ(dim2,p2)],[structure.atom.XYZ(dim3,p1) structure.atom.XYZ(dim3,p2)],colo2,'LineWidth',4.5)
                                                plot3([structure.atom.XYZ(dim1,p1)*0.9+0.1* structure.atom.XYZ(dim1,p2)],[structure.atom.XYZ(dim2,p1) *0.9+0.1*structure.atom.XYZ(dim2,p2)],[structure.atom.XYZ(dim3,p1) *0.9+0.1*structure.atom.XYZ(dim3,p2)],colob,'LineWidth',4.5)
                                                
                                            end
                                            
                                        end
                                        if dist_in_bounds==2%make 3 bonds...
                                            if ((nb_bond_between(dist_in_bounds,ll2,ll1)~=other))
                                                
                                                %   if  (nb_bond_between(dist_in_bounds-0,ll2,other)+nb_bond_between(dist_in_bounds-2,ll2,other))==0
                                                
                                                nb_bond_between(dist_in_bounds+1,ll2,other)=dist_in_bounds+1;
                                                nb_bond_between(dist_in_bounds+1,other,ll2)=dist_in_bounds+1;
                                                
                                                count=count+1;
                                                % if dist_in_bounds==stored_nb_bond_between{dist_in_bounds-1}(ll2,other)+1
                                                %    disp('close loop');
                                                %  nb_bond_between(l2,other)=99;
                                                %   nb_bond_between(other,l2)=99;
                                                %      disp(num2str(nb_bond_between(dist_in_bounds,ll2,other)))
                                                p1=ll2;
                                                p2=other;
                                                if current
                                                    plot3([structure.atom.XYZ(dim1,p1) structure.atom.XYZ(dim1,p2)],[structure.atom.XYZ(dim2,p1) structure.atom.XYZ(dim2,p2)],[structure.atom.XYZ(dim3,p1) structure.atom.XYZ(dim3,p2)],colo2,'LineWidth',4.5)
                                                    plot3([structure.atom.XYZ(dim1,p1)*0.9+0.1* structure.atom.XYZ(dim1,p2)],[structure.atom.XYZ(dim2,p1) *0.9+0.1*structure.atom.XYZ(dim2,p2)],[structure.atom.XYZ(dim3,p1) *0.9+0.1*structure.atom.XYZ(dim3,p2)],colob,'LineWidth',4.5)
                                                    
                                                end
                                            end
                                            
                                        end
                                        if dist_in_bounds==3
                                            if ((nb_bond_between(dist_in_bounds,ll2,other)==0) && (nb_bond_between(dist_in_bounds-2,ll2,other)==0))
                                                
                                                %              nb_bond_between(:,ll2,other)
                                                %   if  (nb_bond_between(dist_in_bounds-0,ll2,other)+nb_bond_between(dist_in_bounds-2,ll2,other))==0
                                                
                                                nb_bond_between(dist_in_bounds+1,ll2,other)=dist_in_bounds+1;
                                                nb_bond_between(dist_in_bounds+1,other,ll2)=dist_in_bounds+1;
                                                
                                                count=count+1;
                                                % if dist_in_bounds==stored_nb_bond_between{dist_in_bounds-1}(ll2,other)+1
                                                %    disp('close loop');
                                                %  nb_bond_between(l2,other)=99;
                                                %   nb_bond_between(other,l2)=99;
                                                %      disp(num2str(nb_bond_between(dist_in_bounds,ll2,other)))
                                                p1=ll2;
                                                p2=other;
                                                if current
                                                    plot3([structure.atom.XYZ(dim1,p1) structure.atom.XYZ(dim1,p2)],[structure.atom.XYZ(dim2,p1) structure.atom.XYZ(dim2,p2)],[structure.atom.XYZ(dim3,p1) structure.atom.XYZ(dim3,p2)],colo2,'LineWidth',4.5)
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
    
    if count==0
        break;
    end
    
    drawnow
end

%stage 2/2 (for missed bond both 1J and 3J... for e.g. in cyclobutane)
% l1-l3          %first  1J
%....l3-l2       %seconc.1J...
%....   l2-l4  %third. 1J... means 3J l1-l3-l2-other
colo='b:';
for l1=1:nn
    for l2=1:nn
        if (l1~=l2) && (nb_bond_between(1,l1,l2)>0)
            for l3=1:nn
                if (l1~=l3) && (l2~=l3) && (nb_bond_between(1,l2,l3)>0)
                    for l4=1:nn
                        if (l1~=l4) && (l2~=l4) && (l3~=l4) && (nb_bond_between(1,l3,l4)>0)   && (nb_bond_between(3,l1,l4)==0)
                            nb_bond_between(3,l1,l4)=l2+0.001*l3;
                            nb_bond_between(3,l4,l1)=l2+0.001*l3;
                            if verbose>2
                                p1=l1;
                                p2=l4;
                                %   structure.atom.XYZ(dim1,:)
                                plot3([structure.atom.XYZ(dim1,p1) structure.atom.XYZ(dim1,p2)],[structure.atom.XYZ(dim2,p1) structure.atom.XYZ(dim2,p2)],[structure.atom.XYZ(dim3,p1) structure.atom.XYZ(dim3,p2)],colo,'LineWidth',6.5)
                            end
                        end
                    end
                end
            end
        end
    end
end

end