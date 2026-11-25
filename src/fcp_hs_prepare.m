function hs = fcp_hs_prepare(mesh,refmesh,Aref,Mref,uref,s)
if mesh.dim==1
    nodes=sort([mesh.nodes;refmesh.nodes]);
    nodes=nodes([true;diff(nodes)>1e-12]);
    count=numel(nodes);
    union=fcp_mesh_data(nodes,[(1:count-1).',(2:count).']);
    union.interior=(2:count-1).';
    union.full_to_interior=zeros(count,1);
    union.full_to_interior(union.interior)=(1:count-2).';
    A=fcp_stiffness_1d(union,s);
    points=nodes(union.interior);
    hs.S=A+full(fcp_mass(union));
    hs.Pc=fcp_evaluate(mesh,points);
    hs.ref=fcp_evaluate(refmesh,points)*uref;
else
    points=refmesh.nodes(refmesh.interior,:);
    hs.S=Aref+Mref;
    hs.Pc=fcp_evaluate(mesh,points);
    hs.ref=uref;
end
end
