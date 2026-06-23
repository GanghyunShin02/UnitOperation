import numpy as np
import matplotlib.pyplot as plt
import dolfinx
from mpi4py import MPI
import gmsh
import ufl
from ufl import grad,TestFunction,TrialFunction,dx
from dolfinx.fem import Constant, Function,functionspace,dirichletbc,locate_dofs_topological,assemble_scalar
from dolfinx import mesh,fem,io
from dolfinx.io import VTXWriter,gmsh as gmshio
from pathlib import Path
from CoolProp.CoolProp import PropsSI
from basix.ufl import element
from petsc4py import PETSc

# 메시생성
gmsh.initialize()
gdim=3

if MPI.COMM_WORLD.rank==0:
    #꽉찬 원통,유체를 생성
    inner_cylinder=gmsh.model.occ.addCylinder(0,0,2,0,0,10,1)
    outer_cylinder=gmsh.model.occ.addCylinder(0,0,2,0,0,10,2)
    fluid_in=gmsh.model.occ.addCylinder(0,0,2,0,0,10,0.9)
    fluid_outt=gmsh.model.occ.addCylinder(0,0,2,0,0,10,1.9)

    shape=[(gdim,inner_cylinder),(gdim,outer_cylinder),(gdim,fluid_in),(gdim,fluid_outt)]
    out,outmap=gmsh.model.occ.fragment(shape,[])
    gmsh.model.occ.synchronize()
    print(f'out{out}')
    print(f'outmap{outmap}')

    fluid_in_tag=fluid_in
    wall_in_tag=outmap[0][0][1]
    fluid_out_tag=outmap[3][0][1]
    wall_out_tag=outmap[1][0][1]

    '''
    fragment로 겹치는걸 다 자름..? 
    out에 대입한 객체의 순서대로 outmap에서 반환
    out[(3, 3), (3, 4), (3, 5), (3, 6)]
    outmap[[(3, 4), (3, 3)], [(3, 5), (3, 6), (3, 4), (3, 3)], [(3, 3)], [(3, 6), (3, 4), (3, 3)]]
    '''


    '''
    # 컷하면 도형들이 리스트안에 튜플로 존재함..
    fluid_out([(3, 5)], [[(3, 5)], []])
    wall_in([(3, 1)], [[(3, 1)], []])
    wall_out([(3, 2)], [[(3, 2)], []])
    두벚ㄴ째 리스트는 outdimtagsmap으로 잘라낸 도형의 tag
    # 컷으로 도형두개 생기면 튜플이 두개.. 이런식, 튜플의 첫번째원소는 차원수 두번째는 tag.
    ,_ 튜플 언팩킹으로 첫번째 리스트만 쓴다..

    컷하면 자꾸 에러남.fragment로 진행
    '''
  


    gmsh.model.addPhysicalGroup(3,[fluid_in_tag],1,'fluid_in')
    gmsh.model.setPhysicalName(3,fluid_in_tag,'fluid_in')
    gmsh.model.addPhysicalGroup(3,[wall_in_tag],2,'wall_in')
    gmsh.model.setPhysicalName(3,wall_in_tag,'wall_in')
    gmsh.model.addPhysicalGroup(3,[fluid_out_tag],3,'fluid_out')
    gmsh.model.setPhysicalName(3,fluid_out_tag,'fluid_out')
    gmsh.model.addPhysicalGroup(3,[wall_out_tag],4,'wall_out')
    gmsh.model.setPhysicalName(3,wall_out_tag,'wall_out')

    gmsh.option.setNumber('Mesh.Algorithm3D',1)
    gmsh.option.setNumber('Mesh.MeshSizeMin',0.05)
    gmsh.option.setNumber('Mesh.MeshSizeMax',0.1)
    gmsh.option.setNumber('Geometry.Tolerance',1e-6)
    gmsh.model.mesh.generate(gdim)

#print(f'modeltomesh{gmshio.model_to_mesh(gmsh.model,MPI.COMM_WORLD,0,gdim)}')
# meshdata는 meshio.Mesh객체, meshdata.mesh는 dolfinx.mesh.Mesh객체, meshdata.facet_tags는 facet_marker, meshdata.cell_tags는 cell_marker
meshdata=gmshio.model_to_mesh(gmsh.model,MPI.COMM_WORLD,0,gdim)
domain=meshdata.mesh
facet_marker=meshdata.facet_tags
cell_marker=meshdata.cell_tags
gmsh.finalize()

with io.XDMFFile(domain.comm, "heat_exchanger_geometry.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    domain.topology.create_connectivity(domain.topology.dim, domain.topology.dim)
    xdmf.write_meshtags(cell_marker, domain.geometry)



'''
일단 병류흐름부터. 질량유량은 둘다 10kg/s로.
내관유체 초기온도는 300K, 외관유체 초기온도는 350K.  관의 초기온도도 동일하게..
'''
r_i,r_O,R_i,R_O=0.9,1.0,1.9,2.0
R_bar=(R_i-r_O)/2 #수력학적 반지름.
A_i=np.pi*r_i**2
A_O=np.pi*(R_i**2-r_O**2)
m=10 #kg/s

def rho_T(T):
    return PropsSI('D','T',T,'P',101325,'Water')
def mu_T(T):
    return PropsSI('V','T',T,'P',101325,'Water')
def nu_T(T):
    return PropsSI('V','T',T,'P',101325,'Water')/rho_T(T)
def k_T(T):
    return PropsSI('L','T',T,'P',101325,'Water')
def cp_T(T):
    return PropsSI('C','T',T,'P',101325,'Water')

# 초기속도
u_i_0=(m/(rho_T(300)*A_i),0,0)
u_O_0=(m/(rho_T(350)*A_O),0,0)

# 속도장계싼에서쓸 물성치
rho_i_0=rho_T(300)
rho_O_0=rho_T(350)
nu_i_0=nu_T(300)
nu_O_0=nu_T(350)
k_i_0=k_T(300)
k_O_0=k_T(350)
cp_i_0=cp_T(300)
cp_O_0=cp_T(350)
mu_i_0=mu_T(300)
mu_O_0=mu_T(350)

# 함수공간
v_cg_1=element('Lagrange',mesh.basix_cell(),2,shape=(mesh.geometry.dim,)) #내관
P_cg_1=element('Lagrange',mesh.basix_cell(),1) #내관
v_cg_2=element('Lagrange',mesh.basix_cell(),2,shape=(mesh.geometry.dim,)) # 외관
P_cg_2=element('Lagrange',mesh.basix_cell(),1) # 외관

V_in=ufl.FunctionSpace(domain,v_cg_1)
V_out=ufl.FunctionSpace(domain,v_cg_2)
P_in=ufl.FunctionSpace(domain,P_cg_1)
P_out=ufl.FunctionSpace(domain,P_cg_2)



# 경계조건.. 관벽에서 속도가 0으로
u_noslip=np.array((0,)*mesh.geometry.dim,dtype=PETSc.ScalarType)
bcu_wall_in=dirichletbc(u_noslip,fem.locate_dofs_topological(V_in,domain.topology.dim-1,facet_marker.find(wall_in_tag)),V_in)
bcu_wall_out=dirichletbc(u_noslip,fem.locate_dofs_topological(V_out,domain.topology.dim-1,facet_marker.find(wall_out_tag)),V_out)
bcu_inlet=dirichletbc(u_i_0,fem.locate_dofs_topological(V_in,domain.topology.dim-1,facet_marker.find(fluid_in_tag)),V_in)
bcu_outlet=dirichletbc(u_O_0,fem.locate_dofs_topological(V_out,domain.topology.dim-1,facet_marker.find(fluid_out_tag)),V_out)
bcu=[bcu_wall_in, bcu_wall_out, bcu_inlet, bcu_outlet]

# 튜토리얼에서는 아웃렛 압력을 0으로고정했는데 일단 그냥해보겠음.


# 이산화한 NS식을 약형식으로..

u_in,v_in=TrialFunction(V_in),TestFunction(V_in)
u_n=ufl.Function(V_in)
u_n1=ufl.Function(V_in)
u_s=ufl.Function(V_in)
p_in,q_in=TrialFunction(P_in),TestFunction(P_in)

u_out,v_out=TrialFunction(V_out),TestFunction(V_out)
p_out,q_out=TrialFunction(P_out),TestFunction(P_out)

dt=fem.Constant(domain,1)

# 1.속도추정
F_in=(rho_i_0/dt)*ufl.dot((u_in-u_n),v_in)*dx
F_in+=rho_i_0*ufl.inner(ufl.grad(1.5*u_n-0.5*u_n1),0.5*ufl.nabla_grad(u_in-u_n),v_in)*dx
F_in+=0.5*mu_i_0*ufl.inner(ufl.nabla_grad(u_in+u_n),ufl.nabla_grad(v_in))*dx
F_in+=ufl.inner(ufl.nabla_grad())

a_in=ufl.lhs(F_in)
L_in=ufl.rhs(F_in)
A_in=fem.petsc.create_vector(a_in)
b_in=fem.petsc.create_vector(ufl.extract_function_space(L_in)) #?

F_out=(rho_O_0/dt)*ufl.dot((u_out-u_n),v_out)*dx
F_out+=rho_O_0*ufl.inner(ufl.grad(1.5*u_n-0.5*u_n1),0.5*ufl.nabla_grad(u_out-u_n),v_out)*dx
F_out+=0.5*mu_O_0*ufl.inner(ufl.nabla_grad(u_out+u_n),ufl.nabla_grad(v_out))*dx
a_out=ufl.lhs(F_out)
L_out=ufl.rhs(F_out)
A_out=fem.petsc.create_vector(a_out)
b_out=fem.petsc.create_vector(ufl.extract_function_space(L_out)) #?

#2. 압력보정

a_p_in=ufl.form(ufl.dot(grad(p_in),grad(q_in))*dx)
L_p_in=ufl.form(-rho_i_0/dt*ufl.dot(ufl.div(u_s)),ufl.grad(q_in)*dx)
A_p_in=fem.assemble_matrix(a_p_in,bcs=[]) # 압력 경계조건...
A_p_in.assemble()
b_p_in=fem.petsc.create_vector(ufl.extract_function_space(L_p_in))

a_p_out=ufl.form(ufl.dot(grad(p_out),grad(q_out))*dx)
L_p_out=ufl.form(-rho_O_0/dt*ufl.dot(ufl.div(u_s)),ufl.grad(q_out)*dx)
A_p_out=fem.assemble_matrix(a_p_out,bcs=[]) # 압력 경계조건...
A_p_out.assemble()
b_p_out=fem.petsc.create_vector(ufl.extract_function_space(L_p_out))

#3. 속도보정
ac=ufl.form(rho_i_0*ufl.dot(u_in,v_in)*dx)
Lc=ufl.form(rho_i_0*ufl.dot(u_s,v_in)*dx-rho_i_0*dt*ufl.dot(grad(p_in),v_in)*dx)
Ac=fem.assemble_matrix(ac,bcs=bcu) # 속도 경계조건
Ac.assemble()
b_c=fem.petsc.create_vector(ufl.extract_function_space(Lc))

acc=ufl.form(rho_O_0*ufl.dot(u_out,v_out)*dx)
Lcc=ufl.form(rho_O_0*ufl.dot(u_s,v_out)*dx-rho_O_0*dt*ufl.dot(grad(p_out),v_out)*dx)
Acc=fem.assemble_matrix(acc,bcs=bcu) # 속도 경계조건
Acc.assemble()
b_cc=fem.petsc.create_vector(ufl.extract_function_space(Lcc))

# Solver for step 1
solver1 = PETSc.KSP().create(mesh.comm)
solver1.setOperators(A_in)
solver1.setType(PETSc.KSP.Type.BCGS)
pc1 = solver1.getPC()
pc1.setType(PETSc.PC.Type.JACOBI)

solver1_1 = PETSc.KSP().create(mesh.comm)
solver1_1.setOperators(A_out)
solver1_1.setType(PETSc.KSP.Type.BCGS)
pc1_1 = solver1_1.getPC()
pc1_1.setType(PETSc.PC.Type.JACOBI)

# Solver for step 2
solver2 = PETSc.KSP().create(mesh.comm)
solver2.setOperators(A_p_in)
solver2.setType(PETSc.KSP.Type.MINRES)
pc2 = solver2.getPC()
pc2.setType(PETSc.PC.Type.HYPRE)
pc2.setHYPREType("boomeramg")

solver2_1 = PETSc.KSP().create(mesh.comm)
solver2_1.setOperators(A_p_out)
solver2_1.setType(PETSc.KSP.Type.MINRES)
pc2_1 = solver2_1.getPC()
pc2_1.setType(PETSc.PC.Type.HYPRE)
pc2_1.setHYPREType("boomeramg")

# Solver for step 3
solver3 = PETSc.KSP().create(mesh.comm)
solver3.setOperators(Ac)
solver3.setType(PETSc.KSP.Type.CG)
pc3 = solver3.getPC()
pc3.setType(PETSc.PC.Type.SOR)

solver3_1 = PETSc.KSP().create(mesh.comm)
solver3_1.setOperators(Acc)
solver3_1.setType(PETSc.KSP.Type.CG)
pc3_1 = solver3_1.getPC()
pc3_1.setType(PETSc.PC.Type.SOR)


# 틀린게 많음....






