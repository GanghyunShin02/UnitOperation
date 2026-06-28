#HeatExchager2.py

import numpy as np
import matplotlib.pyplot as plt
import dolfinx
from dolfinx import mesh,io,fem
from dolfinx.fem import ( Function,
                         FunctionSpace,
                         Constant,
                         form,
                         functionspace
                         )
from dolfinx.io import VTXWriter,gmsh as gmshio
from mpi4py import MPI
from petsc4py import PETSc
import ufl
from ufl import ( grad,
                 nabla_grad,
                 dot,
                 inner,
                 TrialFunction,
                 TestFunction,
                 dx,
                 rhs,
                 lhs,
                 div
                )
from basix.ufl import element
import gmsh
from CoolProp.CoolProp import PropsSI
from dolfinx.fem.petsc import (
    assemble_matrix,
    assemble_vector,
    create_vector,
    apply_lifting,
    set_bc,
    create_matrix
)


'''1. domain 생성.. 
2차원 열교환기 이므로 일단 직사각형을 6개 생성.'''

gmsh.initialize()
gdim=2

if MPI.COMM_WORLD.rank==0:
    wall_1=gmsh.model.occ.addRectangle(0,0,0,10,0.1)
    fluid_1=gmsh.model.occ.addRectangle(0,0.1,0,10,1)
    wall_2=gmsh.model.occ.addRectangle(0,1.1,0,10,0.1)
    fluid_2=gmsh.model.occ.addRectangle(0,1.2,0,10,2)
    wall_3=gmsh.model.occ.addRectangle(0,3.2,0,10,0.1)
    fluid_3=gmsh.model.occ.addRectangle(0,3.3,0,10,1)
    wall_4=gmsh.model.occ.addRectangle(0,4.3,0,10,0.1)
    gmsh.model.occ.synchronize()

    # 내관,외관, 외관유체를 같은 피지컬그룹으로 묶기..
    gmsh.model.addPhysicalGroup(gdim,[fluid_2],1,'InnderFluid')
    gmsh.model.addPhysicalGroup(gdim,[wall_2,wall_3],2,'InnderTube')
    gmsh.model.addPhysicalGroup(gdim,[fluid_1,fluid_3],3,'OuterFluid')
    gmsh.model.addPhysicalGroup(gdim,[wall_1,wall_4],4,'OuterTube')

    # 경계조건을위해 벽들만 묶기..
    gmsh.model.addPhysicalGroup(gdim-1,[3],10,'line1')
    gmsh.model.addPhysicalGroup(gdim-1,[9],11,'line2')
    gmsh.model.addPhysicalGroup(gdim-1,[11],12,'line3')
    gmsh.model.addPhysicalGroup(gdim-1,[17],13,'line4')
    gmsh.model.addPhysicalGroup(gdim-1,[19],14,'line5')
    gmsh.model.addPhysicalGroup(gdim-1,[25],15,'line6')

    # 초기조건을위해 inlet,outlet묶기
    gmsh.model.addPhysicalGroup(gdim-1,[8],20,'fluid1_inlet')
    gmsh.model.addPhysicalGroup(gdim-1,[6],21,'fluid1_outlet')
    gmsh.model.addPhysicalGroup(gdim-1,[16],22,'fluid2_inlet')
    gmsh.model.addPhysicalGroup(gdim-1,[14],23,'fluid1_outlet')
    gmsh.model.addPhysicalGroup(gdim-1,[24],24,'fluid3_inlet')
    gmsh.model.addPhysicalGroup(gdim-1,[22],25,'fluid3_outlet')

    # 모서리도 피지컬그룹만들기..
    '''
    boundary_wall1=gmsh.model.getBoundary([(2,wall_1)],oriented=False)
    boundary_wall2=gmsh.model.getBoundary([(2,wall_2)],oriented=False)
    boundary_wall3=gmsh.model.getBoundary([(2,wall_3)],oriented=False)
    boundary_wall4=gmsh.model.getBoundary([(2,wall_4)],oriented=False)
    boundary_fluid1=gmsh.model.getBoundary([(2,fluid_1)],oriented=False)
    boundary_fluid2=gmsh.model.getBoundary([(2,fluid_2)],oriented=False)
    boundary_fluid3=gmsh.model.getBoundary([(2,fluid_3)],oriented=False)
    
    print(f'wall1{boundary_wall1}')
    print(f'fluid1{boundary_fluid1}')
    print(f'wall2{boundary_wall2}')
    print(f'fluid2{boundary_fluid2}')
    print(f'wall3{boundary_wall3}')
    print(f'fluid3{boundary_fluid3}')
    print(f'wall4{boundary_wall4}')
    
    wall1[(1, 1), (1, 2), (1, 3), (1, 4)]
fluid1[(1, 5), (1, 6), (1, 7), (1, 8)]
wall2[(1, 9), (1, 10), (1, 11), (1, 12)]
fluid2[(1, 13), (1, 14), (1, 15), (1, 16)]
wall3[(1, 17), (1, 18), (1, 19), (1, 20)]
fluid3[(1, 21), (1, 22), (1, 23), (1, 24)]
wall4[(1, 25), (1, 26), (1, 27), (1, 28)]





    def classify_line(tag, tol=1e-8):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(1, tag)
        dx = abs(xmax - xmin)
        dy = abs(ymax - ymin)
        if dx < tol:
            return "vertical", (xmin, ymin, xmax, ymax)
        elif dy < tol:
            return "horizontal", (xmin, ymin, xmax, ymax)
        else:
            return "diagonal/unknown", (xmin, ymin, xmax, ymax)

    boundary_dict = {
        "wall1": [1, 2, 3, 4],
        "fluid1": [5, 6, 7, 8],
        "wall2": [9, 10, 11, 12],
        "fluid2": [13, 14, 15, 16],
        "wall3": [17, 18, 19, 20],
        "fluid3": [21, 22, 23, 24],
        "wall4": [25, 26, 27, 28],
    }

    for name, tags in boundary_dict.items():
        print(f"--- {name} ---")
        for t in tags:
            orientation, coords = classify_line(t)
            print(f"  tag {t}: {orientation}, bbox={coords}")


--- wall1 ---
  tag 1: diagonal/unknown, bbox=(-1e-07, -1e-07, 10.0000001, 1e-07)
  tag 2: diagonal/unknown, bbox=(9.9999999, -1e-07, 10.0000001, 0.10000010000000001)
  tag 3: diagonal/unknown, bbox=(-1e-07, 0.0999999, 10.0000001, 0.10000010000000001)
  tag 4: diagonal/unknown, bbox=(-1e-07, -1e-07, 1e-07, 0.10000010000000001)
--- fluid1 ---
  tag 5: diagonal/unknown, bbox=(-1e-07, 0.0999999, 10.0000001, 0.10000010000000001)
  tag 6: diagonal/unknown, bbox=(9.9999999, 0.0999999, 10.0000001, 1.1000001000000001)
  tag 7: diagonal/unknown, bbox=(-1e-07, 1.0999999, 10.0000001, 1.1000001000000001)
  tag 8: diagonal/unknown, bbox=(-1e-07, 0.09999990000000009, 1e-07, 1.1000001000000001)
--- wall2 ---
  tag 9: diagonal/unknown, bbox=(-1e-07, 1.0999999, 10.0000001, 1.1000001000000001)
  tag 10: diagonal/unknown, bbox=(9.9999999, 1.0999999, 10.0000001, 1.2000001000000002)
  tag 11: diagonal/unknown, bbox=(-1e-07, 1.1999999000000001, 10.0000001, 1.2000001000000002)
  tag 12: diagonal/unknown, bbox=(-1e-07, 1.0999999, 1e-07, 1.2000001000000002)
--- fluid2 ---
  tag 13: diagonal/unknown, bbox=(-1e-07, 1.1999999, 10.0000001, 1.2000001)
  tag 14: diagonal/unknown, bbox=(9.9999999, 1.1999999, 10.0000001, 3.2000001)
  tag 15: diagonal/unknown, bbox=(-1e-07, 3.1999999000000003, 10.0000001, 3.2000001)
  tag 16: diagonal/unknown, bbox=(-1e-07, 1.1999999000000001, 1e-07, 3.2000001)
--- wall3 ---
  tag 17: diagonal/unknown, bbox=(-1e-07, 3.1999999000000003, 10.0000001, 3.2000001)
  tag 18: diagonal/unknown, bbox=(9.9999999, 3.1999999000000003, 10.0000001, 3.3000001)
  tag 19: diagonal/unknown, bbox=(-1e-07, 3.2999999000000004, 10.0000001, 3.3000001)
  tag 20: diagonal/unknown, bbox=(-1e-07, 3.1999999000000003, 1e-07, 3.3000001)
--- fluid3 ---
  tag 21: diagonal/unknown, bbox=(-1e-07, 3.2999999, 10.0000001, 3.3000000999999997)
  tag 22: diagonal/unknown, bbox=(9.9999999, 3.2999999, 10.0000001, 4.3000001)
  tag 23: diagonal/unknown, bbox=(-1e-07, 4.2999998999999995, 10.0000001, 4.3000001)
  tag 24: diagonal/unknown, bbox=(-1e-07, 3.2999999, 1e-07, 4.3000001)
--- wall4 ---
  tag 25: diagonal/unknown, bbox=(-1e-07, 4.2999998999999995, 10.0000001, 4.3000001)
  tag 26: diagonal/unknown, bbox=(9.9999999, 4.2999998999999995, 10.0000001, 4.4000001)
  tag 27: diagonal/unknown, bbox=(-1e-07, 4.399999899999999, 10.0000001, 4.4000001)
  tag 28: diagonal/unknown, bbox=(-1e-07, 4.2999998999999995, 1e-07, 4.4000001)
  '''





gmsh.option.setNumber('Mesh.MeshSizeMax',0.1)
gmsh.option.setNumber('Mesh.MeshSizeMin',0.05)
gmsh.model.mesh.generate(gdim)


mesh_data=gmshio.model_to_mesh(gmsh.model,MPI.COMM_WORLD,0,gdim)
domain=mesh_data.mesh

facet_marker=mesh_data.facet_tags
cell_marker=mesh_data.cell_tags

gmsh.finalize()

'''
with io.XDMFFile(domain.comm, "heat_exchanger2_geometry2.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    domain.topology.create_connectivity(domain.topology.dim, domain.topology.dim)
    xdmf.write_meshtags(cell_marker, domain.geometry)
'''

# 유동장해석

def rho(T):
   return PropsSI('D','T',T,'P',101325,'water')
def mu(T):
  return PropsSI('V','T',T,'P',101325,'water')
def nu(T):
   return mu(T)/rho(T)

rho_test=Constant(domain,PETSc.ScalarType(rho(300)))
mu_test=Constant(domain,PETSc.ScalarType(mu(300)))
nu_test=Constant(domain,PETSc.ScalarType(nu(300)))

# 속도와 압력 함수공간
# element로 요소를 만들어놓으면 서로 섞어서 쓸수있음..f(v,p) 처럼..
v_cg1=element('Lagrange',domain.basix_cell(),2,shape=(domain.geometry.dim,))
v_cg2=element('Lagrange',domain.basix_cell(),2,shape=(domain.geometry.dim,))
#v_cg3=element('Lagrange',domain.basix_cell(),2,shape=(domain.geometry.dim))
q_cg1=element('Lagrange',domain.basix_cell(),1)
q_cg2=element('Lagrange',domain.basix_cell(),1)
#q_cg3=element('Lagrange',domain.basix_cell(),1)

#함수는 mesh전체에서 선언하고 적분만 따로.

v_1=functionspace(domain,v_cg1)
v_2=functionspace(domain,v_cg2)

P_1=functionspace(domain,q_cg1)
P_2=functionspace(domain,q_cg2)

# 경계조건.. 벽에서 속도가 0
# 자유도를 찾아야..

u_noslip=np.zeros(gdim,dtype=PETSc.ScalarType)

print('ff')
inner_facettag=facet_marker.find(2)
outer_facettag=facet_marker.find(4)
print('facet')
fdim=gdim-1


inner_dofs=fem.locate_dofs_topological(v_1,fdim,inner_facettag)
outer_dofs=fem.locate_dofs_topological(v_2,fdim,outer_facettag)
print('1')

bc_inner=fem.dirichletbc(u_noslip,inner_dofs,v_1)
bc_outter=fem.dirichletbc(u_noslip,outer_dofs,v_2)


u_in,v_in=TrialFunction(v_1),TestFunction(v_1)
u_n_in=Function(v_1)
u_n1_in=Function(v_1)
u_s_in=Function(v_1)
P_in,q_in=TrialFunction(P_1),TestFunction(P_1)
P__i,phi_in=Function(P_1),Function(P_1)

u_out,v_out=TrialFunction(v_2),TestFunction(v_2)
u_n_out=Function(v_2)
u_n1_out=Function(v_2)
u_s_out=Function(v_2)
P_out,q_out=TrialFunction(P_2),TestFunction(P_2)
P__O,phi_out=Function(P_2),Function(P_2)


dt=Constant(domain,PETSc.ScalarType(0.1))
dx=ufl.Measure("dx",domain=domain,subdomain_data=cell_marker)

print('2')

#1.속도추정치
F_in=(rho_test/dt)*dot((u_in-u_n_in),v_in)*dx(2)
F_in+=rho_test*inner(dot((1.5*u_n_in-0.5*u_n1_in),nabla_grad(0.5*(u_in+u_n_in))),v_in)*dx(2)
F_in+=0.5*mu_test*inner(nabla_grad(v_in),nabla_grad(u_in+u_n_in))*dx(2)
F_in-=P__i*div(v_in)*dx(2)

a_in=form(lhs(F_in))
L_in=form(rhs(F_in))
A_in=create_matrix(a_in)
b_in=create_vector(fem.extract_function_spaces(L_in))

F_out=(rho_test/dt)*dot((u_out-u_n_out),v_out)*dx(4)
F_out+=rho_test*inner(dot((1.5*u_n_out-0.5*u_n1_out),nabla_grad(0.5*(u_out+u_n_out))),v_out)*dx(4)
F_out+=0.5*mu_test*inner(nabla_grad(v_out),nabla_grad(u_out+u_n_out))*dx(4)
F_out-=P__O*div(v_out)*dx(4)

a_out=form(lhs(F_out))
L_out=form(rhs(F_out))
A_out=create_matrix(a_out)
b_out=create_vector(fem.extract_function_spaces(L_out))

#2 압력보정
a_2_in=form(dot(grad(q_in),grad(P_in))*dx(2))
L_2_in=form((rho_test/dt)*dot(q_in,div(u_s_in))*dx(2))
A2_in=assemble_matrix(a_2_in)
A2_in.assemble()
b2_in=create_vector(fem.extract_function_spaces(L_2_in))

a_2_out=form(dot(grad(q_out),grad(P_out))*dx(4))
L_2_out=form((rho_test/dt)*dot(q_out,div(u_s_out))*dx(4))
A2_out=assemble_matrix(a_2_out)
A2_out.assemble()
b2_out=create_vector(fem.extract_function_spaces(L_2_out))

#3. 속도보정
a_3_in=form(rho_test*dot(u_in,v_in)*dx(2))
L_3_in=form(rho_test*dot(u_n1_in,v_in)*dx(2)-dt*dot(v_in,nabla_grad(phi_in))*dx(2))
A3_in=assemble_matrix(a_3_in)
A3_in.assemble()
b3_in=create_vector(fem.extract_function_spaces(L_3_in))

a_3_out=form(rho_test*dot(u_out,v_out)*dx(4))
L_3_out=form(rho_test*dot(u_n1_out,v_out)*dx(4)-dt*dot(v_out,nabla_grad(phi_out))*dx(4))
A3_out=assemble_matrix(a_3_out)
A3_out.assemble()
b3_out=create_vector(fem.extract_function_spaces(L_3_out))

print('3')

# Solver for step 1
solver1_in = PETSc.KSP().create(domain.comm)
solver1_in.setOperators(A_in)
solver1_in.setType(PETSc.KSP.Type.BCGS)
pc1_in = solver1_in.getPC()
pc1_in.setType(PETSc.PC.Type.JACOBI)

solver1_out = PETSc.KSP().create(domain.comm)
solver1_out.setOperators(A_out)
solver1_out.setType(PETSc.KSP.Type.BCGS)
pc1_out = solver1_out.getPC()
pc1_out.setType(PETSc.PC.Type.JACOBI)

# Solver for step 2
solver2_in = PETSc.KSP().create(domain.comm)
solver2_in.setOperators(A2_in)
solver2_in.setType(PETSc.KSP.Type.MINRES)
pc2_in = solver2_in.getPC()
pc2_in.setType(PETSc.PC.Type.HYPRE)
pc2_in.setHYPREType("boomeramg")

solver2_out = PETSc.KSP().create(domain.comm)
solver2_out.setOperators(A2_out)
solver2_out.setType(PETSc.KSP.Type.MINRES)
pc2_out = solver2_out.getPC()
pc2_out.setType(PETSc.PC.Type.HYPRE)
pc2_out.setHYPREType("boomeramg")


# Solver for step 3
solver3_in = PETSc.KSP().create(domain.comm)
solver3_in.setOperators(A3_in)
solver3_in.setType(PETSc.KSP.Type.CG)
pc3_in = solver3_in.getPC()
pc3_in.setType(PETSc.PC.Type.SOR)

solver3_out = PETSc.KSP().create(domain.comm)
solver3_out.setOperators(A3_out)
solver3_out.setType(PETSc.KSP.Type.CG)
pc3_out = solver3_out.getPC()
pc3_out.setType(PETSc.PC.Type.SOR)

print('4')
from pathlib import Path

folder = Path("results")
folder.mkdir(exist_ok=True)

vtx_u_in = VTXWriter(domain.comm, folder/"u_in.bp", [u_n_in], engine="BP4")
vtx_p_in = VTXWriter(domain.comm, folder/"p_in.bp", [P__i], engine="BP4")

vtx_u_out = VTXWriter(domain.comm, folder/"u_out.bp", [u_n_out], engine="BP4")
vtx_p_out = VTXWriter(domain.comm, folder/"p_out.bp", [P__O], engine="BP4")

t = 0.0
dt_value = float(dt.value)
num_steps = 100

for n in range(num_steps):
    print('A')
    t += dt_value

    # ============================
    # Step 1 : Tentative velocity
    # ============================

    A_in.zeroEntries()
    assemble_matrix(A_in, a_in, bcs=[bc_inner])
    A_in.assemble()
    print('B')
    with b_in.localForm() as loc:
        loc.set(0)

    assemble_vector(b_in, L_in)
    apply_lifting(b_in, [a_in], [[bc_inner]])
    b_in.ghostUpdate(
        addv=PETSc.InsertMode.ADD_VALUES,
        mode=PETSc.ScatterMode.REVERSE
    )
    set_bc(b_in, [bc_inner])

    solver1_in.solve(b_in, u_s_in.x.petsc_vec)
    u_s_in.x.scatter_forward()

    print('C')
    A_out.zeroEntries()
    assemble_matrix(A_out, a_out, bcs=[bc_outter])
    A_out.assemble()

    with b_out.localForm() as loc:
        loc.set(0)

    print('D')
    assemble_vector(b_out, L_out)
    apply_lifting(b_out, [a_out], [[bc_outter]])
    b_out.ghostUpdate(
        addv=PETSc.InsertMode.ADD_VALUES,
        mode=PETSc.ScatterMode.REVERSE
    )
    set_bc(b_out, [bc_outter])

    solver1_out.solve(b_out, u_s_out.x.petsc_vec)
    u_s_out.x.scatter_forward()


    # ============================
    # Step 2 : Pressure correction
    # ============================
    print('E')
    with b2_in.localForm() as loc:
        loc.set(0)

    assemble_vector(b2_in, L_2_in)

    solver2_in.solve(b2_in, phi_in.x.petsc_vec)
    phi_in.x.scatter_forward()

    P__i.x.petsc_vec.axpy(1.0, phi_in.x.petsc_vec)
    P__i.x.scatter_forward()
    print('F')

    with b2_out.localForm() as loc:
        loc.set(0)

    assemble_vector(b2_out, L_2_out)

    solver2_out.solve(b2_out, phi_out.x.petsc_vec)
    phi_out.x.scatter_forward()

    P__O.x.petsc_vec.axpy(1.0, phi_out.x.petsc_vec)
    P__O.x.scatter_forward()

    print('G')
    # ============================
    # Step 3 : Velocity correction
    # ============================

    with b3_in.localForm() as loc:
        loc.set(0)

    assemble_vector(b3_in, L_3_in)

    solver3_in.solve(b3_in, u_n_in.x.petsc_vec)
    u_n_in.x.scatter_forward()

    print('H')
    with b3_out.localForm() as loc:
        loc.set(0)

    assemble_vector(b3_out, L_3_out)

    solver3_out.solve(b3_out, u_n_out.x.petsc_vec)
    u_n_out.x.scatter_forward()

    print('I')
    # ============================
    # Write results
    # ============================

    vtx_u_in.write(t)
    vtx_p_in.write(t)

    vtx_u_out.write(t)
    vtx_p_out.write(t)


    # ============================
    # Update previous solution
    # ============================

    u_n1_in.x.array[:] = u_n_in.x.array
    u_n1_out.x.array[:] = u_n_out.x.array

vtx_u_in.close()
vtx_p_in.close()
vtx_u_out.close()
vtx_p_out.close()



