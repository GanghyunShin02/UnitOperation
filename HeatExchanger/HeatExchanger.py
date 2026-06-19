import numpy as np
import matplotlib.pyplot as plt
import dolfinx
from mpi4py import MPI
import gmsh
import ufl
from ufl import grad,TestFunction,TrialFunction,dx
from dolfinx.fem import Constant, Function,functionspace
from dolfinx import mesh,fem,io
from dolfinx.io import VTXWriter,gmsh as gmshio
from pathlib import Path

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
    gmsh.model.addPhysicalGroup(3,[wall_in_tag],2,'wall_in')
    gmsh.model.addPhysicalGroup(3,[fluid_out_tag],3,'fluid_out')
    gmsh.model.addPhysicalGroup(3,[wall_out_tag],4,'wall_out')

    gmsh.option.setNumber('Mesh.Algorithm3D',1)
    gmsh.option.setNumber('Mesh.MeshSizeMin',0.05)
    gmsh.option.setNumber('Mesh.MeshSizeMax',0.1)
    gmsh.option.setNumber('Geometry.Tolerance',1e-6)
    gmsh.model.mesh.generate(gdim)

#print(f'modeltomesh{gmshio.model_to_mesh(gmsh.model,MPI.COMM_WORLD,0,gdim)}')
meshdata=gmshio.model_to_mesh(gmsh.model,MPI.COMM_WORLD,0,gdim)
domain=meshdata.mesh
facet_marker=meshdata.facet_tags
cell_marker=meshdata.cell_tags
gmsh.finalize()

with io.XDMFFile(domain.comm, "heat_exchanger_geometry.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    domain.topology.create_connectivity(domain.topology.dim, domain.topology.dim)
    xdmf.write_meshtags(cell_marker, domain.geometry)













    