"""
Christian Austin, University of Florida
Part of Research for PhD Thesis, Summer of 2024

Implements the elasto-viscous split stress (EVSS) formulation 
of the Oldroyd 3-parameter subset model. The model itself is 
described by Scott and Girault in 2021, and can be viewed as 
a special case of the Oldroyd 8-Parameter Model (1958). 

The O3 model lacks explicit dissipation when the mixed
formulation is naively implemented. Scott and Girault 
introduced their SRTD formulation to make the momentum 
equation explicitly elliptic, and then also gave an iterative
algorithm for solving the system. 

The EVSS formulation was devised by Rajagopalan, Armstrong,
and Brown in 1990 for constitutive equations which lack 
explicit dissipation. Their original paper only did it for the
UCM and Giesekus model, but they mentioned that this technique
can be modified to many different models. It can, in fact, be 
applied to the O3 Model, and we do exactly that here. 
"""

from fenics import *
#from mshr import *
from meshdata import gen_mesh_jb
import sys
import os

class Results:
    def __init__(self, velocity, deformation, pressure, Sigma, stress_tensor):
        self.velocity = velocity
        self.deformation = deformation
        self.pressure = pressure
        self.Sigma = Sigma
        self.stress_tensor = stress_tensor

# Deformation d is treated as a separate variable to account for the extra derivative introduced

# Journal Bearing Problem 

def oldroyd_3_JB_EVSS(h, rad, ecc, eta, l1, mu1):
    if(rad>=1 or rad<=0 or ecc<0 or rad+ecc>1):
        #throw exception, forgot how lol
        print("Error: Inputs not valid")
    
    meshfile = "meshdata/journal_bearing_h_%.4e.h5"%h
    
    if not os.path.exists(meshfile):
        print("Creating mesh...")
        gen_mesh_jb.main(h, rad, ecc)
    
    #then, simply read the mesh in 
    mesh = Mesh() #empty mesh
    infile = HDF5File(MPI.comm_world, meshfile, 'r')
    infile.read(mesh, '/mesh', True) #for some reason, need this flag to import a mesh?
    infile.close()
    print("Mesh loaded into FEniCS")
    
    
    #boundaries of domain
    class Inner(SubDomain):
        def inside(self, x, on_boundary):
            radius = x[0]*x[0] + (x[1]+ecc)*(x[1]+ecc)
            return (on_boundary and radius <= rad*rad+h)
    
    class Outer(SubDomain):
        def inside(self, x, on_boundary):
            # on_boundary and opposite of inner lol
            radius = x[0]*x[0] + (x[1]+ecc)*(x[1]+ecc)
            return (on_boundary and radius > rad*rad+h)
    
    class TopPoint(SubDomain):
        def inside(self, x, on_boundary):
            return (near(x[0], 0.0) and near(x[1], 1.0))
    
    # Boundary data     
    speed_outer = 0.0 # counter-clockwise tangential speed of outer bearing, 
    speed_inner = 1.0 # clockwise tangential speed of inner bearing, 
    g_inner = Expression(("s*(x[1]+ecc)/r", "-s*x[0]/r"), s=speed_inner, r=rad, ecc=ecc, degree=1) 
    g_outer = Expression(("-s*x[1]", "s*x[0]"), s=speed_outer, degree=1) # no slip on outer bearing 
    f = Constant((0.0, 0.0)) # no body forces
    
    # Element spaces
    V_elem = VectorElement("CG", triangle, 2) # Velocity, degree 2 elements
    P_elem = FiniteElement("CG", triangle, 1) # Pressure, degree 1 elements
    T_elem = TensorElement("CG", triangle, 2, symmetry=True) # Stress tensor, degree 2 elements
    D_elem = TensorElement("CG", triangle, 1, symmetry=True) # Deformation tensor, linear
    
    W_elem = MixedElement([V_elem, P_elem, T_elem, D_elem]) # Mixed element (u, p, T, D)
    
    # Function Spaces    
    W = FunctionSpace(mesh, W_elem) # Only need one, the mixed Function Space
    
    # Interpolate body force and BCs onto velocity FE space
    g_inner = interpolate(g_inner, W.sub(0).collapse())
    g_outer = interpolate(g_outer, W.sub(0).collapse())
    f = interpolate(f, W.sub(0).collapse())

    # Define boundary conditions
    bc_inner = DirichletBC(W.sub(0), g_inner, Inner())
    bc_outer = DirichletBC(W.sub(0), g_outer, Outer())
    bc_press = DirichletBC(W.sub(1), Constant(0.0), TopPoint(), 'pointwise')
    
    # Gather boundary conditions (any others would go here, separated by a comma)
    bcs = [bc_inner, bc_outer, bc_press] 
    
    # Variational Problem: Trial and Test Functions
    w = TrialFunction(W) 
    (u, p, Tau, D) = split(w) #trial/solution functions
    (v, q, S, Phi) = TestFunctions(W) #test functions
    
    # Momentum equation gets velocity TFs v
    # EVSS SF: -eta*Lapl(u) + del(u)*u + del(p) - del.Sigma - f = 0
    momentum = eta*inner(grad(u), grad(v))*dx + inner(dot(grad(u), u) + grad(p) - div(Tau) - f, v)*dx
    
    # Continuity equation gets pressure TFs q
    # SF: del.u = 0
    continuity = div(u)*q*dx
    
    # Constitutive equation gets degree 2 tensor TFs and SUPG (S + h*del(S)*u)
    # CE SF: T + l1*G(T + 2etaE, a) = 0
    h = mesh.hmax()/4
    constitutive = inner( Tau + l1*(dot(grad(Tau+2*eta*D), u) - dot(grad(u), Tau+2*eta*D) - dot(Tau+2*eta*D, grad(u).T)) \
                         +(l1-mu1)*(dot( (grad(u)+grad(u).T)/2, Tau+2*eta*D) + dot(Tau+2*eta*D,(1/2)*(grad(u)+grad(u).T))),\
                         S + h*dot(grad(S), u) )*dx
    
    # Enforcement of deformation tensor gets gets degree 1 tensor TFs Phi
    enforcement = inner(D - (1/2)*(grad(u) + grad(u).T), Phi)*dx
    
    
    # Sum all together
    F = momentum + continuity + constitutive + enforcement
    
    # Solve discrete variational form
    you = Function(W)
    
    F_act = action(F, you) 
    dF = derivative(F_act, you)

    problem = NonlinearVariationalProblem(F_act, you, bcs, dF)
    solver = NonlinearVariationalSolver(problem)
    # need to adjust solver parameters b/c PetSC runs out of memory if default solver is used. 
    prm = solver.parameters
    prm["nonlinear_solver"] = "newton"
    prm["newton_solver"]["linear_solver"] = "mumps" # utilizes parallel processors
    solver.solve()
    
    u1, p1, Sigma, D1 = you.split(deepcopy = True)
    
    #remember, Sigma here is not the stress tensor, but the "modified" tensor Sigma. Return both
    stress = project( Sigma + 2*eta*D1, W.sub(3).collapse()) #projected onto lower degree space
    
    return Results(u1, D1, p1, Sigma, stress)
    
    
# Lid-Driven Cavity Problem
def oldroyd_3_LDC_EVSS(meshsize, corner_ref, total_ref, eta, l1, mu1):
    # boundary data
    g_top = Expression(("30.0*x[0]*x[0]*(1-x[0])*(1-x[0])", "0.0"), degree = 4) # 30x^2(1-x)^2, 30 gives it integral=1
    g_walls = Constant((0.0, 0.0)) #g=0 on walls

    # body forces
    f = Constant((0.0, 0.0)) # no body forces
    
    # Variational problem start: define domain and mesh
    nx = meshsize; ny = meshsize
    mesh = UnitSquareMesh(nx,ny, "crossed")

    # Refine mesh near top corners, where singularities develop for this problem
    for n in range(corner_ref): 
        tl_markers = MeshFunction('bool', mesh, mesh.topology().dim(), False)
        for cell in cells(mesh):
            tl_markers[cell] = cell.contains(Point(0.0, 1.0))
        mesh = refine(mesh, tl_markers)

        tr_markers = MeshFunction('bool', mesh, mesh.topology().dim(), False)
        for cell in cells(mesh):
            tr_markers[cell] = cell.contains(Point(1.0, 1.0))
        mesh = refine(mesh, tr_markers)
    
    for m in range(total_ref):
        mesh = refine(mesh)
    
    # charcteristic mesh size for streamline upwinding
    h = mesh.hmax()/4.0 # for no SUPG, h=0. Not sure how to make this optional yet. 
    
    # Element spaces
    V_elem = VectorElement("CG", triangle, 2) # Velocity, degree 2 elements
    P_elem = FiniteElement("CG", triangle, 1) # Pressure, degree 1 elements
    T_elem = TensorElement("CG", triangle, 2, symmetry=True) # Stress tensor, degree 2 elements
    D_elem = TensorElement("CG", triangle, 1, symmetry=True) # Deformation tensor, linear
    
    W_elem = MixedElement([V_elem, P_elem, T_elem, D_elem]) # Mixed element (u, p, T, D)
    
    # Function Spaces    
    W = FunctionSpace(mesh, W_elem) # Only need one, the mixed Function Space
    
    # Interpolate body force and BCs onto velocity FE space
    g_top = interpolate(g_top, W.sub(0).collapse())
    g_walls = interpolate(g_walls, W.sub(0).collapse())
    f = interpolate(f, W.sub(0).collapse())
    
    # Define boundary conditions
    top_lid = 'near(x[1], 1.0) && on_boundary'
    walls = '(near(x[1], 0.0) || near(x[0], 0.0) || near(x[0], 1.0)) && on_boundary'
    bl_corner = 'near(x[0], 0.0) && near(x[1], 0.0)' #for pressure regulating

    bc_top = DirichletBC(W.sub(0), g_top, top_lid)
    bc_walls = DirichletBC(W.sub(0), g_walls, walls)
    pressure_reg = DirichletBC(W.sub(1), Constant(0.0), bl_corner, 'pointwise')
        
    # Gather boundary conditions (any others would go here, separated by a comma)
    bcs = [bc_top, bc_walls, pressure_reg] 
    
    # Variational Problem: Trial and Test Functions
    w = TrialFunction(W) 
    (u, p, Tau, D) = split(w) #trial/solution functions
    (v, q, S, Phi) = TestFunctions(W) #test functions
    
    
    # Momentum equation gets velocity TFs v
    # EVSS SF: -eta*Lapl(u) + del(u)*u + del(p) - del.Sigma - f = 0
    momentum = eta*inner(grad(u), grad(v))*dx + inner(dot(grad(u), u) + grad(p) - div(Tau) - f, v)*dx
    
    # Continuity equation gets pressure TFs q
    # SF: del.u = 0
    continuity = div(u)*q*dx
    
    # Constitutive equation gets degree 2 tensor TFs and SUPG (S + h*del(S)*u)
    # CE SF: T + l1*G(T + 2etaE, a) = 0
    h = mesh.hmax()/4
    constitutive = inner( Tau + l1*(dot(grad(Tau+2*eta*D), u) - dot(grad(u), Tau+2*eta*D) - dot(Tau+2*eta*D, grad(u).T)) \
                         +(l1-mu1)*(dot( (grad(u)+grad(u).T)/2, Tau+2*eta*D) + dot(Tau+2*eta*D,(1/2)*(grad(u)+grad(u).T))),\
                         S + h*dot(grad(S), u) )*dx
    
    # Enforcement of deformation tensor gets gets degree 1 tensor TFs Phi
    enforcement = inner(D - (1/2)*(grad(u) + grad(u).T), Phi)*dx
    
    # Sum all together
    F = momentum + continuity + constitutive + enforcement
    
    # Solve discrete variational form
    you = Function(W)
    
    F_act = action(F, you) 
    dF = derivative(F_act, you)

    problem = NonlinearVariationalProblem(F_act, you, bcs, dF)
    solver = NonlinearVariationalSolver(problem)
    # need to adjust solver parameters b/c PetSC runs out of memory if default solver is used. 
    prm = solver.parameters
    prm["nonlinear_solver"] = "newton"
    prm["newton_solver"]["linear_solver"] = "mumps" # utilizes parallel processors
    solver.solve()
    
    u1, p1, Sigma, D1 = you.split(deepcopy = True)
    
    #remember, Sigma here is not the stress tensor, but the "modified" tensor Sigma. Return both
    stress = project( Sigma + 2*eta*D1, W.sub(3).collapse()) #projected onto lower degree space
    
    return Results(u1, D1, p1, Sigma, stress)

#post proc stuff here

