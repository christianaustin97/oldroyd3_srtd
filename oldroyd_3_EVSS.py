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
from meshdata import gen_mesh_jb
from meshdata import gen_mesh_ldc
import steady_nse_solver
import os
import matplotlib.pyplot as plt

class Results:
    def __init__(self, converged, velocity, deformation, pressure, Sigma, num_iters):
        self.converged = converged
        self.velocity = velocity
        self.deformation = deformation
        self.pressure = pressure
        self.Sigma = Sigma
        #self.stress_tensor = stress_tensor
        self.num_iters = num_iters

# Deformation d is treated as a separate variable to account for the extra derivative introduced

# Journal Bearing Problem 

def oldroyd_3_JB_EVSS(h, rad, ecc, s, eta, l1, mu1):
    # s is the tangential speed of the bearing 

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
    speed_inner = s # clockwise tangential speed of inner bearing, 
    g_inner = Expression(("s*(x[1]+ecc)/r", "-s*x[0]/r"), s=speed_inner, r=rad, ecc=ecc, degree=1) 
    g_outer = Expression(("-s*x[1]", "s*x[0]"), s=speed_outer, degree=1) # no slip on outer bearing 
    f = Constant((0.0, 0.0)) # no body forces
    
    # Element spaces
    V_elem = VectorElement("CG", triangle, 2) # Velocity, degree 2 elements
    P_elem = FiniteElement("CG", triangle, 1) # Pressure, degree 1 elements
    T_elem = TensorElement("CG", triangle, 2, symmetry=True) # "Stress" tensor (actually Sigma), degree 2 elements
    D_elem = VectorElement("CG", triangle, 1) # Deformation tensor, defined as VectorElement to exploit symmetry

    W_elem = MixedElement([V_elem, P_elem, T_elem, D_elem]) # Mixed element (u, p, T, D)
    
    # Function Spaces
    V_space = FunctionSpace(mesh, V_elem) # need these first 4 for FunctionAssigner
    P_space = FunctionSpace(mesh, P_elem)
    T_space = FunctionSpace(mesh, T_elem)
    D_space = FunctionSpace(mesh, D_elem)
    W = FunctionSpace(mesh, W_elem) # Only need one really, the mixed Function Space
    
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
    (u, p, Tau, D_vec) = split(w) #trial/solution functions
    (v, q, S, Phi_vec) = TestFunctions(W) #test functions

    D = as_tensor([[D_vec[0], D_vec[1]], [D_vec[1], -D_vec[0]]]) # D is symmetric and traceless
    Phi = as_tensor([[Phi_vec[0], Phi_vec[1]], [Phi_vec[1], -Phi_vec[0]]])
    
    # Momentum equation gets velocity TFs v
    # EVSS SF: -eta*Lapl(u) + del(u)*u + del(p) - del.Sigma - f = 0
    momentum = eta*inner(grad(u), grad(v))*dx + inner(dot(grad(u), u) + grad(p) - div(Tau) - f, v)*dx
    
    # Continuity equation gets pressure TFs q
    # SF: del.u = 0
    continuity = div(u)*q*dx
    
    # Constitutive equation gets degree 2 tensor TFs and SUPG (S + h*del(S)*u)
    # CE SF: T + l1*G(T + 2etaE, a) = 0
    constitutive = inner( Tau + l1*(dot(grad(Tau+2*eta*D), u) - dot(grad(u), Tau+2*eta*D) - dot(Tau+2*eta*D, grad(u).T)) \
                         +(l1-mu1)*(dot( (grad(u)+grad(u).T)/2, Tau+2*eta*D) + dot(Tau+2*eta*D,(1/2)*(grad(u)+grad(u).T))),\
                         S + h*dot(grad(S), u) )*dx
    
    # Enforcement of deformation tensor gets gets degree 1 tensor TFs Phi
    enforcement = inner(D - (1/2)*(grad(u) + grad(u).T), Phi)*dx
    
    
    # Sum all together
    F = momentum + continuity + constitutive + enforcement
    
    # Solve discrete variational form
    you = Function(W) #starting guess for Newton

    #get NSE as starting guess for Newton
    nse_solution = steady_nse_solver.navier_stokes_JB(h, rad, ecc, s, eta) 
    print("NSE Solver done")

    u_nse = interpolate(nse_solution.velocity, V_space)
    p_nse = interpolate(nse_solution.pressure, P_space)

    # Found that just using velocity and pressure works best, for some reason
    D_start = interpolate(Constant((0.0, 0.0)), D_space)
    Sigma_start = interpolate( Constant( ((0.0, 0.0), (0.0, 0.0)) ), T_space)

    assigner = FunctionAssigner(W, [V_space, P_space, T_space, D_space])
    assigner.assign(you, [u_nse, p_nse, Sigma_start, D_start])
    
    # done with NSE solution, which is stored in "you" compliant with the mixed form of you
    
    # Now onto EVSS Solve
    F_act = action(F, you) 
    dF = derivative(F_act, you)

    problem = NonlinearVariationalProblem(F_act, you, bcs, dF)
    solver = NonlinearVariationalSolver(problem)
    # need to adjust solver parameters b/c PetSC runs out of memory if default solver is used. 
    prm = solver.parameters
    prm["nonlinear_solver"] = "newton"
    prm["newton_solver"]["linear_solver"] = "mumps" # utilizes parallel processors
    
    try:
        (iters, converged) = solver.solve()
    except:
        print("Newton Solver failed to converge")
        converged = False
        iters = -1
    
    u1, p1, Sigma, D1_vec = you.split(deepcopy = True)
    D1 = as_tensor([[D1_vec[0], D1_vec[1]], [D1_vec[1], -D1_vec[0]]]) # Reshape the strain/velocity gradient tensor

    #remember, Sigma here is not the stress tensor, but the "modified" tensor Sigma. Return both
    #stress = project(Sigma + 2*eta*D1, FunctionSpace(mesh, TensorElement("CG", triangle, 1, symmetry=True))) # project onto deg 1 space
    
    return Results(converged, u1, D1, p1, Sigma, iters)
    
    
# Lid-Driven Cavity Problem

def oldroyd_3_LDC_EVSS(h, s, eta, l1, mu1):
    # s is the average velocity of the top lid 
    nx = round(1/h)
    mesh = UnitSquareMesh(nx, nx)
    print("Mesh loaded into FEniCS")

    # boundary data
    g_top = Expression(("s*16.0*x[0]*x[0]*(1-x[0])*(1-x[0])", "0.0"), s=s, degree = 4)
    g_walls = Constant((0.0, 0.0)) #g=0 on walls

    # body forces
    f = Constant((0.0, 0.0)) # no body forces
    
    # Element spaces
    V_elem = VectorElement("CG", triangle, 2) # Velocity, degree 2 elements
    P_elem = FiniteElement("CG", triangle, 1) # Pressure, degree 1 elements
    T_elem = TensorElement("CG", triangle, 2, symmetry=True) # elastic "Stress" tensor tau, degree 2 elements
    D_elem = VectorElement("CG", triangle, 1) # Deformation tensor, defined as VectorElement to exploit symmetry

    W_elem = MixedElement([V_elem, P_elem, T_elem, D_elem]) # Mixed element (u, p, Sigma, D_Vec)
    
    # Function Spaces
    V_space = FunctionSpace(mesh, V_elem) # need these first 4 for FunctionAssigner
    P_space = FunctionSpace(mesh, P_elem)
    T_space = FunctionSpace(mesh, T_elem)
    D_space = FunctionSpace(mesh, D_elem)
    W = FunctionSpace(mesh, W_elem) # Only need this one for pure EVSS, the mixed Function Space
    
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
    (u, p, Tau, D_vec) = split(w) #trial/solution functions
    (v, q, S, Phi_vec) = TestFunctions(W) #test functions

    D = as_tensor([[D_vec[0], D_vec[1]], [D_vec[1], -D_vec[0]]]) # D is symmetric and traceless...
    Phi = as_tensor([[Phi_vec[0], Phi_vec[1]], [Phi_vec[1], -Phi_vec[0]]]) # ... and so are test funcs for enforcement  
    
    # Momentum equation gets velocity Test Functs v
    # EVSS SF: -eta*Lapl(u) + del(u)*u + del(p) - del.Sigma - f = 0
    momentum = eta*inner(grad(u), grad(v))*dx + inner(dot(grad(u), u) + grad(p) - div(Tau) - f, v)*dx
    
    # Continuity equation gets pressure TFs q
    # SF: del.u = 0
    continuity = div(u)*q*dx
    
    # Constitutive equation gets degree 2 tensor TFs and SUPG (S + h*del(S)*u)
    # CE SF: T + l1*G(T + 2etaE, a) = 0
    constitutive = inner( Tau + l1*(dot(grad(Tau+2*eta*D), u) - dot(grad(u), Tau+2*eta*D) - dot(Tau+2*eta*D, grad(u).T)) \
                         +(l1-mu1)*(dot( (grad(u)+grad(u).T)/2, Tau+2*eta*D) + dot(Tau+2*eta*D,(1/2)*(grad(u)+grad(u).T))),\
                         S + h*dot(grad(S), u) )*dx
    
    # Enforcement of deformation tensor gets gets degree 1 tensor TFs Phi
    enforcement = inner(D - (1/2)*(grad(u) + grad(u).T), Phi)*dx
    
    # Sum all together
    F = momentum + continuity + constitutive + enforcement
    
    # Solve discrete variational form
    you = Function(W) #starting guess for Newton, default Function() is constant zeros
    
    #get NSE as starting guess for Newton
    nse_solution = steady_nse_solver.navier_stokes_LDC(h, s, eta) 
    print("NSE Solver done")

    u_nse = interpolate(nse_solution.velocity, V_space)
    p_nse = interpolate(nse_solution.pressure, P_space)

    # Found that just using velocity and pressure works best, for some reason
    D_start = interpolate(Constant((0.0, 0.0)), D_space)
    Sigma_start = interpolate( Constant( ((0.0, 0.0), (0.0, 0.0)) ), T_space)
    
    assigner = FunctionAssigner(W, [V_space, P_space, T_space, D_space])
    assigner.assign(you, [u_nse, p_nse, Sigma_start, D_start])

    # done with NSE solution, which is stored in "you" compliant with the mixed form of you
    
    F_act = action(F, you) 
    dF = derivative(F_act, you)

    problem = NonlinearVariationalProblem(F_act, you, bcs, dF)
    solver = NonlinearVariationalSolver(problem)
    # need to adjust solver parameters b/c PetSC runs out of memory if default solver is used. 
    prm = solver.parameters
    prm["nonlinear_solver"] = "newton"
    prm["newton_solver"]["linear_solver"] = "mumps" # utilizes parallel processors
    try:
        (iters, converged) = solver.solve()
    except:
        print("Newton Solver failed to converge")
        converged = False
        iters = -1
    
    u1, p1, Sigma, D1_vec = you.split(deepcopy = True)
    D1 = as_tensor([[D1_vec[0], D1_vec[1]], [D1_vec[1], -D1_vec[0]]]) # Reshape the strain/velocity gradient tensor

    #remember, Sigma here is not the stress tensor, but the "modified" tensor Sigma. Return both
    #stress = project(Sigma + 2*eta*D1, FunctionSpace(mesh, TensorElement("CG", triangle, 1, symmetry=True))) # project onto deg 1 space
    
    return Results(converged, u1, D1, p1, Sigma, iters)
    

def oldroyd_3_LDC3D_EVSS(h, s, eta, l1, mu1):
    nx = round(1/h)
    mesh = UnitCubeMesh(nx, nx, nx)
    print("Mesh loaded into FEniCS")

    # boundary data
    g_top = Expression(("s*256.0*x[0]*x[0]*x[1]*x[1]*(1-x[0])*(1-x[0])*(1-x[1])*(1-x[1])", "0.0", "0.0"), s=s, degree = 4) 
    #g_top = Constant((float(s), 0.0))
    g_walls = Constant((0.0, 0.0, 0.0)) #g=0 on walls
    
    # body forces
    f = Constant((0.0, 0.0, 0.0)) # no body forces
    
    # Element spaces
    print("Creating element spaces...")
    V_elem = VectorElement("CG", tetrahedron, 2) # Velocity, degree 2 elements
    P_elem = FiniteElement("CG", tetrahedron, 1) # Pressure, degree 1 elements
    T_elem = TensorElement("CG", tetrahedron, 2, symmetry=True) # elastic "Stress" tensor tau, degree 2 elements
    D_elem = VectorElement("CG", tetrahedron, 1, dim=5) # Deformation tensor, defined as VectorElement to exploit symmetry

    W_elem = MixedElement([V_elem, P_elem, T_elem, D_elem]) # Mixed element (u, p, Sigma, D_Vec)
    
    # Function Spaces
    V_space = FunctionSpace(mesh, V_elem) # need these first 4 for FunctionAssigner
    P_space = FunctionSpace(mesh, P_elem)
    T_space = FunctionSpace(mesh, T_elem)
    D_space = FunctionSpace(mesh, D_elem)
    W = FunctionSpace(mesh, W_elem) # Only need this one for pure EVSS, the mixed Function Space
    print("Done.")

    # Interpolate body force and BCs onto velocity FE space
    print("Interpolating data...")
    g_top = interpolate(g_top, W.sub(0).collapse())
    g_walls = interpolate(g_walls, W.sub(0).collapse())
    f = interpolate(f, W.sub(0).collapse())
    print("Done.")

    # Define boundary conditions
    top_lid = 'near(x[2], 1.0) && on_boundary'
    walls = '(near(x[0], 0.0) || near(x[0], 1.0) || near(x[1], 0.0) || near(x[1], 1.0) || near(x[2], 0.0)) && on_boundary'
    origin = 'near(x[0], 0.0) && near(x[1], 0.0) && near(x[2], 0.0)' #for pressure regulating

    bc_top = DirichletBC(W.sub(0), g_top, top_lid)
    bc_walls = DirichletBC(W.sub(0), g_walls, walls)
    pressure_reg = DirichletBC(W.sub(1), Constant(0.0), origin, 'pointwise')
        
    # Gather boundary conditions (any others would go here, separated by a comma)
    bcs = [bc_top, bc_walls, pressure_reg] 

    # Variational Problem: Trial and Test Functions
    w = TrialFunction(W) 
    (u, p, Tau, D_vec) = split(w) #trial/solution functions
    (v, q, S, Phi_vec) = TestFunctions(W) #test functions

    # D is symmetric and traceless...
    D = as_tensor([[D_vec[0], D_vec[1], D_vec[2]], [D_vec[1], D_vec[3], D_vec[4]], [D_vec[2], D_vec[4], -D_vec[0]-D_vec[3]]])
    # ... and so are test funcs for enforcement  
    Phi = as_tensor([[Phi_vec[0], Phi_vec[1], Phi_vec[2]], [Phi_vec[1], Phi_vec[3], Phi_vec[4]], [Phi_vec[2], Phi_vec[4], -Phi_vec[0]-Phi_vec[3]]])

    # Momentum equation gets velocity Test Functs v
    # EVSS SF: -eta*Lapl(u) + del(u)*u + del(p) - del.Sigma - f = 0
    momentum = eta*inner(grad(u), grad(v))*dx + inner(dot(grad(u), u) + grad(p) - div(Tau) - f, v)*dx
    
    # Continuity equation gets pressure TFs q
    # SF: del.u = 0
    continuity = div(u)*q*dx
    
    # Constitutive equation gets degree 2 tensor TFs and SUPG (S + h*del(S)*u)
    # CE SF: T + l1*G(T + 2etaE, a) = 0
    constitutive = inner( Tau + l1*(dot(grad(Tau+2*eta*D), u) - dot(grad(u), Tau+2*eta*D) - dot(Tau+2*eta*D, grad(u).T)) \
                         +(l1-mu1)*(dot( (grad(u)+grad(u).T)/2, Tau+2*eta*D) + dot(Tau+2*eta*D,(1/2)*(grad(u)+grad(u).T))),\
                         S + h*dot(grad(S), u) )*dx
    
    # Enforcement of deformation tensor gets gets degree 1 tensor TFs Phi
    enforcement = inner(D - (1/2)*(grad(u) + grad(u).T), Phi)*dx
    
    # Sum all together
    F = momentum + continuity + constitutive + enforcement
    
    # Solve discrete variational form
    you = Function(W) #starting guess for Newton, default Function() is constant zeros
    
    #get NSE as starting guess for Newton
    print("Getting Navier-Stokes solution as starting guess...")
    nse_solution = steady_nse_solver.navier_stokes_LDC3D(h, s, eta) 
    print("NSE Solve done")

    #u_nse = nse_solution.velocity
    #p_nse = nse_solution.pressure
    u_nse = interpolate(nse_solution.velocity, V_space)
    p_nse = interpolate(nse_solution.pressure, P_space)

    # Found that just using velocity and pressure works best, for some reason
    D_start = interpolate(Constant((0.0, 0.0, 0.0, 0.0, 0.0)), D_space)
    Sigma_start = interpolate( Constant( ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)) ), T_space)
    
    assigner = FunctionAssigner(W, [V_space, P_space, T_space, D_space])
    assigner.assign(you, [u_nse, p_nse, Sigma_start, D_start])

    # done with NSE solution, which is stored in "you" compliant with the mixed form of you
    
    F_act = action(F, you) 
    dF = derivative(F_act, you)

    problem = NonlinearVariationalProblem(F_act, you, bcs, dF)
    solver = NonlinearVariationalSolver(problem)
    # need to adjust solver parameters b/c PetSC runs out of memory if default solver is used. 
    prm = solver.parameters
    prm["nonlinear_solver"] = "newton"
    prm["newton_solver"]["linear_solver"] = "mumps" # utilizes parallel processors
    try:
        (iters, converged) = solver.solve()
    except:
        print("Newton Solver failed to converge")
        converged = False
        iters = -1
    
    u1, p1, Sigma, D1_vec = you.split(deepcopy = True)
    D1 = as_tensor([[D1_vec[0], D1_vec[1], D1_vec[2]], [D1_vec[1], D1_vec[3], D1_vec[4]], [D1_vec[2], D1_vec[4], -D1_vec[0]-D1_vec[3]]])
    
    #remember, Sigma here is not the stress tensor, but the "modified" tensor Sigma. Return both
    #stress = project(Sigma + 2*eta*D1, FunctionSpace(mesh, TensorElement("CG", tetrahedron, 1, symmetry=True))) # project onto deg 1 space
    
    return Results(converged, u1, D1, p1, Sigma, iters)



#post proc stuff here

