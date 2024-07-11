"""
    Christian Austin, University of Florida
    Part of Research for PhD Thesis, Summer of 2024
    
    Implements the SRTD formulation and recommended iterative 
    algorithm designed by Scott and Girault (2021) for the
    steady flow of a non-Newtonian fluid governed by a certain 
    3-parameter subset of the Oldroyd 8-parameter model 
    (Oldroyd, 1958)
    
    This program utilizes a variational approach and legacy
    (2019) FEniCS. The algorithm is iterative, with each 
    iteration containing 3-stages: the first stage involves 
    solving a Navier-Stokes-like equation for u and the 
    auxiliary pressure pi, then a linear transport equation 
    for the true pressure, and finally a linear transport 
    equation for the stress tensor. 
    
    This version, 'oldroyd_3_SRTD_SUPG', implements streamline
    upwinding Petrov-Galerkin (SU/PG or SUPG) for the pressure 
    and stress transport equations. This is a common technique
    used to help stabilize advection transport equations. 
    
    This file contains built-in functions for solving the 
    lid-driven cavity (ldc) problem and the journal-bearing (jb)
    problem, as their analysis requires tangential Dirichlet 
    boundary conditions. Hopefully more geometries to come soon. 
"""

from fenics import *
from meshdata import gen_mesh_jb
from meshdata import gen_mesh_ldc
import sys
import os

class Results:
    def __init__(self, converged, velocity, aux_pressure, pressure, stress_tensor, num_iters, rel_error):
        self.converged = converged
        self.velocity = velocity
        self.aux_pressure = aux_pressure
        self.pressure = pressure
        self.stress_tensor = stress_tensor
        self.num_iters = num_iters
        self.rel_error = rel_error


# Journal Bearing Problem

def oldroyd_3_JB_SRTD_SUPG(h, rad, ecc, eta, l1, mu1, max_iter, tol):
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
    P_elem = FiniteElement("CG", triangle, 1) #pressure and auxiliary pressure, degree 1 elements
    V_elem = VectorElement("CG", triangle, 2) #velocity, degree 2 elements
    T_elem = TensorElement("CG", triangle, 2, symmetry=True) #tensor, degree 2 elements (only used for outputting stress Tensor)
    
    W_elem = MixedElement([V_elem, P_elem]) # Mixed/Taylor Hood element space for Navier-Stokes type equations

    W = FunctionSpace(mesh, W_elem) # Taylor-Hood/mixed space
    P = FunctionSpace(mesh, P_elem) # true pressure space
    V = FunctionSpace(mesh, V_elem) # velocity space
    T = FunctionSpace(mesh, T_elem) # tensor space
    
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
    w = TrialFunction(W) #will be our NS-like TrialFunction
    (u,pi) = split(w) 
    (v, q) = TestFunctions(W)

    p = TrialFunction(P) #true pressure trial function for auxiliary pressure equation
    r = TestFunction(P)

    tau = TrialFunction(T) #stress trial function for stress tensor equation
    S = TestFunction(T)
    
    # Initial guesses for iterative solve
    p0 = Constant(0.0) # true pressure
    p0 = interpolate(p0, P)

    u0 = Constant((0.0, 0.0)) # velocity
    u0 = interpolate(u0, W.sub(0).collapse())

    T0 = Constant( ((0.0, 0.0),(0.0, 0.0)) ) # stress tensor
    T0 = interpolate(T0, T)
    
    pi1 = Constant(0.0) #not ever used, just so we're never returning an empty var in the event of failure
    pi1 = interpolate(pi1, W.sub(1).collapse())
    
    #LHS of NS-like solve is same throughout, I think this can exist outside the loop
    a_nse = eta*inner(grad(u), grad(v))*dx + dot( dot(grad(u),u), v)*dx - (pi*div(v))*dx + q*div(u)*dx
    
    you = Function(W)
    
    # Begin SRTD iterative solve
    n=0; l2diff = 1.0
    while(n<=max_iter and l2diff > tol):
        # Stage 1: solve NS type equation for u(n) and pi(n) given u(n-1), p(n-1), and T(n-1)
        #   Test functions:  v and q
        #   Trial functions: u and pi
        
        # RHS of NS-like stage is given in section 7 of Girault/Scott paper
        E0 = sym(grad(u0))
        term1 = inner(f, v - l1*dot(grad(v), u0))*dx #orange term
        term2 = (p0*inner(nabla_grad(u0), grad(v)))*dx  #blue term
        term3 = -inner( dot(grad(u0),u0) , dot(grad(v),u0) )*dx #red term
        term4 = inner( dot(grad(u0),T0) , grad(v) )*dx #light green term
        term5 = inner( dot(E0,T0)+dot(T0,E0) , grad(v) )*dx #dark green term
        
        L_nse = term1 - l1*(term2 + term3 + term4) + (l1-mu1)*term5 #mathcal F 
        
        # Nonlinear in u, so must solve a-L==0 and use Newton instead of a==L directly
        F = a_nse - L_nse
        
        #Get Jacobian/Gateaux derivative
        #you = Function(W) # Uncomment this if you want to start the Newton solve from 0 every iter
        F_act = action(F, you) 
        dF = derivative(F_act, you)

        problem = NonlinearVariationalProblem(F_act, you, bcs, dF)
        solver = NonlinearVariationalSolver(problem)
        # need to adjust solver parameters b/c PetSC runs out of memory if default solver is used. 
        prm = solver.parameters
        prm["nonlinear_solver"] = "newton"
        prm["newton_solver"]["linear_solver"] = "mumps" # utilizes parallel processors
        try: 
            solver.solve()
        except: 
            print("Newton Method in the Navier-Stokes-like stage failed to converge")
            return Results(False, u0, pi1, p0, T0, n, None)
        
        #extract solution parts
        u1, pi1 = you.split(deepcopy=True)
        
        
        ###########################################################################
        # These next two are transport equations, meaning that they stand to 
        #   benefit from using Streamline upwinding (SUPG)
        
        # Stage 2: solve scalar transport equation for p(n) given pi(n), u(n)
        #   Trial function: p
        #   Test function: r + h*dot(del(r), u1)
        
        ap = (p + l1*dot(grad(p), u1)) * (r + h*dot(grad(r), u1))* dx 
        Lp = pi1 * (r + h*dot(grad(r), u1))* dx 
        
        p1 = Function(P)
        solve(ap == Lp, p1)
        
        # Step 3: solve Tensor equation for T(n) given u(n), also linear in T 
        #   Trial function: tau 
        #   Test function: S + h*dot(del(S), u1)
        
        E1 = sym(grad(u1))
        R1 = -skew(grad(u1)) #or skew(nabla_grad(u1)), the true vorticity tensor when using physicist notation
        
        aT = inner( tau + l1*(dot(grad(tau),u1) + dot(R1, tau) - dot(tau, R1)) \
                        - mu1*(dot(E1, tau) + dot(tau, E1)) , S + h*dot(grad(S), u1))*dx
        LT = 2.0*eta*inner(E1, S + h*dot(grad(S),u1) )*dx
        
        T1 = Function(T)
        solve(aT == LT, T1)
        
        # End of this iteration
        l2diff = errornorm(u1, u0, norm_type='l2')
        print("SRTD Iteration %d: r = %.4e (tol = %.3e)" % (n, l2diff, tol))
        n = n+1
        
        #update
        p0 = p1
        u0 = u1
        T0 = T1     
        
    # Stuff to do after the iterations are over
    if(l2diff <= tol):
        converged = True
    else:
        converged = False
    return Results(converged, u1, pi1, p1, T1, n-1, l2diff)


# Lid-Driven Cavity Problem

def oldroyd_3_LDC_SRTD_SUPG(h, eta, l1, mu1, max_iter, tol):
    
    meshfile = "meshdata/lid_driven_cavity_h_%.4e.h5"%h
    
    if not os.path.exists(meshfile):
        print("Creating mesh...")
        gen_mesh_ldc.main(h)
    
    #then, simply read the mesh in 
    mesh = Mesh() #empty mesh
    infile = HDF5File(MPI.comm_world, meshfile, 'r')
    infile.read(mesh, '/mesh', True) #for some reason, need this flag to import a mesh?
    infile.close()
    print("Mesh loaded into FEniCS")

    # boundary data
    g_top = Expression(("30.0*x[0]*x[0]*(1-x[0])*(1-x[0])", "0.0"), degree = 4) # 30x^2(1-x)^2, 30 gives it integral=1
    g_walls = Constant((0.0, 0.0)) #g=0 on walls

    # body forces
    f = Constant((0.0, 0.0)) # no body forces
    
    # Element spaces
    P_elem = FiniteElement("CG", triangle, 1) #pressure and auxiliary pressure, degree 1 elements
    V_elem = VectorElement("CG", triangle, 2) #velocity, degree 2 elements
    T_elem = TensorElement("CG", triangle, 2, symmetry=True) #stress tensor, degree 2 elements 
    
    W_elem = MixedElement([V_elem, P_elem]) # Mixed/Taylor Hood element space for Navier-Stokes type equations

    W = FunctionSpace(mesh, W_elem) # Taylor-Hood/mixed space
    P = FunctionSpace(mesh, P_elem) # true pressure space
    V = FunctionSpace(mesh, V_elem) # velocity space
    T = FunctionSpace(mesh, T_elem) # tensor space
    
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
    w = TrialFunction(W) #will be our NS-like TrialFunction
    (u,pi) = split(w) #trial functions
    (v, q) = TestFunctions(W)

    p = TrialFunction(P) #true pressure trial function for auxiliary pressure equation
    r = TestFunction(P)

    tau = TrialFunction(T) #stress trial function for stress tensor equation
    S = TestFunction(T)
    
    # Initial guesses for iterative solve
    p0 = Constant(0.0) # true pressure
    p0 = interpolate(p0, P)

    u0 = Constant((0.0, 0.0)) # velocity
    u0 = interpolate(u0, W.sub(0).collapse())

    T0 = Constant( ((0.0, 0.0),(0.0, 0.0)) ) # stress tensor
    T0 = interpolate(T0, T)
    
    pi1 = Constant(0.0) #not ever used, just so we're never returning an empty var in the event of failure
    pi1 = interpolate(pi1, W.sub(1).collapse())
    
    #LHS of NS-like solve is same throughout, I think this can exist outside the loop
    a_nse = eta*inner(grad(u), grad(v))*dx + dot( dot(grad(u),u), v)*dx - (pi*div(v))*dx + q*div(u)*dx
    
    you = Function(W)
    
    # Begin SRTD iterative solve
    n=1; l2diff = 1.0
    while(n<=max_iter and l2diff > tol):
        # Stage 1: solve NS type equation for u(n) and pi(n) given u(n-1), p(n-1), and T(n-1)
        # Test functions:  v and q
        # Trial functions: u and pi
        
        # RHS of NS-like stage is given in section 7 of Girault/Scott paper
        E0 = sym(grad(u0))
        term1 = inner(f, v - l1*dot(grad(v), u0))*dx #orange term
        term2 = (p0*inner(nabla_grad(u0), grad(v)))*dx  #blue term
        term3 = -inner( dot(grad(u0),u0) , dot(grad(v),u0) )*dx #red term
        term4 = inner( dot(grad(u0),T0) , grad(v) )*dx #light green term
        term5 = inner( dot(E0,T0)+dot(T0,E0) , grad(v) )*dx #dark green term
        
        L_nse = term1 - l1*(term2 + term3 + term4) + (l1-mu1)*term5 #mathcal F 
        
        # Nonlinear in u, so must solve a-L==0 and use Newton instead of a==L directly
        F = a_nse - L_nse
        
        #Get Jacobian/Gateaux derivative
        #you = Function(W) # Uncomment this if you want to start the Newton solve from 0 every iter
        F_act = action(F, you) 
        dF = derivative(F_act, you)

        problem = NonlinearVariationalProblem(F_act, you, bcs, dF)
        solver = NonlinearVariationalSolver(problem)
        # need to adjust solver parameters b/c PetSC can run out of memory if default solver is used
        prm = solver.parameters
        prm["nonlinear_solver"] = "newton"
        prm["newton_solver"]["linear_solver"] = "mumps" # utilizes parallel processors
        try: 
            solver.solve()
        except: 
            print("Newton Method in the Navier-Stokes-like stage failed to converge")
            return Results(False, u0, pi1, p0, T0, n, None)
        
        #extract solution parts
        u1, pi1 = you.split(deepcopy=True)
        
        ###########################################################################
        # These next two are transport equations, meaning that they could possibly
        #   benefit from using Streamline upwinding (SUPG)
        #
        # Stage 2: solve scalar transport equation for p(n) given pi(n), u(n)
        # Trial function: p
        # Test function: r (nonbold 'v' in Scott paper)
        
        ap = (p + l1*dot(grad(p), u1)) * (r + h*dot(grad(r), u1))* dx 
        Lp = pi1 * (r + h*dot(grad(r), u1))* dx 
        
        p1 = Function(P)
        solve(ap == Lp, p1)
        
        # Step 3: solve Tensor equation for T(n) given u(n), also linear in T 
        #   Trial function: tau 
        #   Test function: S + h*dot(del(S), u1)
        
        E1 = sym(grad(u1))
        R1 = -skew(grad(u1)) #or skew(nabla_grad(u1)), the true vorticity tensor when using physicist notation
        
        aT = inner( tau + l1*(dot(grad(tau),u1) + dot(R1, tau) - dot(tau, R1)) \
                        - mu1*(dot(E1, tau) + dot(tau, E1)) , S + h*dot(grad(S), u1))*dx
        LT = 2.0*eta*inner(E1, S + h*dot(grad(S),u1) )*dx
        
        T1 = Function(T)
        solve(aT == LT, T1)
        
        # End of this iteration
        l2diff = errornorm(u1, u0, norm_type='l2')
        print("SRTD Iteration %d: r = %.4e (tol = %.3e)" % (n, l2diff, tol))
        n = n+1
        
        #update
        p0 = p1
        u0 = u1
        T0 = T1     
        
    # Stuff to do after the iterations are over
    if(l2diff <= tol):
        converged = True
    else:
        converged = False
    return Results(converged, u1, pi1, p1, T1, n-1, l2diff)

     
    
#post proc stuff here


