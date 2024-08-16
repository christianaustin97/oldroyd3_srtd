""" 
    Solves the steady-state incompressible Navier Stokes Equations
    for a few few common fluid flow problems using Taylor-Hood
    elements/a mixed formulation
"""

from fenics import *
from mshr import *
import math as math
import matplotlib.pyplot as plt
import sys

class Results: 
    def __init__(self, velocity, pressure, stress_tensor):
        self.velocity = velocity
        self.pressure = pressure
        self.stress_tensor = stress_tensor	


def navier_stokes_ldc(meshsize, eta):
    """
        solves the Steady NSE for the lid-driven cavity problem
    """
    
    # boundary data
    #g_top = Constant((1.0, 0.0)) #g.n=0 on the lid, or top wall
    g_top = Expression(("30.0*x[0]*x[0]*(1-x[0])*(1-x[0])", "0.0"), degree = 4) # 30x^2(1-x)^2, 30 gives it integral=1
    g_walls = Constant((0.0, 0.0)) #g=0 on walls

    # body forces
    f = Constant((0.0, 0.0)) # no body forces
    
    # Variational problem start: define domain and mesh
    nx = meshsize; ny = meshsize
    mesh = UnitSquareMesh(nx,ny, "crossed")

    # Refine mesh near top corners, where singularities develop for this problem
    for n in range(3): #refine 3 times
        tl_markers = MeshFunction('bool', mesh, mesh.topology().dim(), False)
        for cell in cells(mesh):
            tl_markers[cell] = cell.contains(Point(0.0, 1.0))
        mesh = refine(mesh, tl_markers)

        tr_markers = MeshFunction('bool', mesh, mesh.topology().dim(), False)
        for cell in cells(mesh):
            tr_markers[cell] = cell.contains(Point(1.0, 1.0))
        mesh = refine(mesh, tr_markers)
    
    # Element spaces
    P_elem = FiniteElement("CG", triangle, 1) #pressure and auxiliary pressure, degree 1 elements
    V_elem = VectorElement("CG", triangle, 2) #velocity, degree 2 elements
    T_elem = TensorElement("CG", triangle, 2) #tensor, degree 2 elements (only used for outputting stress Tensor)
    
    W_elem = P_elem * V_elem # Mixed/Taylor Hood element space for Navier-Stokes type equations

    W = FunctionSpace(mesh, W_elem) # Taylor-Hood/mixed space
    P = FunctionSpace(mesh, P_elem) # pressure space
    V = FunctionSpace(mesh, V_elem) # velocity space
    T = FunctionSpace(mesh, T_elem) # tensor space
    
    # Interpolate body force and BCs onto velocity FE space
    g_top = interpolate(g_top, W.sub(1).collapse())
    g_walls = interpolate(g_walls, W.sub(1).collapse())
    f = interpolate(f, W.sub(1).collapse())
    
    # Define boundary conditions
    lid    = 'near(x[1], 1.0) && on_boundary'
    walls  = '(near(x[1], 0.0) || near(x[0], 0.0) || near(x[0], 1.0)) && on_boundary'
    corner = 'near(x[0], 0.0) && near(x[1], 0.0)' # for pressure regulating
    
    bc_top   = DirichletBC(W.sub(1), g_top, lid) # driving lid
    bc_walls = DirichletBC(W.sub(1), g_walls, walls) # no slip
    bc_press = DirichletBC(W.sub(0), Constant(0.0), corner, 'pointwise') # pressure regulating
    
    # Gather boundary conditions (any others would go here, separated by a comma)
    bcs = [bc_top, bc_walls, bc_press] #possibly get rid of pressure regulator if it breaks something in NN solve
    
    # Variational Problem: Trial and Test Functions
    w = TrialFunction(W) # nonlinear in w
    (pi,u) = split(w)
    (q, v) = TestFunctions(W)
    
    # eta*<del(u), del(v) > + <del(u).u, v> - <pi, div(v)> + <q, div(u)> 
    a_nse = eta*inner(grad(u), grad(v))*dx + dot( dot(grad(u),u), v)*dx - (pi*div(v))*dx + q*div(u)*dx
    f_nse = inner(f,v)*dx
    
    F = a_nse - f_nse
    
    you = Function(W)
    F_act = action(F, you) 
    dF = derivative(F_act, you)

    problem = NonlinearVariationalProblem(F_act, you, bcs, dF)
    solver = NonlinearVariationalSolver(problem)
    try: 
        solver.solve()
    except: 
        print("Newton Method in the Navier-Stokes-like stage failed to converge")
        return Results(None, None, None)
    
    # Likewise, not sure which is preferred
    p_soln, u_soln = you.split(deepcopy=True)
    """ Or: 
    p_soln = you.sub(0) #pressure part of solution
    u_soln = you.sub(1)
    """
    
    #also return stress tensor value, for reference
    stress_tensor = project( 2*eta*(sym(grad(u_soln))) , T)
    
    print("L2 norm of stress tensor:")
    print(norm(stress_tensor, 'l2'))
    
    return Results(u_soln, p_soln, stress_tensor)
    

def navier_stokes_bearing(meshsize, rad, ecc, eta):
    """
        Solves flow around a journal/plain bearing
        
        Here, the larger circle, called the "bearing," will
        always have radius 1, center/axis (0,0), and the smaller
        circle, called the "journal," will have radius rad and 
        center/axis (0, -ecc)
        
    """
    if(rad>=1 or rad<=0 or ecc<0 or rad+ecc>1):
        #throw exception, forgot how lol
        print("Error: Inputs not valid")
    
    outer_circle = Circle(Point(0.0, 0.0), 1.0, 128) # Using n-gons to approx a circle lol
    inner_circle = Circle(Point(0.0, -ecc), rad, 128)
    domain = outer_circle - inner_circle
    mesh = generate_mesh(domain, meshsize)
    
    #boundaries of domain
    class Inner(SubDomain):
        def inside(self, x, on_boundary):
            radius = x[0]*x[0] + (x[1]+ecc)*(x[1]+ecc)
            return (on_boundary and (radius < rad*rad+DOLFIN_EPS))
    
    class Outer(SubDomain):
        def inside(self, x, on_boundary):
            #radius = x[0]*x[0] + (x[1]+ecc)*(x[1]+ecc)
            # feels like a hack but nothing else would work. Inner() worked fine
            #return (on_boundary and (radius > rad*rad+DOLFIN_EPS)) 
            radius = x[0]*x[0] + x[1]*x[1]
            return (on_boundary and near(radius, 1, 1e-2))
    
    class TopPoint(SubDomain):
        def inside(self, x, on_boundary):
            return (near(x[0], 0.0) and near(x[1], 1.0))
    
    # Boundary data     
    speed_outer = 0.0 # counter-clockwise tangential speed of outer bearing, 
    speed_inner = 1.0 # clockwise tangential speed of inner bearing, 
    g_inner = Expression(("(s*x[1]+ecc)/r", "-s*x[0]/r"), s=speed_inner, r=rad, ecc=ecc, degree=1) 
    g_outer = Expression(("-s*x[1]", "s*x[0]"), s=speed_outer, degree=1) # no slip on outer bearing 
    f = Constant((0.0, 0.0)) # no body forces
    
    # Element spaces
    P_elem = FiniteElement("CG", triangle, 1) #pressure and auxiliary pressure, degree 1 elements
    V_elem = VectorElement("CG", triangle, 2) #velocity, degree 2 elements
    T_elem = TensorElement("CG", triangle, 2) #tensor, degree 2 elements (only used for outputting stress Tensor)
    
    W_elem = P_elem * V_elem # Mixed/Taylor Hood element space for Navier-Stokes type equations

    W = FunctionSpace(mesh, W_elem) # Taylor-Hood/mixed space
    P = FunctionSpace(mesh, P_elem) # pressure space
    V = FunctionSpace(mesh, V_elem) # velocity space
    T = FunctionSpace(mesh, T_elem) # tensor space
    
    # Interpolate body force and BCs onto velocity FE space
    g_inner = interpolate(g_inner, W.sub(1).collapse())
    g_outer = interpolate(g_outer, W.sub(1).collapse())
    f = interpolate(f, W.sub(1).collapse())

    # Define boundary conditions
    bc_inner = DirichletBC(W.sub(1), g_inner, Inner())
    bc_outer = DirichletBC(W.sub(1), g_outer, Outer())
    bc_press = DirichletBC(W.sub(0), Constant(0.0), TopPoint(), 'pointwise')
    
    # Gather boundary conditions (any others would go here, separated by a comma)
    bcs = [bc_inner, bc_outer, bc_press] 
    
    # Variational Problem: Trial and Test Functions
    w = TrialFunction(W) # nonlinear in w
    (pi,u) = split(w)
    (q, v) = TestFunctions(W)
    
    # eta*<del(u), del(v) > + <del(u).u, v> - <pi, div(v)> + <q, div(u)> 
    a_nse = eta*inner(grad(u), grad(v))*dx + dot( dot(grad(u),u), v)*dx - (pi*div(v))*dx + q*div(u)*dx
    f_nse = inner(f,v)*dx
    
    F = a_nse - f_nse
    
    you = Function(W)
    F_act = action(F, you) 
    dF = derivative(F_act, you)

    problem = NonlinearVariationalProblem(F_act, you, bcs, dF)
    solver = NonlinearVariationalSolver(problem)
    try: 
        solver.solve()
    except: 
        print("Newton Method in the Navier-Stokes-like stage failed to converge")
        return Results(None, None, None)
    
    # Likewise, not sure which is preferred
    p_soln, u_soln = you.split(deepcopy=True)
    
    #also return stress tensor value, for reference
    stress_tensor = project( 2*eta*(sym(grad(u_soln))) , T)
    
    print("L2 norm of stress tensor:")
    print(norm(stress_tensor, 'l2'))
    
    return Results(u_soln, p_soln, stress_tensor)
    
    
# Post-processing/new function stuff here
    
  
