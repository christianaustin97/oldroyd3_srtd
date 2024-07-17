"""
Generates much needed plots for the Oldroyd 3 SRTS paper
"""

from fenics import *
import matplotlib.tri
from oldroyd_3_EVSS import *
from oldroyd_3_SRTD import *
from oldroyd_3_SRTD_SUPG import *

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib

# JB geometry parameters
rad = 0.5
ecc = 0.25

# Fluid parameters
eta = 1.0
lambda1 = 1e-2
mu1 = lambda1

# SRTD algorithm parameters
max_srtd_iters = 20
srtd_tol = 1e-9

# LDC
h = 1.25e-2

filename = "lid_driven_cavity_h_%.4e"%h
meshfile = "meshdata/" + filename + ".h5"
svgfile = "plots/" + filename + ".svg" # .svg vector graphics
pdffile = "plots/" + filename + ".pdf"

if not os.path.exists(meshfile):
    print("Creating mesh...")
    gen_mesh_jb.main(h, rad, ecc)
# Read the mesh in 
mesh = Mesh() #empty mesh
infile = HDF5File(MPI.comm_world, meshfile, 'r')
infile.read(mesh, '/mesh', True) #for some reason, need this True flag to import a mesh?
infile.close()

ldc_plot = plot(mesh)
plt.savefig(svgfile, bbox_inches='tight')
plt.savefig(pdffile, bbox_inches='tight')
plt.close()


# solution plots
solnfile = "plots/ldc_ucm_l1=%.3e"
u_srtd = oldroyd_3_LDC_SRTD(h, eta, lambda1, mu1, max_srtd_iters, srtd_tol)
u_evss = oldroyd_3_LDC_EVSS(h, eta, lambda1, mu1)

# SRTD Solution plots
fig = plot(u_srtd.velocity)
plt.title("Velocity")
plt.colorbar(fig, label = "magnitude")
plt.savefig(solnfile + "_srtd_u.pdf", bbox_inches = 'tight')
plt.close()

fig = plot(sqrt(inner(u_srtd.velocity, u_srtd.velocity))) # Log colorbar add argument , norm=matplotlib.colors.SymLogNorm(linthresh)
plt.title("Velocity Magnitude")
plt.colorbar(fig)
plt.savefig(solnfile + "_srtd_u_mag.pdf", bbox_inches = 'tight')
plt.close()

fig = plot(u_srtd.pressure)
plt.title("Pressure")
plt.colorbar(fig)
plt.savefig(solnfile + "_srtd_pressure.pdf", bbox_inches = 'tight')
plt.close()


# EVSS Solution plots
fig = plot(u_evss.velocity)
plt.title("Velocity")
plt.colorbar(fig, label = "magnitude")
plt.savefig(solnfile + "_evss_u.pdf", bbox_inches = 'tight')
plt.close()

fig = plot(sqrt(inner(u_evss.velocity, u_evss.velocity)))
plt.title("Velocity Magnitude")
plt.colorbar(fig)
plt.savefig(solnfile + "_evss_u_mag.pdf", bbox_inches = 'tight')
plt.close()

fig = plot(u_evss.pressure)
plt.title("Pressure")
plt.colorbar(fig)
plt.savefig(solnfile + "_evss_pressure.pdf", bbox_inches = 'tight')
plt.close()


# Trying some streamlines, must calculate the "stream function" psi
Q = FunctionSpace(mesh, "CG", 1)
psi_evss = Function(Q)
psi_srtd = Function(Q)
r = TrialFunction(Q)
s = TestFunction(Q)

# negative Laplacian of stream function is curl of u
velo_srtd = u_srtd.velocity
velo_evss = u_evss.velocity
curl_srtd = project(velo_srtd[1].dx(0) - velo_srtd[0].dx(1), Q)
curl_evss = project(velo_evss[1].dx(0) - velo_evss[0].dx(1), Q)

a = inner(grad(r), grad(s))*dx
L_srtd = curl_srtd*s*dx
L_evss = curl_evss*s*dx

wall = DirichletBC(Q, 0, "on_boundary")
solve(a == L_srtd, psi_srtd, wall)
solve(a == L_evss, psi_evss, wall)

triang = tri.Triangulation(*mesh.coordinates().reshape((-1, 2)).T,
                           triangles=mesh.cells())
Z_srtd = psi_srtd.compute_vertex_values(mesh)
Z_evss = psi_evss.compute_vertex_values(mesh)

plt.figure()
plt.tricontour(triang, Z_srtd, levels = 10, colors="black", linestyles="-")
plt.axis('square')
plt.title("Streamlines")
plt.savefig(solnfile + "_srtd_streamlines.pdf", bbox_inches = 'tight')
plt.close()

plt.figure()
plt.tricontour(triang, Z_evss, levels = 10, colors="black", linestyles="-")
plt.axis('square')
plt.title("Streamlines")
plt.savefig(solnfile + "_evss_streamlines.pdf", bbox_inches = 'tight')
plt.close()









# JB
h=2.5e-2

filename = "journal_bearing_h_%.4e"%h
meshfile = "meshdata/" + filename + ".h5"
svgfile = "plots/" + filename + ".svg" # .svg vector graphics
pdffile = "plots/" + filename + ".pdf"

if not os.path.exists(meshfile):
    print("Creating mesh...")
    gen_mesh_ldc.main(h)
# Read the mesh in 
mesh = Mesh() #empty mesh
infile = HDF5File(MPI.comm_world, meshfile, 'r')
infile.read(mesh, '/mesh', True) #for some reason, need this True flag to import a mesh?
infile.close()

ldc_plot = plot(mesh)
plt.savefig(svgfile, bbox_inches='tight')
plt.savefig(pdffile, bbox_inches='tight')
plt.close()

# solution plots
solnfile = "plots/jb_ucm_l1=%.3e"
u_srtd = oldroyd_3_JB_SRTD(h, rad, ecc, eta, lambda1, mu1, max_srtd_iters, srtd_tol)
u_evss = oldroyd_3_JB_EVSS(h, rad, ecc, eta, lambda1, mu1)

# SRTD Solution plots
fig = plot(u_srtd.velocity)
plt.title("Velocity")
plt.colorbar(fig, label = "magnitude")
plt.savefig(solnfile + "_srtd_u.pdf", bbox_inches = 'tight')
plt.close()

fig = plot(sqrt(inner(u_srtd.velocity, u_srtd.velocity)))
plt.title("Velocity Magnitude")
plt.colorbar(fig)
plt.savefig(solnfile + "_srtd_u_mag.pdf", bbox_inches = 'tight')
plt.close()

fig = plot(u_srtd.pressure)
plt.title("Pressure")
plt.colorbar(fig)
plt.savefig(solnfile + "_srtd_pressure.pdf", bbox_inches = 'tight')
plt.close()


# EVSS Solution plots
fig = plot(u_evss.velocity)
plt.title("Velocity")
plt.colorbar(fig, label = "magnitude")
plt.savefig(solnfile + "_evss_u.pdf", bbox_inches = 'tight')
plt.close()

fig = plot(sqrt(inner(u_evss.velocity, u_evss.velocity)))
plt.title("Velocity Magnitude")
plt.colorbar(fig)
plt.savefig(solnfile + "_evss_u_mag.pdf", bbox_inches = 'tight')
plt.close()

fig = plot(u_evss.pressure)
plt.title("Pressure")
plt.colorbar(fig)
plt.savefig(solnfile + "_evss_pressure.pdf", bbox_inches = 'tight')
plt.close()

# Trying some streamlines, must calculate the "stream function" psi
Q = FunctionSpace(mesh, "CG", 1)
psi_evss = Function(Q)
psi_srtd = Function(Q)
r = TrialFunction(Q)
s = TestFunction(Q)

# negative Laplacian of stream function is curl of u
velo_srtd = u_srtd.velocity
velo_evss = u_evss.velocity
curl_srtd = project(velo_srtd[1].dx(0) - velo_srtd[0].dx(1), Q)
curl_evss = project(velo_evss[1].dx(0) - velo_evss[0].dx(1), Q)

a = inner(grad(r), grad(s))*dx
L_srtd = curl_srtd*s*dx
L_evss = curl_evss*s*dx

wall = DirichletBC(Q, 0, "on_boundary")
solve(a == L_srtd, psi_srtd, wall)
solve(a == L_evss, psi_evss, wall)

triang = tri.Triangulation(*mesh.coordinates().reshape((-1, 2)).T,
                           triangles=mesh.cells())
Z_srtd = psi_srtd.compute_vertex_values(mesh)
Z_evss = psi_evss.compute_vertex_values(mesh)

plt.figure()
plt.tricontour(triang, Z_srtd, levels = 7, colors="black", linestyles="-")
plt.axis('square')
plt.title("Streamlines")
plt.savefig(solnfile + "_srtd_streamlines.pdf", bbox_inches = 'tight')
plt.close()

plt.figure()
plt.tricontour(triang, Z_evss, levels = 7, colors="black", linestyles="-")
plt.axis('square')
plt.title("Streamlines")
plt.savefig(solnfile + "_evss_streamlines.pdf", bbox_inches = 'tight')
plt.close()






