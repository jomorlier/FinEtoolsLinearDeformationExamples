module ESNICE_examples
using Statistics
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using PGFPlotsX
using Test

function ESNICE_energies()
    E = 1e6*phun("PA");
    nu = 0.3;
    L = 2*phun("M");
    hs = collect(linearspace(0.1*L, L, 20))
    mag = 0.001

    PEs = Float64[]
    APEs = Float64[]
    for h = hs
        xs = collect(linearspace(0.0, L, 2))
        ys = collect(linearspace(0.0, h, 2))
        zs = collect(linearspace(0.0, h, 2))
        fens, fes = T4blockx(xs, ys, zs, :a)
        MR = DeforModelRed3D
        geom = NodalField(fens.xyz)
        u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
        numberdofs!(u)

        l1 = selectnode(fens; plane = [1.0 0.0 0.0 0.0], thickness = h/1000)
        for i = l1
            x,y,z = geom.values[i, :]
            u.values[i, 1] = (z - h/2) * mag
        end
        l1 = selectnode(fens; plane = [1.0 0.0 0.0 L], thickness = h/1000)
        for i = l1
            x,y,z = geom.values[i, :]
            u.values[i, 1] = -(z - h/2) * mag
        end

        material = MatDeforElastIso(MR, E, nu)
        femm = FEMMDeforLinear(MR, IntegData(fes, TetRule(1)), material)
        associategeometry!(femm, geom)
        K = stiffness(femm, geom, u)
        U = gathersysvec(u)
        PE = dot(U, K * U) / 2.0
        push!(PEs, PE)
        I = h * h^3 / 12
        APE = 2 * E * I * mag^2 / L
        push!(APEs, APE)
    end

    r = hs ./ L
    @show rPE = PEs ./ APEs
    @pgf a = Axis({
            xlabel = "Aspect ratio",
            ylabel = "Relative Potential Energy",
            grid="major",
            legend_pos  = "north east"
        },
        Plot({mark="circle"}, Table([:x => vec(r), :y => vec(rPE)])))
    display(a)

    # fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)
    # File =  "mt4energy2.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    true
end


function ESNICE_vibration()
    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0*2*pi)^2;
    f = 4450*phun("Hz");
    omega = 2 * pi * f

    MR = DeforModelRed3D
    output = import_ABAQUS("alum_cyl.inp")
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegData(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    @pgf a = Axis({
            xlabel = "Entity",
            ylabel = "Stabilization factor",
            grid="major",
            legend_pos  = "north east"
        },
        Plot({mark="circle"}, Table([:x => vec(1:count(fens)), :y => vec(femm.nphis)])))
    display(a)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")

    vectors = []
    for i = 7:length(fs)
        scattersysvec!(u, v[:, i])
        push!(vectors, ("Mode_$i", deepcopy(u.values)))
    end
    File  =   "alum_cyl_mode_shapes.vtk"
    vtkexportmesh(File, connasarray(fes), fens.xyz, FinEtools.MeshExportModule.T4; vectors = vectors)
    @async run(`"paraview.exe" $File`)


    true
end # function

function allrun()
    println("#####################################################")
    println("# ESNICE_energies ")
    ESNICE_energies()
    println("#####################################################")
    println("# ESNICE_vibration ")
    ESNICE_vibration()
    return true
end # function allrun

end # module LE11NAFEMS_examples
