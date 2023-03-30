using Gridap
using GridapGmsh
using Gmsh
import GeometryBasics: Point, Point2
"""
`koch_subdivide(N, p1, p2)`

Creates a Koch fractal curve of order `N` between points `p1` and `p2`
using the Koch fractal curve. The number of refinements
is given by parameter `N`. The case `N = 0` corresponds
to one straight segment conecting `p1` to `p2`.
"""
function koch_subdivide(N, p1, p2)
    δp = p2-p1
    L = norm(δp)
    if N==0
	return [p1,p2]
    else
	pm1 = p1 + δp/3
	pm2 = pm1 + δp/3
	L₁ = L / 3
	pc = p1 + δp/2
	pv = pc + sqrt(3)/2 * L₁ * Point(-δp[2]/L, δp[1]/L)
	plst1 = koch_subdivide(N-1, p1, pm1)
	plst2 = koch_subdivide(N-1, pm1, pv)
	plst3 = koch_subdivide(N-1, pv, pm2)
	plst4 = koch_subdivide(N-1, pm2, p2)
        
	return [plst1; plst2[2:end]; plst3[2:end]; plst4[2:end]]
    end
    
end

"""
`kochmeshpts(N, L, W, Nw)`

Creates a square with length `L` and width `W` where the
top segment is composed of `Nw` koch curves of order `N`.
"""
function kochmeshpts(N, L, W, Nw)
    p1 = Point(0.0, 0.0)
    p2 = Point(W/Nw, 0.0)
    pp = koch_subdivide(N, p1, p2)
    pts = [Point2{Float64}(W, -L), Point2{Float64}(0.0, -L), Point(0.0, 0.0)]
    for i in 1:Nw
	append!(pts, pp[2:end] .+ Point((W/Nw)*(i-1), 0.0))
    end
    
    return pts
end


function kochmesh(basename, pts, xsize)
	npts = length(pts)
	gmsh.initialize()
	gmsh.model.add("kochrough")
	for (i,p) in enumerate(pts)
		gmsh.model.geo.addPoint(p[1], p[2], 0.0, xsize[i], i)
	end
	for i in 1:(npts-1)
		gmsh.model.geo.addLine(i, i+1, i)
	end
	gmsh.model.geo.addLine(npts, 1, npts)
	gmsh.model.geo.addCurveLoop(collect(1:npts), 1)
	gmsh.model.geo.addPlaneSurface([1], 1)
	gmsh.model.geo.synchronize()
	
	gmsh.model.addPhysicalGroup(0, [1,2], npts+10)
	gmsh.model.addPhysicalGroup(1, [npts], npts + 11)
	gmsh.model.addPhysicalGroup(1, [1], npts + 12)
	gmsh.model.addPhysicalGroup(1, [2], npts + 13)
	gmsh.model.addPhysicalGroup(1, collect(3:npts-1), npts + 14)

	gmsh.model.addPhysicalGroup(2, [1], npts + 15)
	
	#gmsh.model.setPhysicalName(2, 5, "domain")
	gmsh.model.setPhysicalName(0, npts + 10, "inner")
	gmsh.model.setPhysicalName(1, npts + 11, "symmright")
	gmsh.model.setPhysicalName(1, npts + 12, "inner")
	gmsh.model.setPhysicalName(1, npts + 13, "symmleft")
	gmsh.model.setPhysicalName(1, npts + 14, "outer")
	gmsh.model.setPhysicalName(1, npts + 15, "domain")
	
	
	#Gmsh.model.addPhysicalGroup()
	gmsh.model.mesh.generate(2)
	fmesh = "$(basename).msh"
	gmsh.write(fmesh)
	gmsh.finalize()
	return fmesh
end


function calcreflen(pts)
	np = length(pts)
	L = zeros(np)
	for i in 2:np
		L[i-1] = norm(pts[i]-pts[i-1])
	end
	L[end] = norm(pts[end]-pts[1])

	X = zeros(np)

	for i in 2:np
		X[i] = min(L[i-1], L[i])
	end
	X[1] = min(L[end], L[1])
	return X
end

function solve_dirichlet(basename, dvals=[1.0, 0.0];
                         model=model, reffe=reffe, Ω=Ω, dΩ=dΩ)
    
	V = TestFESpace(model, reffe, dirichlet_tags=["inner", "outer"])
	U = TrialFESpace(V, dvals)
	a(u,v) =  ∫(∇(v)⋅∇(u))dΩ
	l(v) = 0
	op = AffineFEOperator(a, l, U, V)
	u = solve(op)
	writevtk(Ω, basename, cellfields=["u"=>u])
	return u
end


function solve_neuman(basename, Ti=1.0, q=1.0;
                      model=model, reffe=reffe, Ω=Ω, dΩ=dΩ, Γ=Γ, dΓ=dΓ)
    
    V = TestFESpace(model, reffe, dirichlet_tags=["inner"])
    U = TrialFESpace(V, [Ti])
    a(u,v) =  ∫(∇(v)⋅∇(u))dΩ
    l(v) = ∫(v*q)*dΓ
    op = AffineFEOperator(a, l, U, V)
    u = solve(op)
    writevtk(Ω, basename, cellfields=["u"=>u])
    return u
end
                      


function poly_area(pts::AbstractVector{Point2{T}}) where {T}
    A = 0.0
    npts = length(pts)
    if npts <= 2
        return 0.0
    end

    for i in 1:npts-1
        A += pts[i][1] * pts[i+1][2]
        A -= pts[i+1][1] * pts[i][2]
    end
    A += pts[end][1]*pts[1][2]
    A -= pts[end][2]*pts[1][1]

    return abs(A/2)
end

function poly_length(pts::AbstractVector{Point2{T}}) where {T}
    L = 0.0
    npts = length(pts)
    for i in 2:length(pts)
	L += norm(pts[i] - pts[i-1])
    end
    L += norm(pts[end]-pts[1])
    return L
end

        
