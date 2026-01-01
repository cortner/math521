#
# This file collects utility functions to 
# solve the Dirichlet problem - Delta u = f 
# subject to homogeneous Dirichlet boundary 
# conditions.
#

using SparseArrays, LinearAlgebra, Ferrite, LaTeXStrings, FerriteViz

function generate_square(N)
    grid = generate_grid(Triangle, (N, N))
    ∂Ω = union(getfaceset.(Ref(grid), ["left", "right", "top", "bottom"])...)
    return grid, ∂Ω
end


function modify_grid(grid, is_in_dom, is_on_bdry)
    # remove all the nodes and elements outside the Lshape
    new_inds = zeros(Int, length(grid.nodes))
    new_nodes = eltype(grid.nodes)[]

    idx = 0 
    for (i_n, n) in enumerate(grid.nodes)
        if is_in_dom(n.x)
            idx += 1
            new_inds[i_n] = idx
            push!(new_nodes, n)
        end        
    end

    new_cells = eltype(grid.cells)[]
    for cell in grid.cells
        rmid = sum(grid.nodes[n_i].x for n_i in cell.nodes ) / 3
        if ( all( new_inds[n_i] > 0   for n_i in cell.nodes )
             && is_in_dom(rmid) )
            new_cell = Triangle(ntuple(t -> new_inds[cell.nodes[t]], 3))
            push!(new_cells, new_cell)
        end
    end
    
    ∂Ω = Set{FaceIndex}()
    for (i_c, c) in enumerate(new_cells)
        nn = c.nodes 
        for (i_f, (i, j)) in enumerate([(nn[1], nn[2]), (nn[2], nn[3]), (nn[1], nn[3])])
            rmid = 0.5 * (new_nodes[i].x + new_nodes[j].x)
            if is_on_bdry(rmid)
                push!(∂Ω, FaceIndex(i_c, i_f))
            end
        end
    end
    
    return Grid(new_cells, new_nodes), ∂Ω    
end
    

function generate_lshape(N)
    
    function is_in_dom(r)
        x = r[1]; y = r[2] 
        t = 1e-14
        return ( (abs(x) <= 1+t) && (abs(y) <= 1+t) && 
                !( (x > t) && (y < - t) ) )
    end
    
    function is_on_bdry(r)
        x = r[1]; y = r[2] 
        t = 1e-14        
        # left and right, top bottom boundaries
        if (abs(x) ≈ 1) || (abs(y) ≈ 1)
            return true 
        end
        #    L vertical part             L horizontal part
        if ((x ≈ 0) && (y <= t)) || ((y ≈ 0) && (x >= - t)) 
            return true 
        end 
        return false 
    end

    # start with a square domain
    grid = generate_grid(Triangle, (N, N))
    return modify_grid(grid, is_in_dom, is_on_bdry)
end
    
function generate_wedge(N)
    
    function is_in_dom(r)
        x = r[1]; y = r[2] 
        t = 1e-14
        return ( (abs(x) <= 1+t) && (abs(y) <= 1+t) && 
                !( (x > t) && (y+x < - t) ) )
    end
    
    function is_on_bdry(r)
        x = r[1]; y = r[2] 
        t = 1e-14        
        # left and right, top bottom boundaries
        if (abs(x) ≈ 1) || (abs(y) ≈ 1)
            return true 
        end
        #    L vertical part             L horizontal part
        if ((x ≈ 0) && (y <= t)) || ((x+y ≈ 0) && (x >= - t)) 
            return true 
        end 
        return false 
    end
    
    # start with a square domain
    grid = generate_grid(Triangle, (N, N))
    return modify_grid(grid, is_in_dom, is_on_bdry)    
end



# function generate_lshape(N)
    
#     # start with a square domain
#     grid = generate_grid(Triangle, (N, N))

#     # remove all the nodes and elements outside the Lshape
#     delete_node = fill(false, length(grid.nodes))
#     new_inds = zeros(Int, length(grid.nodes))
#     new_nodes = eltype(grid.nodes)[]

#     idx = 0 
#     for (i_n, n) in enumerate(grid.nodes)
#         if n.x[1] > 0 && n.x[2] < 0
#             delete_node[i_n] = true
#         else
#             idx += 1
#             new_inds[i_n] = idx
#             push!(new_nodes, n)
#         end        
#     end

#     new_cells = eltype(grid.cells)[]
#     for cell in grid.cells
#         if all( !(delete_node[n_i]) for n_i in cell.nodes )
#             new_cell = Triangle(ntuple(t -> new_inds[cell.nodes[t]], 3))
#             push!(new_cells, new_cell)
#         end
#     end
    
#     # now construct the set of faces on the boundary
#     function on_bdry(r1, r2)
#         # check that it's a vertical or horizontal edge
#         if !( (r1[1] == r2[1]) || (r1[2] == r2[2]) )
#             return false
#         end
#         # left and right boundaries
#         if (abs(r1[1]) == abs(r2[1]) == 1)
#             return true 
#         end
#         # top and bottom boundaries 
#         if (abs(r1[2]) == abs(r2[2]) == 1)
#             return true 
#         end
#         # inside the L vertical part
#         if (r1[1] == r2[1] == 0) && (r1[2] <= 0) && (r2[2] <= 0)
#             return true 
#         end 
#         # inside the L horizontal part
#         if (r1[2] == r2[2] == 0) && (r1[1] >= 0) && (r2[1] >= 0)
#             return true 
#         end 
#         return false 
#     end
    
#     ∂Ω = Set{FaceIndex}()
#     for (i_c, c) in enumerate(new_cells)
#         nn = c.nodes 
#         for (i_f, (i, j)) in enumerate([(nn[1], nn[2]), (nn[2], nn[3]), (nn[1], nn[3])])
#             r1 = new_nodes[i].x 
#             r2 = new_nodes[j].x 
#             if on_bdry(r1, r2)
#                 push!(∂Ω, FaceIndex(i_c, i_f))
#             end
#         end
#     end
    
#     return Grid(new_cells, new_nodes), ∂Ω
# end





function assemble_global!(cellvalues, dh, K, ffun)
   # we pre-allocate the element stiffness matrix and element force vector
   # these will be passed to the element assembly to avoid many allocations
   n_basefuncs = getnbasefunctions(cellvalues)
   Ke = zeros(n_basefuncs, n_basefuncs)
   fe = zeros(n_basefuncs)
    
   # Allocate global force vector f
   f = zeros(ndofs(dh))
    
   # Create an assembler: this object knows how to write 
   # the local arrays Ke, fe into the global arrays K, f 
   assembler = start_assemble(K, f)
    
   # Loop over all cells; this is managed by the DOF handler 
   # since the `cell` comes with information about local 
   # DOFs attached.
   for cell in CellIterator(dh)
       # Reinitialize cellvalues for this cell
       # `cellvalues` has iterators attached that need 
       # to be reset. This seems unnecessary and probably 
       # just a poor code design decision. 
       reinit!(cellvalues, cell)
       # ==========================================
       # Compute element contribution; 
       # this is where the actual work happens
       assemble_element!(Ke, fe, cell, cellvalues, ffun)
       # ==========================================
       # local-to-global assemble Ke and fe into K and f
       assemble!(assembler, celldofs(cell), Ke, fe)
   end
   return K, f
end


function assemble_element!(Ke, fe, cell, cellvalues, ffun)
   # number of local basis functions    
   n_basefuncs = getnbasefunctions(cellvalues)
   # Reset the local arrays 
   fill!(Ke, 0); fill!(fe, 0)
   # precompute the cell coordinates 
   cell_coords = getcoordinates(cell)
    
   # Loop over quadrature points
   for i_q in 1:getnquadpoints(cellvalues)
       # Get the quadrature weight for the current quad point
       # this includes the det(F) terms. It can be thought of 
       # as the volume element (hence dΩ)
       dΩ = getdetJdV(cellvalues, i_q)
        
       # evaluate f at the quadrature point 
       ξ_q = spatial_coordinate(cellvalues, i_q, cell_coords)
       f_q = ffun(ξ_q)
        
       # Loop over test shape functions (basis functions)
       for i in 1:n_basefuncs
           # get the values v = ψ_i(ξ_q), ∇v = ∇ψ_i(ξ_q)
           v  = shape_value(cellvalues, i_q, i)
           ∇v = shape_gradient(cellvalues, i_q, i)
            
           # ∫_K f v dx
           # Add contribution to fe
           fe[i] += f_q * v * dΩ
            
           # Loop over trial shape functions
           for j in 1:n_basefuncs
               ∇u = shape_gradient(cellvalues, i_q, j)
               # Add contribution to Ke
               #  ∫_K ∇v ⋅ ∇u 
               Ke[i, j] += dot(∇v, ∇u) * dΩ
           end
       end
   end
   return Ke, fe
end


function setup_fem(p, grid, ∂Ω; pquad = 2*p)
    dim = 2
    ip = Lagrange{dim, RefTetrahedron, p}()
    ip_geo = Lagrange{dim, RefTetrahedron, 1}()
    qr = QuadratureRule{dim, RefTetrahedron}(pquad)
    cellvalues = CellScalarValues(qr, ip, ip_geo)
    dh = DofHandler(grid)
    add!(dh, :u, 1, ip)
    close!(dh)
    K = create_sparsity_pattern(dh)
    ch = ConstraintHandler(dh)
    dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)
    add!(ch, dbc)
    close!(ch)    
    return cellvalues, dh, ch, K
end

function solve_fem(cellvalues, dh, ch, K, ffun)
    K, f = assemble_global!(cellvalues, dh, K, ξ -> ffun(ξ));
    apply!(K, f, ch)
    u = K \ f;
    return u 
end


"""
- N : # grid pts in each coordinate direction, h = 1/N
- k : polynomial order of the FEM  (Lagrange element)
- ffun : forcing function
"""
function fem_errors(k, N, generate_dom, u_ex)
    grid, ∂Ω = generate_dom(N)
    ffun = x -> - tr( ForwardDiff.hessian(u_ex, x) )
    cellvalues, dh, ch, K = setup_fem(k, grid, ∂Ω)
    u = solve_fem(cellvalues, dh, ch, K, ffun)
    return compute_errors(cellvalues, dh, u, u_ex)
end

function compute_errors(cellvalues, dh, u, u_ex)
    n_basefuncs = getnbasefunctions(cellvalues)
    err_L2 = 0.0 
    err_H1 = 0.0 
    
    # loop over cells (= elements)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)        
        n_basefuncs = getnbasefunctions(cellvalues)
        cell_coords = getcoordinates(cell)
        
        # we also need the local degrees of freedom 
        u_cell = u[cell.dofs]
    
        vK = 0.0
        for i_q in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, i_q)
            ξ_q = spatial_coordinate(cellvalues, i_q, cell_coords)
            u_q   = u_ex(ξ_q)
            ∇u_q  = ForwardDiff.gradient(u_ex, ξ_q)
            uh_q  = function_value(cellvalues, i_q, u_cell)
            ∇uh_q = function_gradient(cellvalues, i_q, u_cell)
            err_L2 += dΩ * (u_q - uh_q)^2
            err_H1 += dΩ * norm(∇u_q - ∇uh_q)^2
        end
    end
    return sqrt(err_L2), sqrt(err_H1)
end