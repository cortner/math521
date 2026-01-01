# Translation of the old Matlab code to Julia 

# Some Standards for Node and Edge numbering
#
#      P1
#      /\
#   E3/  \ E2
#    /    \
#   /______\
# P2   E1   P3
#
#   o E1 is always the refinement edge (newest vertex bisection)
#   o the numbering of the vertices on interior edges is arbitrary, but
#     for boundary edges, the domain always has to be to the left!
#   o aTE(i, k) is the index of the i-th edge of the k-th element
#   o aET(1, i) is the index of the left-hand element neighbour
#     (this means that for boundary edges aET(2, i) == 0)
#     aET(2, i) is the index of the right-hand element neighbour
#

module AFEM 

using StaticArrays

const TAET = typeof( (left=1, right=2) )

mutable struct Triangulation{T}
   nodes::Vector{SVector{2, T}}
   triangles::Vector{Tuple{Int, Int, Int}}
   edges::Vector{Tuple{Int, Int}}
   aTE::Vector{Tuple{Int, Int, Int}} 
   aET::Vector{TAET}
   mark::Vector{Bool}
end 

function Triangulation(nodes::AbstractVector{SVector{2, T}}, 
                       triangles::Vector{Tuple{Int, Int, Int}}, 
                       fixbdry = nothing) where {T} 

   # compute and store the edge information
   edges, aTE, aET = edges_init(triangles) 
   # edge refinement markers
   mark = fill(false, length(aET))

   # initialize the finite element method
   # afem = afem.fem.initialize(afem);

   return Triangulation(nodes, triangles, edges, aTE, aET, mark)
end

function refine_uniform(mesh, nref)
   for n = 1:nref
      mesh.mark[:] .= true
      refine!(mesh)
   end
   return mesh 
end 

edge(i1::Integer, i2::Integer) = extrema((i1, i2))
edge(i12::Tuple) = edge(i12...)


"""
Computes edges and adjacency relations for the mesh. This routine
is only used for the initial mesh. Afterwards these relations are
all automatically updated during refinement.
"""
function edges_init(elements)
   # All edges with the interior edges duplicated
   E = Tuple{Int, Int}[] 
   for T in elements 
      for i = 1:3 
         push!(E, edge(T[i], T[mod1(i+1, 3)]))
      end
   end
   edges = sort(unique(E))
   nE = length(edges) 
   nT = length(elements)

   # make a lookup table for the edges 
   edge2index = Dict{Tuple{Int, Int}, Int}()
   for (i, e) in enumerate(edges)
      edge2index[e] = i 
   end   

   # find the adjacency relation between elements and edges 
   aTE = Tuple{Int, Int, Int}[]
   aET_pre = zeros(Int, (2, nE))
   for (i, T) in enumerate(elements)
      # edges are ordered such that the local edge j is opposite local node j
      ee = ((T[2], T[3]), (T[3], T[1]), (T[1], T[2]))
      ee_srt = edge.(ee)
      ee_idx = ntuple(a -> edge2index[ee_srt[a]], 3)
      push!(aTE, ee_idx)

      # Find the adjacency relation between E and T. The order
      # defines also whether an element is considered "inside" the
      # edge, or "outside" the edge. Boundary edges are always "inside".
      # The mesh normal (DG) always points outside).
      for j = 1:3
         # case 1: left
         if edges[aTE[i][j]] == T[mod1(j+1, 3)]
            aET_pre[1, aTE[i][j]] = i
         # case two: right
         else
            aET_pre[2, aTE[i][j]] = i;
         end
      end
   end
   aET = [ (aET_pre[1, i], aET_pre[2, i]) for i = 1:size(aET_pre, 2) ]

   # Finally, fix the direction of the edges. For all boundary edges we need 
   # to check whether they are oriented in such a way that the (unique)
   # neighbouring element is to the left.
   iE = findall(aET_pre[1, :] .== 0)
   for i in iE 
      edges[i] = (edges[i][2], edges[i][1])
      aET[i] = (aET[i][2], 0)
   end

   return edges, aTE, aET
end

end

##
using StaticArrays
X = [ SA[0.0, 0.0], SA[1.0, 0.0], SA[1.0, 1.0], SA[0.0, 1.0] ] 
T = [ (1, 2, 4), (4, 2, 3) ]
mesh = AFEM.Triangulation(X, T)

mesh

##

# %
# % Basic refinement routine
# %
# function afem = refine(afem)

# % do nothing if nothing needs to be done.
# if (nnz(afem.mark) == 0)
#   return;
# end

# % at most one mark allowed!
# afem.mark = min(afem.mark, 1);

# % get an admissible mesh marking by marking the refinement edge
# % in every element which is marked for refinement but for which
# % the refinement edge itself is not marked.
# rerun = true;
# while (rerun)
#   rerun = false;
#   markel = afem.mark(afem.aTE);
#   summark = sum(markel, 1);
#   iFix = find( (markel(1, :) == 0) & (summark > 0) );
#   if ~isempty(iFix)
#     rerun = true;
#     afem.mark(afem.aTE(1, iFix)) = 1;
#   end
# end

# % Now do plain bisections until no further elements need to be refined
# iEref = find(afem.mark > 0);
# while ~isempty(iEref)
#   % separate iEref into boundary and interior edges
#   iEref1 = iEref(afem.aET(2, iEref) > 0);
#   iEref2 = iEref(afem.aET(2, iEref) == 0);
#   % first the list of interior edges which can be refined.
#   if ~isempty(iEref1)
#     iEref1 = iEref1( (afem.aTE(1, afem.aET(1, iEref1)) == iEref1) & ...
#       ((afem.aTE(1, afem.aET(2, iEref1)) == iEref1)) );
#   end
#   % list of boundary edges which can be refined
#   if ~isempty(iEref2)
#     iEref2 = iEref2( afem.aTE(1, afem.aET(1, iEref2)) == iEref2 );
#   end
#   % if no edges are left for refinement, then something bad has happened
#   if (length(iEref1) + length(iEref2) == 0)
#     error(['[afemm2:refine()] Some edges remain to be refined, but ' ...
#       'none are admissible. This happens when the initial mesh does ' ...
#       'not satisfy the xxx property. Use xxx to fix this.']);
#   end
#   % bisect the edges which we have found
#   afem = bisection(afem, iEref1, iEref2);
#   % get list of remaining edfes to refine
#   iEref = find(afem.mark > 0);
# end

# end


# %
# % iEref1: interior refinement edges
# % iEref2: boundary refinement edges
# %
# %
# % Definition of the bisection:
# %
# %        P1
# %        /|\
# %     E3/ | \ E2
# %      /  e  \
# %     /   |   \
# %    /  R | L  \
# %   /     |     \
# %  /_____\|/_____\
# %  P2     Q      P3
# %        E1     
# %
# % o The edge e points from P1 to Q and in this relation we understand
# %   the left and right child.
# % o The left child will overwrite the parent element. The right child will
# %   be added at the end of the T array.
# % o For the edge bisection, we always keep the direction of the edge on
# %   on both portions
# %
# function afem = bisection(afem, iEref1, iEref2)

# % all refinement edges
# iEref = [iEref1, iEref2]; nEref = length(iEref);

# % quick sanity check.
# if isempty(iEref)
#   error('[afemm2:bisection()] iEref is empty but it shouldn''t be!');
# end

# nEref1 = length(iEref1); nEref2 = length(iEref2);
# nEnew = 3 * nEref1 + 2 * nEref2;
# nTnew = 2 * nEref1 + nEref2;
# nPnew = nEref1 + nEref2;

# %% I. Create the new vertices
# %%
# % create the nodes: just add the new nodes at the end of the P array
# % in iPnew, store the indices of the newly added nodes
# afemnew = empty_afemm2(afem.fem, afem.fixbdry);
# afemnew.nP0 = afem.nP0; afemnew.nT0 = afem.nT0;
# afemnew.nP = afem.nP + nPnew;
# afemnew.P = [ afem.P, ...
#   0.5 * (afem.P(:, afem.E(1, iEref)) + afem.P(:, afem.E(2, iEref))) ];
# iPnew = (afem.nP+1):(afemnew.nP);

# %% II. Create the edges
# %%
# % create the edges: first connect the new midpoints with the second
# % point and add them at the end, then connect the first point with the
# % new midpoint and overwrite the old edge with it.
# afemnew.nE = afem.nE + nEnew;
# afemnew.E = [ ...
#   ... % 1. the old edges (here we will later add one of the refinements)
#   [afem.E] ...
#   ... % 2. the second part of the bisected edge (keep the direction!)
#   [iPnew; afem.E(2, iEref)] ...
#   ... % 3. the left-hand edge bisection the element, make it point towards
#   ... % the new vertex!
#   [afem.T(1, afem.aET(1, iEref)); iPnew] ...
#   ... % 4. only for the interior refinement edges: add the edge bisection
#   ... % the right-hand element neighbour
#   [afem.T(1, afem.aET(2, iEref1)); iPnew(1:nEref1)] ];
# % fix the old edges, which are now bisected.
# afemnew.E(:, iEref) = [afem.E(1, iEref); iPnew];
# % update the refinement markers
# afemnew.mark = [afem.mark, zeros(1, nEnew)];
# afemnew.mark(iEref) = 0;
# % store some info for later
# iEnew_edge = [(afem.nE+1):(afem.nE+nEref)];
# iEnew_left = [(afem.nE+nEref+1):(afem.nE+2*nEref)];
# iEnew_right = [(afem.nE+2*nEref+1):(afem.nE+2*nEref+nEref1)];

# %% III. Create the new elements 
# %%
# iL = afem.aET(1, iEref);
# iR = afem.aET(2, iEref1);
# iTref = [iL, iR];
# afemnew.nT = afem.nT + nTnew;
# % create new elements
# afemnew.T = [ ...
#   ... % 1. the old elements (will later be corrected to the left child)
#   [afem.T] ...
#   ... % 2. the right child of the element
#   [[iPnew, iPnew(1:nEref1)]; afem.T(1, iTref); afem.T(2, iTref)] ];
# % 4. correct the left child
# afemnew.T(:, iTref) = ...
#   [[iPnew, iPnew(1:nEref1)]; afem.T(3, iTref); afem.T(1, iTref)];
# % store some index info for later
# iTnew_left = (afem.nT+1):(afem.nT+length(iL));
# iTnew_right = (afem.nT+length(iL)+1):(afemnew.nT);

# %% IV. Neighbourhood Information
# %%
# % aTE relation, which is safe to compute this way.
# afemnew.aTE = [ afem.aTE, ...
#   [afem.aTE(3, iL); iEref; iEnew_left], ...
#   [afem.aTE(3, iR); iEnew_edge(1:nEref1); iEnew_right] ];
# afemnew.aTE(:, iL) = [ afem.aTE(2, iL); iEnew_left; iEnew_edge ];
# if ~isempty(iR)
#   afemnew.aTE(:, iR) = [ afem.aTE(2, iR); iEnew_right; iEref1 ];
# end

# % the aET relation will be recovered from aTE: 
# afemnew.aET = zeros(2, afemnew.nE);
# for j = 1:3
#   % myfind returns the "find" indices in the first argument and
#   % the complement in the second argument.
# %   [ileft, iright] = ...
# %     myfind((afemnew.E(1, afemnew.aTE(j, :)) == afemnew.T(mod(j,3)+1, :)));
# %   afemnew.aET(1, afemnew.aTE(j, ileft)) = ileft;
# %   afemnew.aET(2, afemnew.aTE(j, iright)) = iright;
#   % The "pure Matlab" version.
#   ileft = find(afemnew.E(1, afemnew.aTE(j,:)) == afemnew.T(mod(j,3)+1,:));
#   afemnew.aET(1, afemnew.aTE(j, ileft)) = ileft;
#   iright = find(afemnew.E(1, afemnew.aTE(j,:)) ~= afemnew.T(mod(j,3)+1,:));
#   afemnew.aET(2, afemnew.aTE(j, iright)) = iright;
# end

# %% V. Refinement information and Projection:
# %%
# if ~isempty(afem.fem)
#   % We can think of each edge being subdivided into two edges and one
#   % vertex. We store this information in refE2PE
#   %   [ indices of edges in afem which are refined; ...
#   %     indices of new vertices on those edges; ...
#   %     indices of first edge built from old ones (=1st row); ...
#   %     indices of second edge built from old ones ];
#   refE2PE = [iEref; iPnew; iEref; iEnew_edge];
#   % furthermore, each element is subdivided into one edge and two elements
#   % we store this information in refT2ET
#   %   [ indices of original elements which were refined; ...
#   %     indices of the bisection edges;  ...
#   %     indices of the left child elements (=1st row); ...
#   %     indices of the right child elements ]
#   refT2ET = [ ...
#     [afem.aET(1,iEref); iEnew_left; afem.aET(1,iEref); iTnew_left], ...
#     [afem.aET(2,iEref1); iEnew_right; afem.aET(2,iEref1); iTnew_right] ];

#   % Projection is done by the FEM datastructure.
#   afem = afem.fem.interpolate(afem, afemnew, refE2PE, refT2ET);
# else
#   afem = afemnew;
# end

# if ~isempty(afem.fixbdry)
#   afem = afem.fixbdry(afem);
# end

# end





# ======================================================================
# below is some additional functionality that we 
# don't need for now. 




# %
# % helper function for creating new afem objects from old ones.
# % links all the functions but doesn't create the actual mesh.
# %
# function afem = empty_afemm2(fem, fixbdry)

# % store the finite element data structure
# afem.fem = fem;

# % utility functions
# afem.plotMesh = @plotMesh;
# afem.refine = @refine;
# afem.refineUniform = @refineUniform;
# afem.fixbdry = fixbdry;
# afem.coarsenAll = @coarsenAll;
# afem.coarsenNodes = @coarsenNodes;
# afem.getCoarseNodes = @getCoarseNodes;
# afem.randomRefine = @randomRefine;
# afem.midpoints = @midpoints;
# afem.dorfler_marking = @dorfler_marking;
# afem.dorfler_marking_T = @dorfler_marking_T;
# afem.mark_biggest = @mark_biggest;
# afem.mark_T = @mark_T;

# end




# %
# % afem = coarsenAll(afem)
# % coarsens all nodes which can be coarsened.
# %
# function afem = coarsenAll(afem)

# % for safety check whether there are any vertices at all
# if afem.nP == afem.nP0
#   return;
# end

# % get all free nodes
# iP = getCoarseNodes(afem);

# % coarsen them
# afem = coarsenNodes(afem, iP);

# end



# function afem = mark_T(afem, iMark)
# for i = 1:3
#   afem.mark(afem.aTE(i, iMark)) = 1;
# end
# end

# %
# % return midpoints of elements
# %
# function Q = midpoints(afem)
# Q = ( afem.P(:, afem.T(1, :)) + afem.P(:, afem.T(2,:)) ...
#   + afem.P(:, afem.T(3, :)) ) / 3;
# end


# %
# % Dorfler Marking: takes list of edge indicators
# %    and marks a minimal number of edges such that
# %    total \sum_{marked} eta \geq \theta \sum_{all} eta
# %
# % TODO: rename to dorfler_marking_E
# %
# function afem = dorfler_marking(afem, eta, theta)
# % sort indicators
# [eta, I] = sort(eta, 'descend');
# % compute total sum of indicators
# sum_all = sum(eta);
# % in a loop find out how many we need to take.
# sum_marked = 0; n = 0;
# while sum_marked < theta * sum_all
#   n = n + 1;
#   sum_marked = sum_marked + eta(n);
# end
# % mark edges with the n largest indicators
# afem.mark(I(1:n)) = 1;
# end

# %
# % Dorfler Marking of Elements instead of edges
# %
# function afem = dorfler_marking_T(afem, eta, theta)
# [eta, I] = sort(eta, 'descend');
# sum_eta = sum(eta);
# sum_marked = 0; n = 0;
# while sum_marked < theta * sum_eta
#   n = n + 1;
#   sum_marked = sum_marked + eta(n);
# end
# iM = I(1:n);
# afem = mark_T(afem, iM);
# end

# %
# % Mark Biggest: marks the theta-fraction of the edges with
# %    biggest indicators (and mark at least 10)
# %
# function afem = mark_biggest(afem, eta, theta)
# [eta, I] = sort(eta, 'descend');
# n = min(afem.nE, max(10, ceil(theta * afem.nE)));
# afem.mark(I(1:n)) = 1;
# end

# %
# % Plot just the mesh itself.
# %
# function plotMesh(afem, col)
# if nargin < 2
#   col = 'b';
# end
# trisurf(afem.T', afem.P(1, :)', afem.P(2, :)', 0*afem.P(1,:)', ...
#   'facecolor', 'none', 'edgecolor', col);
# view(2);
# end


# %
# % function J = getCoarseNodes(afem)
# % returns all the nodes which are free for coarsening
# %
# function J = getCoarseNodes(afem)

# % for safety check whether there are any vertices at all
# if afem.nP == afem.nP0
#   J = [];
#   return;
# end

# %% step 1: find the vertices which can be refined
# % valence of vertices (in how many elements do they appear?)
# valence = accumarray(afem.T(:), ones(3*afem.nT, 1), [afem.nP, 1]);
# % find the nodes which have valence 2 or 4 and which were not
# % in the original mesh!
# J = unique(afem.T(1, afem.T(1, :) > afem.nP0));
# J = J( ((valence(J) == 2) | (valence(J) == 4)) );

# end



# %
# % node coarsening routine
# %
# function afemnew = coarsenNodes(afem, iP)

# % for safety check whether there are any vertices at all
# if afem.nP == afem.nP0
#   afemnew = afem;
#   return;
# end

# %% check also whether all nodes iP CAN be coarsened! (TODO)
# iP0 = getCoarseNodes(afem);
# iP = intersect(iP, iP0);

# % copy the fem datastructure
# afemnew = afem;

# %% preliminary version of coarsening: loop through all options ...
# marker = zeros(afemnew.nP, 1); marker(iP) = 1;
# for iT = 1:afemnew.nT
#   p0 = afemnew.T(1, iT);
#   if p0 > 0
#     if marker(p0) == 1
#       % the left child comes first so here is the right "neighbour"
#       neig = afemnew.aET(2, afemnew.aTE(2, iT));
#       % restore the parent element
#       afemnew.T(:, iT) = [afem.T(3, iT); afem.T(3, neig); afem.T(2, iT)];
#       % kill the neighbour
#       afemnew.T(1, neig) = 0;
#     end
#   end
# end

# %% clean up T and P arrays
# % delete the removed elements from the mesh.
# afemnew.T(:, afemnew.T(1, :) == 0) = [];
# afemnew.nT = size(afemnew.T, 2);
# % make an array of P indices, giving each vertex in the
# % old mesh a new index in the new mesh
# indP = ones(afemnew.nP, 1); indP(iP) = 0; indP = cumsum(indP);
# afemnew.T = indP(afemnew.T);
# % remove the coarsened vertices
# afemnew.P(:, iP) = [];
# afemnew.nP = size(afemnew.P, 2);

# %% restriction of the FE DOFs
# afemnew = afem.fem.restrict(afemnew, afem, iP);

# %% update edges and adjacency arrays...
# % do the cheap version ... !!!! can't do this for edge-based DOFs !!!!
# afemnew = edges_init(afemnew);
# afemnew.mark = zeros(1, afemnew.nE);

# end



# function dorfler_marking2_T()
# dummy = 1;
# end

