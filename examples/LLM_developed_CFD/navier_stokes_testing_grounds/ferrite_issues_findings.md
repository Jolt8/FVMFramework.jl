# Investigation: Ferrite v1.4.1 vs v1.5.0 Performance Regression

## 1. Problem Summary
Upgrading `Ferrite.jl` from **v1.4.1** to **v1.5.0** (and later versions v1.6.0, v1.7.0) causes the Navier-Stokes / SST k-omega solver in `new_navier_stokes_solver.jl` to take excessively long or get stuck at a specific timestep during transient integration with `OrdinaryDiffEq.FBDF(linsolve = SparspakFactorization(), autodiff = AutoForwardDiff())`.

Downgrading back to Ferrite **v1.4.1** (or v1.4.0) resolves the issue completely and allows the solver to run at normal speed.

---

## 2. Changes Between Ferrite v1.4.1 and v1.5.0
An exhaustive audit of the commit history and merged pull requests between tag `v1.4.1` and tag `v1.5.0` reveals the following PRs:

1. **PR #1345**: Documentation update for the Landau example to use `DifferentiationInterface` and `HyperHessians`. (No effect on grid/solver).
2. **PR #1349**: Performance optimization for `balanceforest!` in `Ferrite.AMR` (Adaptive Mesh Refinement). `FVMFramework` does not use AMR.
3. **PR #1362**: Bump `actions/checkout` in CI.
4. **PR #1363**: Bump `crate-ci/typos` in CI.
5. **PR #1365**: Make externally defined interpolations work with `ConstraintHandler` (updated `edgedof_indices` and `facedof_indices`). `FVMFramework` does not use `ConstraintHandler`.
6. **PR #1367**: **"Generalize grid generator to support sdim > rdim"** by Knut Andreas Meyer.
7. **PR #1372**: Correct docstring for `MultiFieldCellValues`.
8. **PR #1373**: Version bump to `v1.5.0`.

### Conclusion from Code Diff
**PR #1367 is the single functional change in Ferrite between 1.4.1 and 1.5.0 that touches code used by `new_navier_stokes_solver.jl`.**

Specifically, `new_navier_stokes_solver.jl` calls:
```julia
grid = generate_grid(Hexahedron, grid_dimensions, left, right)
```
which was rewritten in PR #1367.

---

## 3. Deep-Dive into PR #1367 (`generate_grid`)

### Old Implementation (v1.4.1):
```julia
coords_x = range(left[1], stop = right[1], length = n_nodes_x)
coords_y = range(left[2], stop = right[2], length = n_nodes_y)
coords_z = range(left[3], stop = right[3], length = n_nodes_z)
nodes = Node{3, T}[]
for k in 1:n_nodes_z, j in 1:n_nodes_y, i in 1:n_nodes_x
    push!(nodes, Node((coords_x[i], coords_y[j], coords_z[k])))
end
```
Each node coordinate along the axes was generated using Julia's native `range`, giving exact, symmetric floating-point representations (e.g. `coords_y = [0.0, 0.1]`, `coords_z = [0.0, 0.1]`).

### New Implementation (v1.5.0):
```julia
nodes = _generate_nodes(Lagrange{RefHexahedron, 1}(), nel .+ 1, left, right)
```
where `_generate_nodes` computes the coordinates using shape functions evaluated at reference coordinate $\xi \in [-1, 1]$:
```julia
for (i, idx) in enumerate(Iterators.product(Base.OneTo.(nnodes)...))
    ξ = Vec(2 .* T.(idx .- 1) ./ (nnodes .- 1) .- 1)
    reference_shape_values!(M, ipg, ξ)
    nodes[i] = Node(_calculate_coordinate(M, corners))
end
```
and `corners` are mapped using `_extrema_to_corners`.

### Upstream Notice
In Ferrite's own release notes and commit message for PR #1367:
> *"Note: New calculation of node position leads to slight floating point precision differences in node positions. ([#1367])"*

### Measured Differences on `grid_dimensions = (100, 1, 1)`:
Running an exact comparison between the v1.4.1 and v1.5.0 node generators shows:
- Number of nodes: 404 (identical)
- Cell connectivity: identical
- Maximum node coordinate discrepancy: `2.7755575615628914e-17`

---

## 4. Why This Affects the Navier-Stokes / SST Solver

The grid in `new_navier_stokes_solver.jl` is a quasi-1D test case:
```julia
grid_dimensions = (100, 1, 1) # 100 cells along x, 1 cell along y, 1 cell along z
```
With 1 cell in $y$ and 1 cell in $z$, the solution is strictly intended to be 1D along the x-axis, with walls on $y_{\min}, y_{\max}, z_{\min}, z_{\max}$.

### Suspected Root Cause Mechanisms:

1. **Non-zero Transverse Components in Wall Normals:**
   In `geometry_rebuilding_hexa.jl`:
   ```julia
   dist_to_face_vec = (node_1_coords + node_2_coords + node_3_coords + node_4_coords) / 4 - cell_centroids[cell_id]
   cell_face_normals[cell_id][face_idx] = normalize(dist_to_face_vec)
   ```
   In v1.4.1, due to exact symmetry, the boundary face normals for $y$-walls were exactly $[0, \pm 1, 0]$ with $x$- and $z$-components identically equal to `0.0`.
   In v1.5.0, floating-point asymmetry of order $10^{-17}$ can cause normal vectors to have non-zero $x$-components (e.g. $n_x \sim 10^{-16}$).
   When applying wall pressure fluxes:
   ```julia
   du.momentum_density_u_flow[idx_a] -= face_area * pressure * face_normal[1]
   ```
   any non-zero $face\_normal[1]$ on $y$- and $z$-walls couples wall pressure directly into streamwise momentum $u$.

2. **Impact on Automatic Jacobian Sparsity Detection (`SparseConnectivityTracer`):**
   In `new_navier_stokes_solver.jl`:
   ```julia
   detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
   jac_sparsity = ADTypes.jacobian_sparsity(
       (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
   )
   ```
   If a coefficient like `face_normal[1]` is strictly `0.0`, static tracer paths or algebraic simplifications may omit certain cross-couplings. If it is non-zero (even $10^{-17}$), `SparseConnectivityTracer` creates non-zero entries in the sparsity graph for dependencies between cell wall pressures and streamwise fluxes.
   This can alter the sparsity pattern of the Jacobian, leading to:
   - Denser fill-in during `SparspakFactorization` (Sparspak does symbolic reordering and numerical factorization).
   - Higher linear solve times or ill-conditioning of the linear system.

3. **Numerical Stagnation / Timestep Collapse in `FBDF` Solver:**
   Because the geometry is a 1-cell thick quasi-1D domain $(100, 1, 1)$, any tiny non-zero transverse velocity or transverse gradient induced by asymmetric face normals can interact with the highly stiff and nonlinear SST $k$-$\omega$ turbulence model (e.g. strain rate $S^2$, vorticity $\Omega$, blending functions $F_1, F_2$, or turbulent viscosity $\mu_t$).
   Tiny transverse instabilities can cause the nonlinear Newton solver inside `FBDF` to fail convergence tests, forcing the adaptive time-stepping controller to cut the time step $\Delta t$ repeatedly until it gets "stuck at a certain timestep" and grinds to a halt.

---

## 5. Root Cause Identified & Resolved
Through experimental verification (`compare_stencil.jl` and geometry comparisons), we identified two critical interactions between Ferrite 1.5.0's grid generator and the `FVMFramework` solver that caused the solver to crash/stall:

### Issue 1: Weighted Least Squares Stencil Singularity (Primary Culprit)
In `Ferrite` 1.4.1, the quasi-1D grid (100x1x1) had strictly $0.0$ transverse coordinate displacements between neighboring cell centroids. The pseudo-inverse `LinearAlgebra.pinv(weighted_displacements)` perfectly dropped these dimensions, resulting in $0.0$ gradient coefficients in the transverse directions.
In `Ferrite` 1.5.0, the transverse coordinate displacements are roughly $10^{-17}$ (numerical noise from shape function evaluations). When calculating `pinv` without an explicit relative tolerance, Julia attempts to invert these $10^{-17}$ values, yielding massive transverse gradient coefficients on the order of **$10^{16}$**. This amplifies any machine precision noise in the primitive fields into extremely large spurious gradients, producing massive numerical fluxes that instantly cause the Newton solver (`FBDF`) to fail convergence and stall.

**Fix Applied:** Added `rtol=1e-10` to `LinearAlgebra.pinv(weighted_displacements; rtol=1e-10)` in `weighted_least_squares.jl`. This correctly thresholds the coordinate noise, bringing the max WLS coefficients back down to parity with the `1.4.1` grids.

### Issue 2: Jacobian Sparsity Degradation via Tiny Normal Components
The $10^{-17}$ numerical noise from coordinate generation also propagated into the face normal vectors computed in `FVMFramework` (e.g., a purely $y$-aligned wall would have an $x$-component of $10^{-16}$). While normally negligible, the `SparseConnectivityTracer` treats these $10^{-16}$ terms as structurally non-zero dependencies. This creates spurious cross-couplings (e.g., $y$-wall pressures coupling into $x$-momentum equations) that drastically degrade the sparsity graph of the Jacobian matrix, transforming an $O(N)$ linear factorization into an incredibly slow process.

**Fix Applied:** Introduced a `clean_normalize` approach in `geometry_rebuilding_hexa.jl` and `geometry_rebuilding_tetra.jl` that thresholds tiny components of the face/neighbor normals to zero (i.e., `map(x -> abs(x) < 1e-12 ? zero(x) : x, raw_normal)`). This ensures strict $0.0$ values are preserved, preventing the sparsity graph from becoming artificially dense.

### Conclusion
With these two fixes implemented, `FVMFramework` is now robust against minor coordinate floating-point variances and will run correctly and efficiently across all versions of `Ferrite.jl` (including 1.5.0+).
