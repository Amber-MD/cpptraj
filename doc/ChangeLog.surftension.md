# Changelog: `surftension` Action

Branch: `action-surface_tension`
Started: 2026-08-31
Target version: V7.11.0
Author: Nathan D Levinzon <ndlevinzon@gmail.com>

Working log for the capillary-wave surface-tension Action.

---

## ChangeLog.v7.md (landed)

New command `surftension` added under New Commands. Internal version `V7.11.0`.

```
surftension [<name>] <mask> [mask2 <mask2>] temp <T>
  [normal {x|y|z}] [nsurf {1|2}] [side {upper|lower}]
  [interface {willard|itim}]
  [gridspacing <d>] [dz <d> | dnormal <d>]
  [sigmaxy <d>] [sigmaz <d> | sigmanormal <d>]
  [bulkhalfwidth <d>] [threshold <frac>]
  [qmin <q>] [qmax <q>] [lx <Lx>] [ly <Ly>] [lz <Lz>]
  [nblock <frames>] [dt <ps>] [blocktime <ps>]
  [spectrumout <file>] [roughout <file>] [blockout <file>]
  [summaryout <file>]
  [spectrumagr <file>] [roughagr <file>] [blockagr <file>]
  [spectrumgnu <file>] [roughgnu <file>] [blockgnu <file>]
```

---

## Files changed

| File | Change |
|------|--------|
| `src/Action_SurfaceTension.h` | New Action class; MPI `SyncAction` / `trajComm_` |
| `src/Action_SurfaceTension.cpp` | Init / Setup / DoAction / Print / SyncAction + helpers |
| `src/Command.cpp` | `#include` and `AddCmd(..., "surftension")` |
| `src/cpptrajfiles` | `Action_SurfaceTension.cpp` in `COMMON_SOURCES` |
| `src/cpptrajheaders` | Install header for libcpptraj |
| `src/Version.h` | `V7.10.0` → `V7.11.0` |
| `doc/ChangeLog.v7.md` | New Commands entry |
| `test/Test_SurfTension/RunTest.sh` | Help + smoke test |
| `test/Makefile` | `test.surftension` target |

---

## Development log

### 2026-08-31

- Created this changelog on branch `action-surface_tension` before any Action source.
- Implemented `Action_SurfaceTension`: z-slab capillary-wave γ from Gaussian-smoothed
  density interfaces, 2D DFT height spectrum, q-shell averaging, optional `*out` files.
- Wired `surftension` in `Command.cpp` and `cpptrajfiles`.
- Bumped version to V7.11.0; documented in `ChangeLog.v7.md`.
- Added `test/Test_SurfTension` smoke test (help + 1-frame Init/Setup).
- NVT (fixed Lx, Ly). No plots.
- MPI `SyncAction`: packed `ReduceMaster` SUM of the three |h_q|² spectra
  (Radial-style). Frame counts / leftover `nblock` frames AllReduce SUM;
  nx, ny, |q|, Lx, Ly AllReduce/ReduceMaster MAX so empty ranks (zeros) do
  not clobber. Roughness and block DataSets keep default concat sync.
  `nblock` is per-rank. Print block SEM uses `block_gamma_->Size()` after
  DataSet sync, not the master's local `n_blocks_`.
- `interface {willard|itim}`: default Willard-Chandler is the existing Gaussian
  density isosurface. `itim` is per-column min/max of `<mask>`, split at Lz/2
  after circular recenter (empty half-column skips the frame). `dz` / Gaussian
  / `threshold` apply only to Willard-Chandler.
- OpenMP on the 2-D DFT (`ST_HeightPower`): flattened (kx, ky) loop with
  `schedule(dynamic)`, same style as `radial` / `rms2d` (no `collapse()`).
- OpenMP on the 3-D Gaussian filter: one parallel region, three separable
  passes; each 1-D line is an independent `omp for` with thread-local buffers.
- Grace (xmgrace) and gnuplot DataFile writers: `*out` uses the extension
  (`.agr`/`.xmgr`, `.gnu`); `spectrumagr`/`spectrumgnu` (and rough/block)
  force the format. Spectrum meshes get xlabel `q (Ang^-1)`.
- Bending modulus κ: Helfrich linear fit of `1/(q² S)` vs `q²` on `[qmin,qmax]`
  (≥ 3 shells). Reports κ in kT (and J); `kappaq` / `kappaqtop` / `kappaqbot`
  vs q; `bkappa` with `nblock`. Plateau γ is unchanged (Python).
- `normal {x|y|z}` (default z): permute Cartesian coordinates into
  (lateral 1, lateral 2, normal) before wrap / circular recenter / Willard /
  ITIM / DFT. Freeze lateral lengths only (NVT). Area = Lt1×Lt2.
  `lx`/`ly`/`lz` override Cartesian box lengths. `dnormal` ≡ `dz`,
  `sigmanormal` ≡ `sigmaz` (error if both given and they differ).
  Init/Print report e.g. `Slab normal: z (interface plane x-y)`. Warns on
  non-orthogonal boxes.
-   Height-field 2-D FFT uses cpptraj `PubFFT` (row-column 1-D FFTs), then
  divides by `nx*ny` (numpy `fft2` convention). Same `S(q)` as the old
  direct DFT. Per-block γ/κ print in `Print()`, not during the frame loop.
- Default `qmin` is the smallest fundamental wavevector from the unit-cell
  laterals, `2π/max(Lt1,Lt2)`, set on the first good frame. Both `2π/Lt1`
  and `2π/Lt2` are printed. `summaryout` writes a parseable key/value file
  (γ, κ, roughness, q window, box, frames). Skip-frame warnings are capped
  at 5. `blocktime` (ps) with `dt` (analyzed-frame spacing) sets `nblock`.
  Output directories are not created; the user must give existing paths.
- `nsurf {1|2}` (default 2, vacuum or a second phase on both sides of the
  film). `nsurf 1` with `side {upper|lower}` uses one interface. `mask2`
  is a second mask for the lower surface (leaflet / liquid–liquid); the
  upper surface then comes from `<mask>` with no mid-box split. Both
  masks share one circular recenter so the film is not split apart.

---

## Not in this change

- ITIM probe-sphere radius (min/max is the probe → 0 limit)
- Numeric regression vs a water-slab trajectory (no slab traj in the test suite yet)
