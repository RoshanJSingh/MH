# MinimaHop-Stacker

A global optimization tool that finds the most stable way to stack layers of 2D materials such as COFs, graphene and MOFs. It uses the Minima Hopping algorithm, with energies calculated by xtb.

## The problem

When you place one molecular layer on top of another, there are many possible positions, and most of them are only locally stable. The energy landscape is rough and full of shallow valleys, so an ordinary optimizer tends to settle in the first one it finds.

Minima Hopping gets around this by treating the search like a physical system that can jump over energy barriers. It changes the size of its jumps (its "kinetic energy") and how willing it is to accept a worse position (its "temperature") as it goes, so it does not get stuck.

## What is being optimized

The position of the top layer relative to the bottom layer is described by 7 parameters:

| Parameter | Meaning |
|-----------|---------|
| `Tx`, `Ty`, `Tz` | shift in X and Y, and the vertical gap Z |
| `Cx`, `Cy` | the point the layer rotates around |
| `cos_like`, `sin_like` | the rotation angle, stored as two normalized parts |

## How one hop works

Each hop has three steps. A simple way to picture it is a hiker walking blindfolded through hills.

1. **Kick.** Shake the stack at random to get out of the current valley.
2. **Slide.** Nudge the stack step by step until it settles at the bottom of the new valley.
3. **Decide.** Either stay in the new valley or go back to the old one.

## How that looks in the code

The core code is in `scripts/main.py`.

### 1. The kick

To get out of a local minimum, random Gaussian noise is added to all 7 parameters. How big the kick is depends on the current kinetic energy `ke`.

```python
# from the run_minima_hopping() loop
# step size depends on the current kinetic energy (higher KE, bigger step)
step = _step_scales(base, ke, cfg.mh.ke_ref_kj, lo, hi)

# x: current coordinates
# rng.normal: the random kick
x_prop = np.clip(x + rng.normal(0.0, step), lows, highs)
```

### 2. The slide

After the kick the structure is in a high energy state. A Hooke-Jeeves pattern search brings it down to the nearest local minimum. It needs no gradients. It just tries moving each parameter a little up and a little down, and keeps any change that lowers the energy.

```python
# from pattern_search_refine()
for d in range(D):           # for each of the 7 parameters
    for sgn in (+1, -1):     # try a step up and a step down
        cand = x.copy()
        cand[d] += sgn * step[d]

        # calculate the energy
        f_cand = float(eval_fn(cand))

        # keep the change if the energy went down
        if f_cand < fx - 1e-15:
            x, fx = cand, f_cand
            improved = True
```

### 3. The decision

Once the stack has settled into a new minimum, the Metropolis rule decides whether to keep it. A lower energy is always accepted. A higher energy is sometimes accepted, depending on the temperature `T`.

```python
# from run_minima_hopping()
delta = f_ref - fx  # new energy minus old energy

# _metropolis returns True if the move should be accepted
acc = _metropolis(delta, T, rng)

if acc:
    # accepted: raise KE a little to explore further
    ke = min(cfg.mh.ke_max_kj, max(cfg.mh.ke_min_kj, ke * cfg.mh.beta_up))
else:
    # rejected: lower KE and try a smaller jump next time
    ke = min(cfg.mh.ke_max_kj, max(cfg.mh.ke_min_kj, ke * cfg.mh.alpha_down))
```

## Installation

You need:

- Python 3.8 or newer
- NumPy (`pip install numpy`)
- [xtb](https://github.com/grimme-lab/xtb), installed and available on your `PATH`. All energy calculations go through it.

Then clone the repo:

```bash
git clone https://github.com/RoshanJSingh/MH.git
cd MH
```

## Usage

Put your single-layer (monomer) structure in an `.xyz` file. `BTA.xyz` in the repo is an example.

Create a starting config for it:

```bash
python scripts/run_mh.py --make-default --xyz BTA.xyz
```

This writes `BTA_config.json`. Edit it if you want to change the settings, then run:

```bash
python scripts/run_mh.py --config BTA_config.json
```

`input.json` in the repo is another example config. Set its `input_xyz` to your own file before using it.

When the run finishes it prints the best energy, the number of distinct minima found and the output folder. Results from eight earlier runs on BTA are in `bta_runs/`.

## Outputs

Each run writes a folder (`mh_out` by default) with:

- `local_minima.json`: every distinct stable structure found, ranked, with its energy and parameters
- `optimization_results.txt`: a summary of the best structure
- `energy_log.csv`: energy, temperature and kinetic energy at every hop
- `swarm_trace.csv`: the coordinates of every attempt
- `mh_localmin_XX_3layers.xyz`: 3D structures of the best minima, which you can open in VESTA, Avogadro or Ovito

## Tuning

The behaviour of the search is set in the `MHConfig` class:

```python
@dataclass
class MHConfig:
    n_hops: int = 20          # how many hops to try
    ke_start_kj: float = 2.3  # starting kick size
    alpha_down: float = 0.9   # how much to shrink the kick after a rejection
    beta_up: float = 1.15     # how much to grow the kick after an acceptance
```
