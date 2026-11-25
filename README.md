# A Galerkin Finite Element Method for the Fractional Calderón Problem

This repository contains the MATLAB code accompanying the paper **A Galerkin Finite Element Method for the Fractional Calderón Problem** by **Mukul Dwivedi, Jesse Railo, and Andreas Rupp**.

The examples reconstruct the interior solution and the potential from one partial exterior measurement. They use the two-stage Galerkin method described in the paper.

## Running the examples

Extract the archive and open its main folder in MATLAB. To run all three examples, enter

```matlab
run_all
```

Each example can also be run separately from the same folder.

| Script | Example and figures |
| --- | --- |
| `examples/run_example41_1d_smooth.m` | Example 4.1. The one-dimensional Lipschitz potential, the recovered solution, the error decomposition, and the noise study. |
| `examples/run_example42.m` | Example 4.2. The quadratic reconstruction of the discontinuous potential. |
| `examples/run_example42_tv.m` | Example 4.2. The TV reconstruction and the reconstruction after the interval refit. |
| `examples/run_example43_2d.m` | Example 4.3. The true and recovered potentials, the recovered solution, the state error decomposition, and the potential noise study. |

For example,

```matlab
run('examples/run_example43_2d.m')
```

The scripts begin with `clear all; clc; close all;`. The individual drivers leave their computed results in `result`. The combined driver leaves them in `results{1}`, `results{2}`, and `results{3}`.

The code uses MATLAB without additional toolboxes. Example 4.3 forms dense matrices with 6,241 interior unknowns and needs several gigabytes of memory. Its fine reference calculation is the most expensive part of the run.

## Figures and numerical output

The figures are saved as PDFs in `output/figures`, including the TV reconstruction before the interval refit. A repeated run overwrites the corresponding figures. The two potential surfaces and the line plots use vector graphics. The two-dimensional solution surfaces use 400 dpi.

Parameters, mesh diameters, residuals, reconstruction errors, and fitted constants are printed in the MATLAB command window. The drivers do not save data files, logs, or LaTeX files. The dashed logarithmic fit and its legend use `C` rounded to one significant digit. The fitted values of `C` and `gamma` are stored at full precision in the results and printed to six significant digits in the command window.

PDF export uses the `Append` option of `exportgraphics`, available from MATLAB R2021b. The original experiments used MATLAB R2025b Update 5. If figure creation fails, the numerical results remain in memory. The figures can be regenerated without repeating the solve, for example

```matlab
fcp_plot_results(results{3},results{3}.plot_config)
```

After running an individual example, the corresponding command is `fcp_plot_results(result,result.plot_config)`. These commands use results from the current MATLAB session.

## Parameters and source files

The interior and observation meshes are constructed in `src/fcp_mesh.m` and `src/fcp_observation_mesh.m`, respectively. The mesh subdivision counts, exterior data, noise levels, and regularization parameters are specified in `src/fcp_parameters.m`. The parameter rules are evaluated in `src/fcp_parameter_rule.m`, and `src/fcp_run_experiment.m` assembles and runs each experiment.

| Setting | Example 4.1 | Example 4.2 | Example 4.3 |
| --- | --- | --- | --- |
| Fractional order | 0.6 | 0.6 | 0.5 |
| Reconstruction mesh | 800 intervals | 800 intervals | 80 by 80 squares, each split into two triangles |
| Reference mesh | 4,000 graded intervals | 4,000 graded intervals | 640 by 640 squares, each split into two triangles |
| Coefficient degree | Discontinuous linear | Piecewise constant | Discontinuous linear |
| Noise study | 30 levels, 10 realizations per level | Four representative noise levels | 20 levels, 15 realizations per level |

The synthetic data are computed on the reference mesh. The exact potential is used to generate these data and evaluate errors. The reconstruction uses the noisy data and the prescribed parameter rules. Noise is added directly to the preprocessed data and normalized in the observation mass norm. The seed is fixed for reproducibility.

Example 4.2 distinguishes the TV minimizer from the subsequent interval refit. The refit assumes a potential that is constant on one interval, zero elsewhere, and bounded between 0 and 1. The interval is detected by thresholding at 0.5, and its constant value is found by regularized least squares with the same coefficient penalty.

The remaining files in `src` contain matrix assembly, quadrature, forward and inverse solvers, error evaluation, and plotting. The figure settings are in `src/fcp_plot_original_style.m`. The potential surfaces use fixed page dimensions, and the other figures are exported by `src/save_figure_pdf.m`.

The command window also reports the mesh condition indicator `h^beta/sqrt(alpha)` with `beta = 0.45` and `h` equal to the largest interior, observation, or coefficient element diameter. The sufficient mesh inequality in Lemma 3.6 requires this indicator to be at most one. The present fixed-mesh noise studies use smaller regularization parameters and do not verify the joint limit in Theorem 3.7. A convergence study would also need the reference-data error to decrease. The indicator is reported without changing the reconstruction parameters.
