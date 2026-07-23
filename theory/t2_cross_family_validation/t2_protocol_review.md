# T2 Protocol Review

This review separates three questions that the first depth-4 aggregate can conflate: data/fitting space, canonical depth-2 behavior, and whether depth-4 winners are interpretable or polynomial-like overfits.

## Live Search Semantics

- Generated candidates are sorted by RMSE for reporting.
- The diversity beam is the selected vocabulary carried to the next step; it is not the same as the top-50-by-RMSE list.
- The manuscript's step-2 winner language corresponds to the lowest-RMSE step-2 candidate, while diversity selection only controls what expressions remain available for deeper search.

## Synthetic Data/Fitting Space

### T2a Gaussian -> Euclidean

- generated data space: sorted rank-frequency curve formed from a two-component `Gaussian` density/PMF evaluated on a fixed support grid.
- fitted model space: `nonlinear_least_squares_on_sorted_log_density`.
- SR residual: `log_frequency - fitted_log_frequency` on the sorted rank-frequency curve.
- fitted parameters: `{'log_amp': 0.6189920211097156, 'mu': 0.002866994315395568, 'sigma': 1.0339335610784388}`.
- bounds hit: `{'mu': False, 'sigma': False}`.
- mixture curvature diagnostic: `{'component_corr': -0.2403354910528297, 'component_relative_l2': 1.161570620761769, 'max_abs_second_difference': 0.029259591116467476, 'mean_abs_second_difference': 0.0031154425117400778, 'sign_changes_second_difference': 859}`.

### T2b Poisson -> generalized KL

- generated data space: sorted rank-frequency curve formed from a two-component `Poisson` density/PMF evaluated on a fixed support grid.
- fitted model space: `nonlinear_least_squares_on_sorted_log_pmf`.
- SR residual: `log_frequency - fitted_log_frequency` on the sorted rank-frequency curve.
- fitted parameters: `{'log_amp': 1.4991777214494428, 'lambda': 23.123137374329787, 'mixture_mean_mle_lambda': 15.0}`.
- bounds hit: `{'lambda': False}`.
- mixture curvature diagnostic: `{'component_corr': -0.15567069042850176, 'component_relative_l2': 1.2216381573617072, 'max_abs_second_difference': 0.32094750002199923, 'mean_abs_second_difference': 0.05297847643773922, 'sign_changes_second_difference': 26}`.

### T2c Gamma -> Itakura-Saito

- generated data space: sorted rank-frequency curve formed from a two-component `Gamma` density/PMF evaluated on a fixed support grid.
- fitted model space: `expected_mle_on_unsorted_gamma_grid_plus_log_amplitude_alignment`.
- SR residual: `log_frequency - fitted_log_frequency` on the sorted rank-frequency curve.
- fitted parameters: `{'log_amp': -1.2723647516023069, 'shape': 0.8727198490161334, 'rate': 0.414261567520883}`.
- bounds hit: `{'shape': False, 'rate': False}`.
- mixture curvature diagnostic: `{'component_corr': 0.023421708096459083, 'component_relative_l2': 1.0363898155512987, 'max_abs_second_difference': 0.04984924305706989, 'mean_abs_second_difference': 0.0001225109035905885, 'sign_changes_second_difference': 1}`.

