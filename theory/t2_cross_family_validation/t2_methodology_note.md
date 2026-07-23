# T2 Methodology Note

T2 uses the live Section 2.4 enumerative search implementation from `eml_zipf_enriched_search.py` with `beam_width=50`, `max_steps=4`, `keep_all_until_step=2`, semantic hashing rounded to ten decimals, `EXP_CLAMP=30`, and `VALUE_ABS_LIMIT=1e6`.

The depth-4 revision is required because the Poisson generalized-KL generator `x*log(x)-x+1` is not reachable at step 2. It first becomes operationally reachable at step 3 as the equivalent expression `sub[add[1,neg[x]],mul[neg[x],log[x]]]` when the target residual is exactly KL-shaped. Euclidean and Itakura-Saito are reachable at step 2.

For each synthetic family, the generated frequencies are sorted descending to form a rank-frequency curve, the corresponding single-family curve is fitted, and the residual `log_frequency - fitted_log_frequency` is passed to the parameter-free search on the normalized rank axis `x = 0.05 + 0.95*log(r)/log(V)`.
