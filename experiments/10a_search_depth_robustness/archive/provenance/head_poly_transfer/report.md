# Zipf Head Polynomial Transfer

- Zero-shot test: fit the degree-5 head polynomial on one corpus, keep those coefficients fixed, and evaluate on the other corpus with that corpus's own ZM baseline plus the universal step-2 term.

## In-Domain Reference

### Shakespeare
- ZM RMSE: `0.184430425006`
- Step-2 RMSE: `0.182851543528`
- Step-2 + in-domain poly5 RMSE: `0.150549619124`
- Step-10 monster RMSE: `0.159871852997`

### War and Peace
- ZM RMSE: `0.182982301215`
- Step-2 RMSE: `0.180436822151`
- Step-2 + in-domain poly5 RMSE: `0.145812713187`
- Step-10 monster RMSE: `0.145347404652`

## Zero-Shot Transfer

### Shakespeare poly5 -> War and Peace
- Target ZM RMSE: `0.182982301215`
- Target Step-2 RMSE: `0.180436822151`
- Target in-domain poly5 RMSE: `0.145812713187`
- Target step-10 RMSE: `0.145347404652`
- Zero-shot transferred poly5 RMSE: `0.151277470943`

### War and Peace poly5 -> Shakespeare
- Target ZM RMSE: `0.184430425006`
- Target Step-2 RMSE: `0.182851543528`
- Target in-domain poly5 RMSE: `0.150549619124`
- Target step-10 RMSE: `0.159871852997`
- Zero-shot transferred poly5 RMSE: `0.155652365381`
