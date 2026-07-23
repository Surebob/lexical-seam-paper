# Angle 2: Structural Predictors of lambda_k

- Best soft-k `lambda_k` values were taken from the existing soft-k sweep.
- Structural metadata are coarse hand-coded estimates meant to test whether lambda varies with corpus heterogeneity at all.

## Correlations With log10(lambda_k)

- structure code: Pearson `-0.073491`, Spearman `-0.023816`
- log unit count: Pearson `-0.087508`, Spearman `-0.053956`
- log author count: Pearson `-0.267295`, Spearman `-0.174043`
- log era span: Pearson `-0.100705`, Spearman `-0.063112`
- heterogeneity score: Pearson `-0.157873`, Spearman `-0.046349`
- token count: Pearson `-0.011303`, Spearman `-0.058536`
- vocabulary size: Pearson `0.028974`, Spearman `0.049501`

## By Structure

- single_work: n=`15`, median lambda_k=`0.0003`
- framed_multi_unit: n=`1`, median lambda_k=`0.001`
- single_author_collection: n=`4`, median lambda_k=`0.002`
- multi_source_collection: n=`3`, median lambda_k=`0.0003`
- multi_author_composite: n=`2`, median lambda_k=`0.00165`

## Highest lambda_k

- War and Peace: lambda_k=`0.03`, structure=`single_work`, units=`1`, authors=`1`, span=`4`
- Les Miserables: lambda_k=`0.03`, structure=`single_work`, units=`5`, authors=`1`, span=`17`
- Emile: lambda_k=`0.01`, structure=`single_work`, units=`5`, authors=`1`, span=`1`
- Collected Poe: lambda_k=`0.01`, structure=`single_author_collection`, units=`70`, authors=`1`, span=`20`
- Federalist Papers: lambda_k=`0.003`, structure=`multi_author_composite`, units=`85`, authors=`3`, span=`1`

## Lowest lambda_k

- Complete Works of Shakespeare: lambda_k=`0.0001`, structure=`single_author_collection`, units=`37`, authors=`1`, span=`24`
- Moby Dick: lambda_k=`0.0001`, structure=`single_work`, units=`1`, authors=`1`, span=`1`
- Grimm's Fairy Tales: lambda_k=`0.0001`, structure=`multi_source_collection`, units=`200`, authors=`50`, span=`45`
- Wealth of Nations: lambda_k=`0.0001`, structure=`single_work`, units=`5`, authors=`1`, span=`1`
- Ulysses: lambda_k=`0.0001`, structure=`single_work`, units=`18`, authors=`1`, span=`8`

## Bible Decomposition Context

- whole-Bible failure is structurally real: per-book soft-k beats MOE on `45` / `66` books, and median per-book soft-k minus MOE is `-0.002803782823`.
- aggregate per-book held-out avg NLL is `5.604219509299` vs whole-Bible single-fit `6.016015736045`.
