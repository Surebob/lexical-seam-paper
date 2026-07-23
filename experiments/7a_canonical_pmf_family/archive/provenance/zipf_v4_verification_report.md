# v4 Verification Note

## POS Exponents

- automated all-25 POS forced alpha: `0.5449219781561354`
- automated all-25 POS 95% CI: `[0.5324369101316629, 0.5574070461806079]`
- automated all-25 POS mean-alpha t-test p-value vs 0.50: `5.350442172058531e-08`
- smooth free-fit statistical alpha: `0.5214222874234132`
- smooth free-fit statistical 95% CI: `[0.4914985326409799, 0.5513460422058466]`
- smooth free-fit mean-alpha t-test p-value vs 0.50: `0.1421598294352621`

## Simulation Recovery

- overall smooth replicate exact-winner match rate: `0.484000`
- overall single-ZM replicate exact-winner match rate: `0.200000`
- overall smooth majority-winner match count: `12` / 25
- overall single-ZM majority-winner match count: `5` / 25
- high-c / IS block size: `11`
- high-c smooth exact-match mean: `1.000000`
- high-c ZM exact-match mean: `0.000000`
- low-c / exp block size: `14`
- low-c smooth exact-match mean (full winner): `0.078571`
- low-c ZM exact-match mean (full winner): `0.357143`
- low-c smooth modal-vs-empirical-top100 match rate: `0.714286`
- low-c ZM modal-vs-empirical-top100 match rate: `0.571429`
- low-c smooth modal-vs-empirical-full match rate: `0.071429`
- low-c ZM modal-vs-empirical-full match rate: `0.357143`

## Pending Artifact Status

- four-way PMF table materialized here: `table_a_fourway_pmf.csv`
- four-way winner counts: `{'softk': 10, 'zm': 11, 'moe': 4}`
- soft-k lambda distribution: `{0.0003: 7, 0.001: 5, 0.01: 2, 0.003: 4, 0.0001: 5, 0.03: 2}`
- 66-book Bible table already exists in `../zipf_angle6_bible_books/bible_books_table.csv`
- extended widened-grammar low-c / multilingual follow-up has not been run yet in a saved bundle

