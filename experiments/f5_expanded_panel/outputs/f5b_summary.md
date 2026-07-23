# F5b — s-law on the expanded panel

| corpus | V | s | s/V | rmse |
|---|---:|---:|---:|---:|
| lang_swedish | 7906 | 35429.9 | 4.4814 | 0.0832 |
| lang_portuguese | 8840 | 41551.8 | 4.7004 | 0.0804 |
| lang_german | 10898 | 95.8 | 0.0088 | 0.1024 |
| lang_polish | 11730 | 37407.5 | 3.1890 | 0.0788 |
| lang_italian | 17050 | 159.1 | 0.0093 | 0.1070 |
| lang_finnish | 21795 | 164.4 | 0.0075 | 0.1032 |
| brown | 42256 | 514.5 | 0.0122 | 0.1113 |
| wikitext_1M | 45905 | 549.6 | 0.0120 | 0.1102 |
| cornell_dialogs | 52280 | 630.4 | 0.0121 | 0.1144 |
| census_surnames | 162254 | 4316.6 | 0.0266 | 0.0094 |

- panel (language corpora) median s/V: 0.0121
- combined regression (25 English f2 + 9 panel): beta=0.6085 95% CI [-0.3032,1.5202] R2=0.0508