# f15b — cross-language width-law panel (full depth; reconstructed from run log)

| corpus | lang | tokens | V | tok/V | s/V | reps |
|---|---|---:|---:|---:|---:|---:|
| da_pg2144 | da | 494,964 | 20,924 | 23.7 | **0.0118** | 2 |
| da_pg2143 | da | 170,218 | 9,466 | 18.0 | **0.0117** | 2 |
| de_pg40739 | de | 211,545 | 23,871 | 8.9 | **0.0093** | 2 |
| de_pg65662 | de | 159,570 | 24,345 | 6.6 | **0.0080** | 2 |
| de_pg69944 | de | 159,548 | 21,977 | 7.3 | **0.0090** | 2 |
| el_pg33396 | el | 156,937 | 25,868 | 6.1 | **0.0082** | 2 |
| es_pg2000 | es | 383,633 | 22,943 | 16.7 | **0.0104** | 2 |
| es_pg78068 | es | 264,719 | 23,452 | 11.3 | **0.0102** | 2 |
| es_pg61851 | es | 177,773 | 14,582 | 12.2 | **0.0103** | 2 |
| fi_pg70428 | fi | 267,610 | 52,228 | 5.1 | **0.0081** | 2 |
| fi_pg69850 | fi | 186,984 | 40,282 | 4.6 | **0.0079** | 2 |
| fi_pg67616 | fi | 158,848 | 28,815 | 5.5 | **0.0079** | 2 |
| fr_pg25335 | fr | 282,117 | 20,668 | 13.6 | **0.0109** | 2 |
| fr_pg13951 | fr | 239,953 | 13,950 | 17.2 | **0.0110** | 2 |
| fr_pg14287 | fr | 209,844 | 14,224 | 14.8 | **0.0111** | 2 |
| hu_pg66903 | hu | 217,243 | 37,754 | 5.8 | **0.0074** | 2 |
| hu_pg41504 | hu | 158,451 | 31,555 | 5.0 | **0.0077** | 2 |
| it_pg65391 | it | 444,080 | 25,423 | 17.5 | **0.0109** | 2 |
| it_pg78619 | it | 202,263 | 23,491 | 8.6 | **0.0094** | 2 |
| it_pg26961 | it | 166,146 | 10,890 | 15.3 | **0.0096** | 2 |
| nl_pg28120 | nl | 196,197 | 16,537 | 11.9 | **0.0109** | 2 |
| nl_pg22968 | nl | 184,600 | 26,576 | 6.9 | **0.0084** | 2 |
| nl_pg23759 | nl | 173,560 | 12,786 | 13.6 | **0.0115** | 2 |
| pt_pg31552 | pt | 1,767,544 | 179,958 | 9.8 | **0.0090** | 2 |
| pt_pg62383 | pt | 797,817 | 26,962 | 29.6 | **0.0121** | 2 |
| pt_pg40409 | pt | 221,439 | 20,056 | 11.0 | **0.0107** | 2 |
| sv_pg2100 | sv | 736,812 | 24,297 | 30.3 | **0.0121** | 2 |
| sv_pg57357 | sv | 183,205 | 9,915 | 18.5 | **0.0116** | 2 |
| sv_pg13100 | sv | 152,634 | 18,085 | 8.4 | **0.0099** | 2 |

- da: n=2, median s/V 0.0118, median tok/V 20.8
- de: n=3, median s/V 0.0090, median tok/V 7.3
- el: n=1, median s/V 0.0082, median tok/V 6.1
- es: n=3, median s/V 0.0103, median tok/V 12.2
- fi: n=3, median s/V 0.0079, median tok/V 5.1
- fr: n=3, median s/V 0.0110, median tok/V 14.8
- hu: n=2, median s/V 0.0076, median tok/V 5.4
- it: n=3, median s/V 0.0096, median tok/V 15.3
- nl: n=3, median s/V 0.0109, median tok/V 11.9
- pt: n=3, median s/V 0.0107, median tok/V 11.0
- sv: n=3, median s/V 0.0116, median tok/V 18.5

- panel: 29 corpora, 11 languages; median s/V 0.0102
- corr(log tokens-per-type, s/V) = +0.938
- deep corpora (tok/V >= 12): n=13, median s/V 0.0111  |  shallow (tok/V < 12): n=16, median 0.0087
- EN reference 0.0118 (median tok/V ~ 19). One corpus (pt_pg31552) has a single replicate (process died before its second).