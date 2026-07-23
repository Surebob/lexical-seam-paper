# F14b — nulls: what pins the 1.2%?

| variant | run | M | V | s | s/V |
|---|---:|---:|---:|---:|---:|
| classic_simon | 0 | 0.4M | 20011 | 201.8 | 0.01008 |
| classic_simon | 0 | 1.6M | 80114 | 811.0 | 0.01012 |
| classic_simon | 1 | 0.4M | 20126 | 208.2 | 0.01034 |
| classic_simon | 1 | 1.6M | 80237 | 814.9 | 0.01016 |
| ga_reference | 0 | 0.4M | 25825 | 305.3 | 0.01182 |
| ga_reference | 0 | 1.6M | 58938 | 730.3 | 0.01239 |
| ga_reference | 1 | 0.4M | 25942 | 311.2 | 0.012 |
| ga_reference | 1 | 1.6M | 58991 | 729.5 | 0.01237 |
| pure_simon_decay | 0 | 0.4M | 20985 | 240.3 | 0.01145 |
| pure_simon_decay | 0 | 1.6M | 52837 | 634.1 | 0.012 |
| pure_simon_decay | 1 | 0.4M | 20999 | 245.5 | 0.01169 |
| pure_simon_decay | 1 | 1.6M | 52640 | 638.2 | 0.01212 |

- ga_reference: s/V 0.0118-0.0124 (median 0.0122)
- pure_simon_decay: s/V 0.0115-0.0121 (median 0.0118)
- classic_simon: s/V 0.0101-0.0103 (median 0.0101)