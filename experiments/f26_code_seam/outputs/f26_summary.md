# f26 — the code seam

Reference fingerprints: language s/V 0.0118 (lam ~20-26 by depth), surnames 0.0266, belt 0.0166, classic Simon 0.0101.

| corpus | tokens | V | tok/V | b | c | lambda | s/V |
|---|---:|---:|---:|---:|---:|---:|---:|
| python_sitepackages | 22,355,917 | 513,736 | 43.5 | 1.2738 | 37.86 | 54.51 | 0.0132 |
| js_projects | 430,553 | 18,178 | 23.7 | 1.4651 | 175.28 | 27.84 | 0.50482 |

Caveats: identifier tokenization (no snake/camel splitting — a raw 'code word' definition); site-packages mixes library code + vendored assets; depth differs from book regime — compare via the depth curves, not raw constants.