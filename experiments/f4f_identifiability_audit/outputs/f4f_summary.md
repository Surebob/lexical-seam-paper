# F4f — identifiability audit

## A. Model-free dial: corr(sd(logN|N>=5), f8 sd_T) = -0.530
- FI Seitseman veljesta: sd_logcount_ge5=0.853 (f8 sd_T=2.076)
- ZH Three Kingdoms: sd_logcount_ge5=0.876 (f8 sd_T=2.478)
- AR 1001 Nights: sd_logcount_ge5=0.925 (f8 sd_T=2.371)
- RU War and Peace: sd_logcount_ge5=0.927 (f8 sd_T=1.810)
- EN Federalist: sd_logcount_ge5=0.999 (f8 sd_T=1.682)
- EN Ulysses: sd_logcount_ge5=1.007 (f8 sd_T=2.102)
- IT Promessi Sposi: sd_logcount_ge5=1.019 (f8 sd_T=2.245)
- PT Dom Casmurro: sd_logcount_ge5=1.021 (f8 sd_T=1.845)
- DE Wahlverwandtschaften: sd_logcount_ge5=1.025 (f8 sd_T=2.065)
- EN Moby Dick: sd_logcount_ge5=1.060 (f8 sd_T=1.721)
- EN War and Peace (tr): sd_logcount_ge5=1.060 (f8 sd_T=1.732)
- EN Pride and Prejudice: sd_logcount_ge5=1.103 (f8 sd_T=1.739)
- EN Shakespeare: sd_logcount_ge5=1.126 (f8 sd_T=1.940)
- SV Gosta Berling: sd_logcount_ge5=1.142 (f8 sd_T=1.939)

## B. Landscape (per corpus: solutions within dNLL<=5 of best)
- Shakespeare: 3 near-solutions; sd_T range 1.88-3.47; n_total range 109637-600659; best-NLL sd_T=1.88 (pred_err 0.488); best-pred sd_T=3.37 (pred_err 0.227)
- War and Peace: 4 near-solutions; sd_T range 2.71-3.71; n_total range 42187-78707; best-NLL sd_T=2.72 (pred_err 0.659); best-pred sd_T=3.71 (pred_err 0.402)
- Moby Dick: 6 near-solutions; sd_T range 2.17-3.23; n_total range 52621-112120; best-NLL sd_T=2.33 (pred_err 0.248); best-pred sd_T=3.23 (pred_err 0.026)
- Pride and Prejudice: 15 near-solutions; sd_T range 2.43-3.69; n_total range 20668-45192; best-NLL sd_T=2.72 (pred_err 0.831); best-pred sd_T=2.43 (pred_err 0.428)
- Federalist Papers: 8 near-solutions; sd_T range 2.68-4.01; n_total range 23722-71435; best-NLL sd_T=4.01 (pred_err 0.738); best-pred sd_T=3.59 (pred_err 0.118)

## C. Joint multi-depth fits
- EN Shakespeare@65k: joint sd_T=3.056, near-solution spread 0.044
- EN Moby@65k: joint sd_T=3.840, near-solution spread 1.914
- EN P&P@65k: joint sd_T=2.571, near-solution spread 0.051
- EN Ulysses@65k: joint sd_T=2.769, near-solution spread 1.580
- RU WarPeace@65k: joint sd_T=1.955, near-solution spread 1.459
- ZH ThreeKingdoms@65k: joint sd_T=3.311, near-solution spread 1.186
- FI Seitseman@65k: joint sd_T=2.533, near-solution spread 0.164
- IT Promessi@65k: joint sd_T=2.046, near-solution spread 1.142