# f20 — lambda as fingerprint axis / depth invariance

Reference: language free-fit lambda median 20.6 (EN 25, f3), 21.9 (7 languages, f19).

## A. Simulated mechanisms (free lambda)
- ga_reference: median lam 24.1 (n=4, range [22.6, 24.8])
- pure_simon_decay: median lam 24.0 (n=4, range [22.2, 25.6])
- classic_simon: median lam 26.4 (n=4, range [23.9, 28.3])

## B. 12-language panel concats
- median lam 26.4, IQR [25.1, 26.8], range [24.0, 35.5]
- corr(lam, tok/V) across panel: -0.240
- lambda-ZM improves fit on 11/11

## C. Lambda vs depth (within-corpus thinning)
- Complete Works of Shakespeare: tok/V 40: lam 23.5  tok/V 26: lam 21.8  tok/V 26: lam 21.9  tok/V 18: lam 20.5  tok/V 17: lam 20.1  tok/V 12: lam 18.3  tok/V 12: lam 18.6  tok/V 8: lam 16.1  tok/V 8: lam 16.3
- War and Peace: tok/V 33: lam 23.7  tok/V 21: lam 21.6  tok/V 21: lam 21.7  tok/V 14: lam 19.0  tok/V 14: lam 19.3  tok/V 10: lam 17.4  tok/V 10: lam 17.3  tok/V 7: lam 16.2  tok/V 7: lam 16.0
- Les Miserables: tok/V 25: lam 23.7  tok/V 17: lam 22.9  tok/V 17: lam 22.6  tok/V 11: lam 21.2  tok/V 11: lam 21.3  tok/V 8: lam 19.7  tok/V 8: lam 19.8  tok/V 6: lam 18.6  tok/V 6: lam 19.0
- King James Bible: tok/V 63: lam 23.8  tok/V 40: lam 21.6  tok/V 40: lam 21.5  tok/V 26: lam 18.6  tok/V 26: lam 18.8  tok/V 17: lam 16.1  tok/V 17: lam 15.7  tok/V 12: lam 14.2  tok/V 12: lam 14.6
- pooled corr(lam, tok/V): +0.629