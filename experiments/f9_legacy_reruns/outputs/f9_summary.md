# F9 — legacy reruns

## A. 2b same-exponent control (alpha=1.5/1.5, 3 reps)
- rep0: ZM rmse 0.2597, winner is, unit-amplitude helps: True
- rep1: ZM rmse 0.2600, winner is, unit-amplitude helps: True
- rep2: ZM rmse 0.2474, winner is, unit-amplitude helps: True

## B. E2-lite: winner stability under weighting — stable on 14/25
- CHANGES: Complete Works of Shakespeare: ols=is rank_w=exp freq_w=exp
- CHANGES: War and Peace: ols=is rank_w=exp freq_w=exp
- CHANGES: King James Bible: ols=is rank_w=exp freq_w=exp
- CHANGES: Federalist Papers: ols=is rank_w=exp freq_w=exp
- CHANGES: The Iliad: ols=is rank_w=exp freq_w=exp
- CHANGES: Democracy in America: ols=is rank_w=exp freq_w=exp
- CHANGES: Origin of Species: ols=is rank_w=exp freq_w=exp
- CHANGES: Wealth of Nations: ols=is rank_w=exp freq_w=exp
- CHANGES: Les Miserables: ols=is rank_w=exp freq_w=exp
- CHANGES: Decline and Fall Vol 1: ols=is rank_w=exp freq_w=exp
- CHANGES: Principia Ethica: ols=is rank_w=exp freq_w=exp

## C. 10a transfers
- War and Peace → King James Bible: transfer 0.1498 vs ZM 0.1882 vs in-domain 0.1494 (beats ZM: True)
- Moby Dick → King James Bible: transfer 0.1799 vs ZM 0.1882 vs in-domain 0.1494 (beats ZM: True)
- Complete Works of Shakespeare → War and Peace: transfer 0.1513 vs ZM 0.1830 vs in-domain 0.1458 (beats ZM: True)