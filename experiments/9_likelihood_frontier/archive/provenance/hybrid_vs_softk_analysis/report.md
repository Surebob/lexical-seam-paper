# Experiment C: Hybrid vs Soft-k Head-to-Head

- hybrid wins on `16` / 25 corpora
- soft-k wins on `9` / 25 corpora

## Correlations With delta = hybrid minus soft-k

- anthology-like: Pearson `-0.036296`, Spearman `-0.337691`
- IS-head phase: Pearson `-0.332554`, Spearman `-0.301718`
- log(1 + ZM c): Pearson `-0.157888`, Spearman `-0.154615`
- raw ZM c: Pearson `-0.154209`, Spearman `-0.154615`
- log vocabulary: Pearson `0.109866`, Spearman `0.064615`
- log token count: Pearson `-0.288220`, Spearman `-0.161538`
- soft-k transition fraction: Pearson `0.391177`, Spearman `0.143846`
- hybrid transition fraction: Pearson `0.206100`, Spearman `0.110769`
- log10 soft-k lambda: Pearson `-0.258221`, Spearman `-0.065215`

## By Winner

- hybrid-win corpora: median ZM c `2.905112`, anthology fraction `0.312`, IS fraction `0.562`
- soft-k-win corpora: median ZM c `2.849243`, anthology fraction `0.111`, IS fraction `0.222`

## Strongest Hybrid Wins

- Grimm's Fairy Tales: delta `-0.003727`, ZM c `4.423357`, bregman `exp`, anthology `1`
- King James Bible: delta `-0.002695`, ZM c `3.423773`, bregman `IS`, anthology `1`
- Jane Eyre: delta `-0.002415`, ZM c `3.749111`, bregman `exp`, anthology `0`
- Don Quixote: delta `-0.002026`, ZM c `4.877938`, bregman `exp`, anthology `0`
- Federalist Papers: delta `-0.001940`, ZM c `0.872445`, bregman `IS`, anthology `1`

## Strongest Soft-k Wins

- Ulysses: delta `0.019944`, ZM c `1.743818`, bregman `exp`, anthology `0`
- Aesop's Fables: delta `0.014239`, ZM c `1.177734`, bregman `exp`, anthology `1`
- Dubliners: delta `0.004266`, ZM c `2.849243`, bregman `exp`, anthology `0`
- Pride and Prejudice: delta `0.004199`, ZM c `8.000332`, bregman `exp`, anthology `0`
- Moby Dick: delta `0.002309`, ZM c `2.143819`, bregman `exp`, anthology `0`
