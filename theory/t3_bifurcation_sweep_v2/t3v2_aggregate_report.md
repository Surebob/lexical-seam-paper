# T3v2 Aggregate Report

T3v2 tests bifurcation structure using a break-rank sweep in an explicitly two-regime rank-frequency construction.

## Preflight

- status: `pass`
- c delta K=50 to K=5000: `74.2309`
- endpoint winners differ: `False`
- K `50`: c `15.6365`, b `1.52225`, winner `eml[sub[x,1],eml[x,1]]`, closest `exp`
- K `5000`: c `89.8674`, b `0.998332`, winner `eml[sub[x,1],eml[x,1]]`, closest `exp`

## Primary Sweep

- transition type: `sharp_IS_exp_bifurcation_candidate`
- winner identity changes: `2`
- closest-generator changes: `2`
- fitted c range: `15.6365` to `412.463`
- c-band 66-79 hits: `0`

## Winner Sequence

- K `50`: c `15.6365`, winner `eml[sub[x,1],eml[x,1]]`, closest `exp`, cos(IS) `0.9107`, cos(exp) `1.0000`
- K `100`: c `29.0241`, winner `eml[sub[x,1],eml[x,1]]`, closest `exp`, cos(IS) `0.9107`, cos(exp) `1.0000`
- K `200`: c `54.7906`, winner `eml[sub[x,1],eml[x,1]]`, closest `exp`, cos(IS) `0.9107`, cos(exp) `1.0000`
- K `500`: c `132.465`, winner `eml[sub[x,1],eml[x,1]]`, closest `exp`, cos(IS) `0.9107`, cos(exp) `1.0000`
- K `1000`: c `255.749`, winner `eml[sub[x,1],eml[x,1]]`, closest `exp`, cos(IS) `0.9107`, cos(exp) `1.0000`
- K `2000`: c `412.463`, winner `sub[sub[x,1],log[x]]`, closest `IS`, cos(IS) `1.0000`, cos(exp) `0.9107`
- K `5000`: c `89.8674`, winner `eml[sub[x,1],eml[x,1]]`, closest `exp`, cos(IS) `0.9107`, cos(exp) `1.0000`

## Interpretation

Primary shows a discrete IS/exp-related transition. Dense sweep around the transition is warranted next.
