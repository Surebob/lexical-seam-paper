# Extended Widened-Grammar Low-c Manifold Check

- Selection note: English side uses the previously established low-c family (original step-2 exp-Bregman winner) because the literal saved single-ZM c<2.5 threshold only selects Ulysses and does not match the named low-c set.
- Protocol matches the original widened diagnostic through step 2: deterministic enumerative beam search, beam 50, keep-all-until-step-2, seeds `{x, 1}`, no coefficient fitting, extended grammar with trig/hyperbolic/erf/gamma/J0.
- Implementation note: this extension materializes only steps `1..2`, because later steps do not affect the step-2 beam under `keep-all-until-step-2` and the requested outputs are step-2-only.
- This extension uses raw ZM residuals only.

- manifold holds across `17/17` low-c corpora under the strict `R^2 > 0.9` rule.
- verdict counts: `yes=17`, `partial=0`, `no=0`

| corpus | ZM c | widened step-2 winner | cos_vs_exp | cos_vs_xx | 2D_span_R^2 | manifold_verdict |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| One Thousand and One Nights (Arabic, Wikisource) | 0.000000 | `sub[pow[x,x],sqrt[x]]` | 0.972777 | 1.000000 | 1.000000 | `yes` |
| Aesop's Fables | 3.687693 | `sub[erf[x],sin[x]]` | -0.942769 | -0.986976 | 0.977716 | `yes` |
| Arabian Nights (Vol 1) | 6.425403 | `eml[sub[x,1],eml[x,1]]` | 1.000000 | 0.974221 | 1.000000 | `yes` |
| Canterbury Tales | 15.250535 | `sub[erf[x],sin[x]]` | -0.970966 | -0.994007 | 0.988364 | `yes` |
| Collected Poe | 6.957528 | `sub[erf[x],sin[x]]` | -0.965352 | -0.993203 | 0.986457 | `yes` |
| Complete Sherlock Holmes | 10.294480 | `sub[erf[x],sin[x]]` | -0.951031 | -0.990290 | 0.981545 | `yes` |
| Critique of Pure Reason | 57.958304 | `sub[sin[x],erf[x]]` | 0.953213 | 0.990332 | 0.981864 | `yes` |
| Don Quixote | 65.663942 | `sub[erf[x],sin[x]]` | -0.970331 | -0.993833 | 0.987932 | `yes` |
| Dubliners | 3.263282 | `eml[sub[x,1],eml[x,1]]` | 1.000000 | 0.967881 | 1.000000 | `yes` |
| Emile | 51.480679 | `sub[erf[x],sin[x]]` | -0.962554 | -0.992451 | 0.984986 | `yes` |
| Grimm's Fairy Tales | 45.982513 | `sub[sin[x],erf[x]]` | 0.937146 | 0.985007 | 0.975952 | `yes` |
| Jane Eyre | 16.307585 | `sub[erf[x],sin[x]]` | -0.963677 | -0.992558 | 0.985171 | `yes` |
| Moby Dick | 10.497932 | `sub[erf[x],sin[x]]` | -0.972718 | -0.994718 | 0.989855 | `yes` |
| Pride and Prejudice | 34.119421 | `sub[erf[x],sin[x]]` | -0.940556 | -0.986442 | 0.976390 | `yes` |
| Ulysses | 0.632415 | `eml[sub[x,1],eml[x,1]]` | 1.000000 | 0.973988 | 1.000000 | `yes` |
| Romance of the Three Kingdoms (Chinese, Gutenberg 23950) | 0.000000 | `sub[pow[x,x],sqrt[x]]` | 0.974192 | 1.000000 | 1.000000 | `yes` |
| War and Peace (Russian, Wikisource) | 0.000000 | `sub[pow[x,x],sqrt[x]]` | 0.974310 | 1.000000 | 1.000000 | `yes` |

## One Thousand and One Nights (Arabic, Wikisource)

- language: `Arabic`
- ZM c: `0.000000000000`
- original-grammar step-2 winner: `sub[pow[x,x],sqrt[x]]`
- widened step-2 winner: `sub[pow[x,x],sqrt[x]]`
- widened step-2 math: `(pow(x,x)-sqrt(x))`
- widened step-2 RMSE: `0.203593864943`
- widened matches original grammar winner: `True`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `0.972777142021`
- weighted centered cosine vs x^x-sqrt(x): `1.000000000000`
- weighted centered cosine vs IS-Bregman: `0.990180445326`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `1.000000000000`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.203593864943 | False |
| 2 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.204883634906 | False |
| 3 | `mul[sub[1,x],sub[1,x]]` | `((1-x)*(1-x))` | 0.204888719820 | False |
| 4 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.205081093156 | False |
| 5 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.205371943512 | False |
| 6 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.206985790728 | True |
| 7 | `mul[sub[x,1],log[x]]` | `((x-1)*log(x))` | 0.207525743558 | False |
| 8 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.208202352136 | True |
| 9 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.210585513014 | False |
| 10 | `mul[sub[x,1],mul[x,x]]` | `((x-1)*(x*x))` | 0.211548749968 | False |

## Aesop's Fables

- language: `English`
- ZM c: `3.687692675564`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[erf[x],sin[x]]`
- widened step-2 math: `(erf(x)-sin(x))`
- widened step-2 RMSE: `0.168174716300`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `-0.942769218794`
- weighted centered cosine vs x^x-sqrt(x): `-0.986976267529`
- weighted centered cosine vs IS-Bregman: `-0.978687611223`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.977716140945`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.168174716300 | True |
| 2 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.168324944818 | True |
| 3 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.169229540562 | False |
| 4 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.171598714015 | False |
| 5 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.172399601262 | False |
| 6 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.174794575752 | False |
| 7 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.174939224236 | False |
| 8 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.174971606065 | False |
| 9 | `sub[div[x,1],pow[x,x]]` | `((x/1)-pow(x,x))` | 0.175495674851 | False |
| 10 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.177143891100 | False |

## Arabian Nights (Vol 1)

- language: `English`
- ZM c: `6.425403186436`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 math: `EML((x-1),EML(x,1))`
- widened step-2 RMSE: `0.182786889897`
- widened matches original grammar winner: `True`
- widened winner is Bregman (`IS` or `exp`): `True`
- weighted centered cosine vs exp-Bregman: `1.000000000000`
- weighted centered cosine vs x^x-sqrt(x): `0.974220931487`
- weighted centered cosine vs IS-Bregman: `0.937911539691`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `1.000000000000`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.182786889897 | False |
| 2 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.182826378906 | True |
| 3 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.182849716229 | True |
| 4 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.184569963574 | False |
| 5 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.185598857047 | False |
| 6 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.186249532567 | False |
| 7 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.186944177029 | False |
| 8 | `mul[sub[1,x],sub[1,x]]` | `((1-x)*(1-x))` | 0.187324523018 | False |
| 9 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.187676049126 | False |
| 10 | `sub[div[x,1],pow[x,x]]` | `((x/1)-pow(x,x))` | 0.187692874179 | False |

## Canterbury Tales

- language: `English`
- ZM c: `15.250535352833`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[erf[x],sin[x]]`
- widened step-2 math: `(erf(x)-sin(x))`
- widened step-2 RMSE: `0.167147186914`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `-0.970965595501`
- weighted centered cosine vs x^x-sqrt(x): `-0.994007444177`
- weighted centered cosine vs IS-Bregman: `-0.976801244531`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.988364429434`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.167147186914 | True |
| 2 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.167382215552 | True |
| 3 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.167824228316 | False |
| 4 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.170167788125 | False |
| 5 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.170570692565 | False |
| 6 | `sub[div[x,1],pow[x,x]]` | `((x/1)-pow(x,x))` | 0.170896120853 | False |
| 7 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.170936239689 | False |
| 8 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.171137046263 | False |
| 9 | `add[sub[1,x],log[x]]` | `((1-x)+log(x))` | 0.171557667877 | False |
| 10 | `mul[sub[1,x],sub[x,1]]` | `((1-x)*(x-1))` | 0.172329928701 | False |

## Collected Poe

- language: `English`
- ZM c: `6.957528471369`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[erf[x],sin[x]]`
- widened step-2 math: `(erf(x)-sin(x))`
- widened step-2 RMSE: `0.173153722601`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `-0.965351987117`
- weighted centered cosine vs x^x-sqrt(x): `-0.993202896440`
- weighted centered cosine vs IS-Bregman: `-0.980148002628`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.986456807437`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.173153722601 | True |
| 2 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.173639311534 | True |
| 3 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.174648939637 | False |
| 4 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.175450873254 | False |
| 5 | `sub[div[x,1],pow[x,x]]` | `((x/1)-pow(x,x))` | 0.177225929308 | False |
| 6 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.177809687212 | False |
| 7 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.178794149594 | False |
| 8 | `mul[sub[1,x],sub[x,1]]` | `((1-x)*(x-1))` | 0.178843440321 | False |
| 9 | `add[sub[1,x],log[x]]` | `((1-x)+log(x))` | 0.179031418767 | False |
| 10 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.179114214436 | False |

## Complete Sherlock Holmes

- language: `English`
- ZM c: `10.294479612917`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[erf[x],sin[x]]`
- widened step-2 math: `(erf(x)-sin(x))`
- widened step-2 RMSE: `0.168026734170`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `-0.951030875490`
- weighted centered cosine vs x^x-sqrt(x): `-0.990289811541`
- weighted centered cosine vs IS-Bregman: `-0.978065455092`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.981545249621`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.168026734170 | True |
| 2 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.168317252441 | True |
| 3 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.169219698739 | False |
| 4 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.171087576614 | False |
| 5 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.172091578429 | False |
| 6 | `sub[div[x,1],pow[x,x]]` | `((x/1)-pow(x,x))` | 0.173102097105 | False |
| 7 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.173873441627 | False |
| 8 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.174037304517 | False |
| 9 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.174300779585 | False |
| 10 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.174513848556 | False |

## Critique of Pure Reason

- language: `English`
- ZM c: `57.958303594664`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[sin[x],erf[x]]`
- widened step-2 math: `(sin(x)-erf(x))`
- widened step-2 RMSE: `0.176414607275`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `0.953213005226`
- weighted centered cosine vs x^x-sqrt(x): `0.990331740637`
- weighted centered cosine vs IS-Bregman: `0.979149139702`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.981864346518`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.176414607275 | True |
| 2 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.176473125495 | False |
| 3 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.176654474198 | True |
| 4 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.177089563255 | False |
| 5 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.177471112049 | False |
| 6 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.179747369537 | False |
| 7 | `mul[sub[1,x],sub[1,x]]` | `((1-x)*(1-x))` | 0.182014390896 | False |
| 8 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.182260580714 | False |
| 9 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.182361019849 | False |
| 10 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.182551617616 | False |

## Don Quixote

- language: `English`
- ZM c: `65.663941885764`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[erf[x],sin[x]]`
- widened step-2 math: `(erf(x)-sin(x))`
- widened step-2 RMSE: `0.175393453160`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `-0.970330717246`
- weighted centered cosine vs x^x-sqrt(x): `-0.993832824102`
- weighted centered cosine vs IS-Bregman: `-0.977477137830`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.987931518070`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.175393453160 | True |
| 2 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.175554057352 | True |
| 3 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.176335684153 | False |
| 4 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.176927874023 | False |
| 5 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.177782754952 | False |
| 6 | `sub[div[x,1],pow[x,x]]` | `((x/1)-pow(x,x))` | 0.178871591007 | False |
| 7 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.179168857322 | False |
| 8 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.179627786536 | False |
| 9 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.180139476641 | False |
| 10 | `mul[sub[1,x],sub[x,1]]` | `((1-x)*(x-1))` | 0.180398644290 | False |

## Dubliners

- language: `English`
- ZM c: `3.263282326596`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 math: `EML((x-1),EML(x,1))`
- widened step-2 RMSE: `0.175189799539`
- widened matches original grammar winner: `True`
- widened winner is Bregman (`IS` or `exp`): `True`
- weighted centered cosine vs exp-Bregman: `1.000000000000`
- weighted centered cosine vs x^x-sqrt(x): `0.967880939522`
- weighted centered cosine vs IS-Bregman: `0.927391860475`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `1.000000000000`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.175189799539 | False |
| 2 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.175441469879 | True |
| 3 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.175722721509 | True |
| 4 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.176329148371 | False |
| 5 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.178698650196 | False |
| 6 | `mul[sub[1,x],sub[1,x]]` | `((1-x)*(1-x))` | 0.180770489689 | False |
| 7 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.180842299637 | False |
| 8 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.181085624717 | False |
| 9 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.181437086892 | False |
| 10 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.181777848554 | False |

## Emile

- language: `English`
- ZM c: `51.480679153812`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[erf[x],sin[x]]`
- widened step-2 math: `(erf(x)-sin(x))`
- widened step-2 RMSE: `0.179430495661`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `-0.962554070882`
- weighted centered cosine vs x^x-sqrt(x): `-0.992450864114`
- weighted centered cosine vs IS-Bregman: `-0.977908016341`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.984986486968`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.179430495661 | True |
| 2 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.179566983788 | True |
| 3 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.180338373379 | False |
| 4 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.181539251240 | False |
| 5 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.182204486164 | False |
| 6 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.182904259564 | False |
| 7 | `sub[div[x,1],pow[x,x]]` | `((x/1)-pow(x,x))` | 0.183649620892 | False |
| 8 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.183873689823 | False |
| 9 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.184425551075 | False |
| 10 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.185101418388 | False |

## Grimm's Fairy Tales

- language: `English`
- ZM c: `45.982512741972`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[sin[x],erf[x]]`
- widened step-2 math: `(sin(x)-erf(x))`
- widened step-2 RMSE: `0.164796636043`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `0.937146489224`
- weighted centered cosine vs x^x-sqrt(x): `0.985007132136`
- weighted centered cosine vs IS-Bregman: `0.977695225465`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.975951978697`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.164796636043 | True |
| 2 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.164991436277 | True |
| 3 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.165504850724 | False |
| 4 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.167056462899 | False |
| 5 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.167335168714 | False |
| 6 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.170095413830 | False |
| 7 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.170260637427 | False |
| 8 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.171385602051 | False |
| 9 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.171868001662 | False |
| 10 | `mul[sub[1,x],sub[1,x]]` | `((1-x)*(1-x))` | 0.172998752660 | False |

## Jane Eyre

- language: `English`
- ZM c: `16.307585257849`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[erf[x],sin[x]]`
- widened step-2 math: `(erf(x)-sin(x))`
- widened step-2 RMSE: `0.161839671517`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `-0.963677396459`
- weighted centered cosine vs x^x-sqrt(x): `-0.992557513077`
- weighted centered cosine vs IS-Bregman: `-0.976557532634`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.985170888403`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.161839671517 | True |
| 2 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.162250103404 | True |
| 3 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.163330107971 | False |
| 4 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.163956239344 | False |
| 5 | `sub[div[x,1],pow[x,x]]` | `((x/1)-pow(x,x))` | 0.165375523703 | False |
| 6 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.166783505173 | False |
| 7 | `mul[sub[1,x],sub[x,1]]` | `((1-x)*(x-1))` | 0.166934720135 | False |
| 8 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.167072557922 | False |
| 9 | `add[sub[1,x],log[x]]` | `((1-x)+log(x))` | 0.167205257156 | False |
| 10 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.167450260797 | False |

## Moby Dick

- language: `English`
- ZM c: `10.497932011426`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[erf[x],sin[x]]`
- widened step-2 math: `(erf(x)-sin(x))`
- widened step-2 RMSE: `0.162126475177`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `-0.972717871054`
- weighted centered cosine vs x^x-sqrt(x): `-0.994718304972`
- weighted centered cosine vs IS-Bregman: `-0.979686982521`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.989854823468`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.162126475177 | True |
| 2 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.162497224226 | True |
| 3 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.163345757725 | False |
| 4 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.164335235061 | False |
| 5 | `sub[div[x,1],pow[x,x]]` | `((x/1)-pow(x,x))` | 0.165273131933 | False |
| 6 | `add[sub[1,x],log[x]]` | `((1-x)+log(x))` | 0.166190141470 | False |
| 7 | `mul[sub[1,x],sub[x,1]]` | `((1-x)*(x-1))` | 0.166644482153 | False |
| 8 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.166703594345 | False |
| 9 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.166844607438 | False |
| 10 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.167205907572 | False |

## Pride and Prejudice

- language: `English`
- ZM c: `34.119420842756`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `sub[erf[x],sin[x]]`
- widened step-2 math: `(erf(x)-sin(x))`
- widened step-2 RMSE: `0.166759034021`
- widened matches original grammar winner: `False`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `-0.940556346681`
- weighted centered cosine vs x^x-sqrt(x): `-0.986441621684`
- weighted centered cosine vs IS-Bregman: `-0.975512601277`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `0.976389870575`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.166759034021 | True |
| 2 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.166872198494 | True |
| 3 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.167839264182 | False |
| 4 | `sub[sqrt[x],pow[x,x]]` | `(sqrt(x)-pow(x,x))` | 0.170269310529 | False |
| 5 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.170407186984 | False |
| 6 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.171361271351 | False |
| 7 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.172476794415 | False |
| 8 | `sub[div[x,1],pow[x,x]]` | `((x/1)-pow(x,x))` | 0.172612187992 | False |
| 9 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.172786892910 | False |
| 10 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.173283694728 | False |

## Ulysses

- language: `English`
- ZM c: `0.632415212219`
- original-grammar step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 winner: `eml[sub[x,1],eml[x,1]]`
- widened step-2 math: `EML((x-1),EML(x,1))`
- widened step-2 RMSE: `0.183018885952`
- widened matches original grammar winner: `True`
- widened winner is Bregman (`IS` or `exp`): `True`
- weighted centered cosine vs exp-Bregman: `1.000000000000`
- weighted centered cosine vs x^x-sqrt(x): `0.973987778958`
- weighted centered cosine vs IS-Bregman: `0.936315689422`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `1.000000000000`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.183018885952 | False |
| 2 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.183579227995 | False |
| 3 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.183986063447 | True |
| 4 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.184192053333 | False |
| 5 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.184496783476 | True |
| 6 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.184706126540 | False |
| 7 | `mul[sub[1,x],sub[1,x]]` | `((1-x)*(1-x))` | 0.184937892886 | False |
| 8 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.188096808392 | False |
| 9 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.189045527066 | False |
| 10 | `add[sub[1,x],log[x]]` | `((1-x)+log(x))` | 0.189076660663 | False |

## Romance of the Three Kingdoms (Chinese, Gutenberg 23950)

- language: `Mandarin`
- ZM c: `0.000000000000`
- original-grammar step-2 winner: `sub[pow[x,x],sqrt[x]]`
- widened step-2 winner: `sub[pow[x,x],sqrt[x]]`
- widened step-2 math: `(pow(x,x)-sqrt(x))`
- widened step-2 RMSE: `0.234078632609`
- widened matches original grammar winner: `True`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `0.974192152314`
- weighted centered cosine vs x^x-sqrt(x): `1.000000000000`
- weighted centered cosine vs IS-Bregman: `0.992023342036`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `1.000000000000`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.234078632609 | False |
| 2 | `mul[sub[1,x],sub[1,x]]` | `((1-x)*(1-x))` | 0.235271982657 | False |
| 3 | `mul[sub[x,1],log[x]]` | `((x-1)*log(x))` | 0.235640510080 | False |
| 4 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.235974708385 | False |
| 5 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.235984330519 | False |
| 6 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.237531066030 | False |
| 7 | `mul[sub[x,1],mul[x,x]]` | `((x-1)*(x*x))` | 0.239758163001 | False |
| 8 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.240230574904 | True |
| 9 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.242106148320 | True |
| 10 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.243239965651 | False |

## War and Peace (Russian, Wikisource)

- language: `Russian`
- ZM c: `0.000000000000`
- original-grammar step-2 winner: `sub[pow[x,x],sqrt[x]]`
- widened step-2 winner: `sub[pow[x,x],sqrt[x]]`
- widened step-2 math: `(pow(x,x)-sqrt(x))`
- widened step-2 RMSE: `0.196378600312`
- widened matches original grammar winner: `True`
- widened winner is Bregman (`IS` or `exp`): `False`
- weighted centered cosine vs exp-Bregman: `0.974310084379`
- weighted centered cosine vs x^x-sqrt(x): `1.000000000000`
- weighted centered cosine vs IS-Bregman: `0.991338120410`
- weighted centered span R^2 in span{exp, x^x-sqrt(x)}: `1.000000000000`
- manifold verdict: `yes`

| rank | expr | math | RMSE | new-op? |
| ---: | --- | --- | ---: | --- |
| 1 | `sub[pow[x,x],sqrt[x]]` | `(pow(x,x)-sqrt(x))` | 0.196378600312 | False |
| 2 | `eml[sub[x,1],eml[x,1]]` | `EML((x-1),EML(x,1))` | 0.197156484789 | False |
| 3 | `add[neg[x],pow[x,x]]` | `((-x)+pow(x,x))` | 0.197448065822 | False |
| 4 | `mul[sub[1,x],sub[1,x]]` | `((1-x)*(1-x))` | 0.197736852490 | False |
| 5 | `sub[sub[x,1],log[x]]` | `((x-1)-log(x))` | 0.198061463241 | False |
| 6 | `sub[sin[x],erf[x]]` | `(sin(x)-erf(x))` | 0.198762595176 | True |
| 7 | `sub[erf[x],sin[x]]` | `(erf(x)-sin(x))` | 0.199713478429 | True |
| 8 | `mul[sub[x,1],log[x]]` | `((x-1)*log(x))` | 0.201119368614 | False |
| 9 | `div[sub[x,1],eml[1,x]]` | `((x-1)/EML(1,x))` | 0.202544956386 | False |
| 10 | `div[sub[1,x],eml[1,x]]` | `((1-x)/EML(1,x))` | 0.204102881816 | False |
