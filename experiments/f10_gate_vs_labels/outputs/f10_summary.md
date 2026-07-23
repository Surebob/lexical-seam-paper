# F10 — gate vs labels (the identity capstone)

- corpora: 10 (all EXP03 label-complete)
- corr(log k_gate, log k_lab) = 0.1973; median k_gate/k_lab = 3.104
- corr(w_gate, w_lab) = -0.4275; median w_gate/w_lab = 1.769
- median gate-vs-pi curve RMSE = 0.2341 (pi is a share in [0,1])
- labeled k_POS ~ V^beta: beta = 0.298 95% CI [0.076, 0.519] (paper claim 0.545; n=10)
- labeled ablation: unit step-2 helps on full corpus 3/10 -> on content-only 4/10

## Inter-annotator pilot (Principia, 3 taggers vs manual)
- tagger1_principia_ethica vs manual: exact tag 91.8%, closed/open 97.0% (n=500)
- tagger2_principia_ethica vs manual: exact tag 91.8%, closed/open 96.6% (n=500)
- tagger5_principia_ethica vs manual: exact tag 89.9%, closed/open 93.3% (n=119)