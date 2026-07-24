# f22 — likelihood-space PMF bake-off (held-out NLL per token)

80/20 binomial splits x3; support/ranks from train; all models normalized over [1, V_train]. ZM, frozen lambda-ZM, MOEZipf carry 2 free params each; free lambda-ZM carries 3.

- lzm_frozen beats zm: 6/75 (median NLL delta -0.05361 nats/token)
- lzm_frozen beats moe: 3/75 (median NLL delta -0.06312 nats/token)
- lzm_free beats zm: 69/75 (median NLL delta +0.00609 nats/token)
- lzm_free beats moe: 60/75 (median NLL delta +0.00170 nats/token)
- lzm_free beats lzm_frozen: 72/75 (median NLL delta +0.06615 nats/token)