# Experiment A: Mechanism-Preserving Hybrid

- Hybrid head-ZM + tail-MOE PMF refit with an added in-sample step-2 gain penalty.
- Objective adds `lambda_mech * max(0, train_step2_gain)^2` on top of the existing hybrid k-regularization.
- Held-out step-2 help is evaluated on held-out residuals using the train-fitted PMF and held-out positive-count ranks.

## lambda_mech = 0.01

- held-out wins vs MOE: `24` / 25
- held-out wins vs soft-k: `22` / 25
- held-out wins vs ZM: `16` / 25
- four-way winner counts: `{"zipf": 0, "zm": 9, "moe": 0, "mech": 16}`
- median mech minus MOE held-out avg NLL: `-0.003355932154`
- median mech minus ZM held-out avg NLL: `-0.002270694371`
- median mech minus soft-k held-out avg NLL: `-0.001753708762`
- held-out step-2 help count: `25` / 25

## lambda_mech = 0.1

- held-out wins vs MOE: `20` / 25
- held-out wins vs soft-k: `20` / 25
- held-out wins vs ZM: `13` / 25
- four-way winner counts: `{"zipf": 0, "zm": 9, "moe": 3, "mech": 13}`
- median mech minus MOE held-out avg NLL: `-0.003966316227`
- median mech minus ZM held-out avg NLL: `-0.000047463895`
- median mech minus soft-k held-out avg NLL: `-0.001723305944`
- held-out step-2 help count: `25` / 25

## lambda_mech = 1

- held-out wins vs MOE: `16` / 25
- held-out wins vs soft-k: `15` / 25
- held-out wins vs ZM: `14` / 25
- four-way winner counts: `{"zipf": 0, "zm": 10, "moe": 3, "mech": 12}`
- median mech minus MOE held-out avg NLL: `-0.001905136380`
- median mech minus ZM held-out avg NLL: `-0.001882157623`
- median mech minus soft-k held-out avg NLL: `-0.001044850161`
- held-out step-2 help count: `25` / 25

## lambda_mech = 10

- held-out wins vs MOE: `15` / 25
- held-out wins vs soft-k: `15` / 25
- held-out wins vs ZM: `15` / 25
- four-way winner counts: `{"zipf": 0, "zm": 9, "moe": 4, "mech": 12}`
- median mech minus MOE held-out avg NLL: `-0.001858563700`
- median mech minus ZM held-out avg NLL: `-0.001193117046`
- median mech minus soft-k held-out avg NLL: `-0.001231101222`
- held-out step-2 help count: `25` / 25

## lambda_mech = 100

- held-out wins vs MOE: `9` / 25
- held-out wins vs soft-k: `11` / 25
- held-out wins vs ZM: `11` / 25
- four-way winner counts: `{"zipf": 0, "zm": 11, "moe": 7, "mech": 7}`
- median mech minus MOE held-out avg NLL: `0.002295201896`
- median mech minus ZM held-out avg NLL: `0.001667842911`
- median mech minus soft-k held-out avg NLL: `0.004726844033`
- held-out step-2 help count: `25` / 25

- No lambda_mech satisfied the decision rule `(held-out step-2 help <= 4/25) AND (beats MOE >= 18/25)`.
