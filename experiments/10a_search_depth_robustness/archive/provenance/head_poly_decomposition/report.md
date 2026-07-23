# Zipf Head Polynomial Decomposition

- Model tested: `ZM + ((x - 1) - log(x)) + poly_d(x)` fit on the full corpus.
- `poly_d(x)` is least-squares fit on the remaining full-corpus residual after ZM plus the fixed step-2 term.

## Shakespeare

- ZM alone RMSE: `0.184430425006`
- Step-2 only RMSE: `0.182851543528`
- Step-10 monster RMSE: `0.159871852997`
- Step-2 + poly degree 3 RMSE: `0.168711743486`
  coefficients: `-15.763455504217, 34.167830373167, -23.267018182264, 4.796014266327`
- Step-2 + poly degree 4 RMSE: `0.168704763154`
  coefficients: `-1.436826925293, -11.908838214828, 30.508260642776, -21.846285693606, 4.615282012436`
- Step-2 + poly degree 5 RMSE: `0.150549619124`
  coefficients: `286.438070995346, -912.078297072259, 1080.512379675886, -574.144816867051, 126.396316734489, -7.097188429179`

## War and Peace

- ZM alone RMSE: `0.182982301215`
- Step-2 only RMSE: `0.180436822151`
- Step-10 monster RMSE: `0.145347404652`
- Step-2 + poly degree 3 RMSE: `0.156221916805`
  coefficients: `-19.122079531307, 41.067854518246, -27.655789526727, 5.620136649880`
- Step-2 + poly degree 4 RMSE: `0.154030149021`
  coefficients: `-22.640939099529, 41.068041186545, -15.398628198977, -6.101606640490, 2.952689759501`
- Step-2 + poly degree 5 RMSE: `0.145812713187`
  coefficients: `173.536913384938, -570.260445775193, 691.935073800698, -371.297343023335, 79.651249636619, -3.623287366774`
