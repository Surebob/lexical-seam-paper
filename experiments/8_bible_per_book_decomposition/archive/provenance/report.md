# Angle 6: King James Bible by Book

- Each canonical book was fit independently with the soft-k Seam-Mandelbrot PMF.
- books analyzed: `66`
- per-book step-2 help count: `6` / `66`
- per-book soft-k beats ZM on held-out avg NLL: `33` / `66`
- per-book soft-k beats MOE on held-out avg NLL: `45` / `66`
- median per-book soft-k minus ZM held-out avg NLL: `-0.000133398377`
- median per-book soft-k minus MOE held-out avg NLL: `-0.002803782823`
- median per-book step-2 gain: `-0.003613180110`

- whole-Bible single-fit soft-k held-out avg NLL: `6.016015736045`
- aggregate per-book soft-k held-out avg NLL: `5.604219509299`

| book | soft-k - ZM | soft-k - MOE | step-2 gain | step-2 expr |
| --- | ---: | ---: | ---: | --- |
| 2 Timothy | -0.014556 | -0.037388 | -0.008146 | eml[sub[x,1],eml[x,1]] |
| 1 Peter | -0.014447 | -0.027491 | -0.005852 | eml[sub[x,1],eml[x,1]] |
| 1 Thessalonians | -0.005834 | -0.024783 | -0.007426 | eml[sub[x,1],eml[x,1]] |
| Philemon | -0.004238 | -0.023806 | -0.014455 | eml[sub[x,1],eml[x,1]] |
| Lamentations | -0.013188 | -0.022151 | -0.005328 | eml[sub[x,1],eml[x,1]] |
| 2 Corinthians | -0.000517 | -0.020931 | -0.001712 | sub[pow[x,x],sqrt[x]] |
| 2 Thessalonians | -0.009691 | -0.019928 | -0.009442 | eml[sub[x,1],eml[x,1]] |
| Hosea | -0.010711 | -0.019178 | -0.001243 | sub[pow[x,x],sqrt[x]] |
| James | -0.009138 | -0.019090 | -0.006964 | eml[sub[x,1],eml[x,1]] |
| 3 John | -0.006959 | -0.018594 | -0.008686 | sub[pow[x,x],sqrt[x]] |
| Galatians | -0.016741 | -0.018565 | 0.000724 | sub[pow[x,x],sqrt[x]] |
| 1 Timothy | 0.001128 | -0.015417 | -0.001207 | sub[pow[x,x],sqrt[x]] |
| Micah | -0.009599 | -0.014708 | -0.007010 | eml[sub[x,1],eml[x,1]] |
| 1 Corinthians | -0.007388 | -0.012536 | -0.001856 | sub[pow[x,x],sqrt[x]] |
| Genesis | -0.023117 | -0.010072 | -0.002529 | eml[sub[x,1],eml[x,1]] |
| Numbers | -0.018029 | -0.008899 | -0.003867 | eml[sub[x,1],eml[x,1]] |
| Malachi | 0.002093 | -0.008727 | -0.001879 | sub[pow[x,x],sqrt[x]] |
| Obadiah | 0.000250 | -0.008353 | -0.012499 | eml[sub[x,1],eml[x,1]] |
| Romans | -0.003209 | -0.007764 | 0.000793 | sub[pow[x,x],sqrt[x]] |
| Philippians | 0.017320 | -0.007600 | -0.007007 | eml[sub[x,1],eml[x,1]] |
| Joshua | -0.009083 | -0.006965 | -0.004585 | eml[sub[x,1],eml[x,1]] |
| Revelation | -0.011143 | -0.006927 | -0.004696 | eml[sub[x,1],eml[x,1]] |
| 2 Samuel | -0.010729 | -0.006003 | -0.003351 | eml[sub[x,1],eml[x,1]] |
| Job | 0.001544 | -0.005791 | -0.003125 | eml[sub[x,1],eml[x,1]] |
| Zephaniah | 0.004092 | -0.005501 | -0.007625 | eml[sub[x,1],eml[x,1]] |
| Amos | 0.002116 | -0.004412 | -0.005662 | eml[sub[x,1],eml[x,1]] |
| Hebrews | 0.004615 | -0.003888 | 0.000016 | sub[pow[x,x],sqrt[x]] |
| 1 Samuel | -0.010964 | -0.003814 | -0.003016 | eml[sub[x,1],eml[x,1]] |
| Jeremiah | -0.011929 | -0.003419 | -0.001938 | eml[sub[x,1],eml[x,1]] |
| Luke | -0.005033 | -0.003074 | -0.001553 | eml[sub[x,1],eml[x,1]] |
| Judges | -0.005423 | -0.002880 | -0.003601 | eml[sub[x,1],eml[x,1]] |
| Ezekiel | -0.011456 | -0.002838 | -0.002173 | eml[sub[x,1],eml[x,1]] |
| Acts | -0.006979 | -0.002832 | -0.002715 | eml[sub[x,1],eml[x,1]] |
| Ecclesiastes | 0.006875 | -0.002776 | 0.000761 | sub[pow[x,x],sqrt[x]] |
| Leviticus | -0.010725 | -0.002559 | -0.003715 | eml[sub[x,1],eml[x,1]] |
| 1 Kings | -0.007107 | -0.002465 | -0.002550 | eml[sub[x,1],eml[x,1]] |
| 2 Kings | -0.004617 | -0.002166 | -0.004215 | eml[sub[x,1],eml[x,1]] |
| Matthew | -0.003776 | -0.002017 | -0.002202 | eml[sub[x,1],eml[x,1]] |
| Deuteronomy | -0.004110 | -0.001926 | -0.001805 | eml[sub[x,1],eml[x,1]] |
| 2 Chronicles | -0.002718 | -0.001898 | -0.003630 | eml[sub[x,1],eml[x,1]] |
| Isaiah | -0.001593 | -0.001730 | -0.001451 | eml[sub[x,1],eml[x,1]] |
| Ephesians | 0.005206 | -0.001377 | 0.000260 | sub[pow[x,x],sqrt[x]] |
| 1 Chronicles | 0.000542 | -0.001287 | -0.003843 | eml[sub[x,1],eml[x,1]] |
| Exodus | -0.006775 | -0.001236 | -0.002601 | eml[sub[x,1],eml[x,1]] |
| Nehemiah | 0.000601 | -0.001091 | -0.005348 | eml[sub[x,1],eml[x,1]] |
| Psalms | 0.001458 | 0.000148 | -0.002707 | eml[sub[x,1],eml[x,1]] |
| 1 John | 0.007165 | 0.000201 | -0.005316 | sub[pow[x,x],sqrt[x]] |
| Colossians | 0.018091 | 0.000835 | -0.002388 | sub[pow[x,x],sqrt[x]] |
| Joel | 0.006725 | 0.000876 | -0.007310 | eml[sub[x,1],eml[x,1]] |
| Ruth | 0.003784 | 0.001126 | -0.007192 | eml[sub[x,1],eml[x,1]] |
| Ezra | 0.004089 | 0.002002 | -0.006162 | eml[sub[x,1],eml[x,1]] |
| Mark | 0.003545 | 0.002981 | -0.003625 | eml[sub[x,1],eml[x,1]] |
| Esther | 0.004977 | 0.004828 | -0.005436 | eml[sub[x,1],eml[x,1]] |
| John | 0.004766 | 0.005348 | -0.004739 | eml[sub[x,1],eml[x,1]] |
| Titus | 0.040424 | 0.005505 | -0.003358 | sub[pow[x,x],sqrt[x]] |
| Proverbs | 0.014838 | 0.007218 | -0.001446 | sub[pow[x,x],sqrt[x]] |
| Haggai | 0.012673 | 0.008414 | -0.010643 | eml[sub[x,1],eml[x,1]] |
| Zechariah | 0.018984 | 0.010816 | -0.001899 | sub[pow[x,x],sqrt[x]] |
| Daniel | 0.016177 | 0.013395 | -0.003584 | sub[pow[x,x],sqrt[x]] |
| 2 Peter | 0.030298 | 0.016957 | -0.001148 | sub[pow[x,x],sqrt[x]] |
| Nahum | 0.020913 | 0.020086 | -0.009175 | eml[sub[x,1],eml[x,1]] |
| Song of Solomon | 0.037063 | 0.020399 | -0.000923 | sub[pow[x,x],sqrt[x]] |
| Habakkuk | 0.030283 | 0.023675 | 0.005712 | sub[pow[x,x],sqrt[x]] |
| Jonah | 0.044156 | 0.031162 | -0.009883 | eml[sub[x,1],eml[x,1]] |
| 2 John | 0.047763 | 0.035796 | -0.006864 | sub[pow[x,x],sqrt[x]] |
| Jude | 0.056021 | 0.052761 | -0.012119 | eml[sub[x,1],eml[x,1]] |
