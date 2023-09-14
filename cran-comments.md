# occumb (version 1.0.0)

## Test environments

- local: Linux Mint 20.2 R 4.3.1
- win-builder (release)
- mac-builder (release)

## R CMD check results

### local

0 errors | 1 warning | 0 notes

* checking data for ASCII and uncompressed saves ...
WARNING
‘qpdf’ is needed for checks on size reduction of PDFs

This warning appears locally but not in win-builder and mac-builder. Hopefully it does not matter when compiled on CRAN.

### win-builder

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... [7s] NOTE
Maintainer: 'Keiichi Fukaya <fukayak99@gmail.com>'
New submission
Possibly misspelled words in DESCRIPTION:
  Metabarcoding (2:54)
  metabarcoding (10:5)
  multispecies (9:19)

These words are not misspelled.

### mac-builder

0 errors | 0 warnings | 0 note

