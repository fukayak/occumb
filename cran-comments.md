# occumb (version 1.1.0)

## Test environments

- local: Linux Mint 20.2 R 4.3.3
- win-builder (release, devel)
- mac-builder (release)

## R CMD check results

### local

0 errors | 1 warning | 0 notes

* checking data for ASCII and uncompressed saves ...
WARNING
‘qpdf’ is needed for checks on size reduction of PDFs

This warning appears locally but not in win-builder and mac-builder.

### win-builder (release)

0 errors | 0 warnings | 0 note

### win-builder (devel)

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... [12s] NOTE
Maintainer: 'Keiichi Fukaya <fukayak99@gmail.com>'

### mac-builder

0 errors | 1 warnings | 0 note

* checking whether package ‘occumb’ can be installed ... [1s/1s] WARNING
Found the following significant warnings:
  Warning: package ‘jagsUI’ was built under R version 4.3.1



# occumb (version 1.0.3)

## Test environments

- local: Linux Mint 20.1 R 4.3.2
- win-builder (release, devel)
- mac-builder (release)

## R CMD check results

### local

0 errors | 1 warning | 0 notes

* checking data for ASCII and uncompressed saves ...
WARNING
‘qpdf’ is needed for checks on size reduction of PDFs

This warning appears locally but not in win-builder and mac-builder.

### win-builder

0 errors | 0 warnings | 0 note

### mac-builder

0 errors | 0 warnings | 0 note



# occumb (version 1.0.2)

## Test environments

- local: Linux Mint 20.1 R 4.3.1
- win-builder (release, devel)
- mac-builder (release)

## R CMD check results

### local

0 errors | 1 warning | 1 notes

* checking data for ASCII and uncompressed saves ...
WARNING
‘qpdf’ is needed for checks on size reduction of PDFs

This warning appears locally but not in win-builder and mac-builder.

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable

This is likely a temporary issue because version 4.4.0 of V8 package is available in CRAN.

### win-builder

0 errors | 0 warnings | 0 note

### mac-builder

0 errors | 0 warnings | 0 note



# occumb (version 1.0.1)

## Resubmission

This is a resubmission. In this version I have:

- Added a reference describing the methods in the package to the Description field of the `DESCRIPTION` file.

- Added the Value field to `.Rd` files regarding exported methods (i.e., `plot` and `summary`).

- Removed `\dontrun` in cases where the execution time of example codes was < 5 sec and otherwise replaced `\dontrun` with `\donttest`. Also, fixed bugs in the example codes found during this process.

## Test environments

- local: Linux Mint 20.1 R 4.3.1
- win-builder (release, devel)
- mac-builder (release)

## R CMD check results

### local

0 errors | 1 warning | 1 notes

* checking data for ASCII and uncompressed saves ...
WARNING
‘qpdf’ is needed for checks on size reduction of PDFs

This warning appears locally but not in win-builder and mac-builder. Hopefully it does not matter when compiled on CRAN.

* checking package dependencies ... NOTE
Package suggested but not available for checking: ‘spelling’

This is likely a temporary issue because version 2.2.1 of spelling package is available in CRAN.

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

