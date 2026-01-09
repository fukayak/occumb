# Fish eDNA metabarcoding dataset

A dataset of fish eDNA metabarcoding collected in the Kasumigaura
watershed, Japan.

## Usage

``` r
fish_raw
```

## Format

A list containing an array of sequence read count, `y`, a vector of the
total number of mismatched bases in the priming region of the forward
and reverse primers for each species, `mismatch`, and a factor
indicating whether the riverbank at each site lacked aquatic and
riparian vegetation, `riverbank`. Sequence reads were obtained from
three replicates (collected from the center of the river and near the
left and right riverbanks) from 50 sites across the watershed, of which
read counts from six samples were missing. The resulting sequence counts
of 50 freshwater fish taxa detected were recorded.

## Source

K. Fukaya, N. I. Kondo, S. S. Matsuzaki, T. Kadoya (2021) Data from:
Multispecies site occupancy modeling and study design for spatially
replicated environmental DNA metabarcoding. Dryad Digital Repository.
[doi:10.5061/dryad.3bk3j9kkm](https://doi.org/10.5061/dryad.3bk3j9kkm)
