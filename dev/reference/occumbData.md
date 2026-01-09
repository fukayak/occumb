# Constructor for occumbData data class.

`occumbData()` creates a data list compatible with the model fitting
function
[`occumb()`](https://fukayak.github.io/occumb/dev/reference/occumb.md).
The element (i.e., covariate) names for `spec_cov`, `site_cov`, and
`repl_cov` must all be unique. If `y` has a `dimnames` attribute, it is
retained in the resulting `occumbData` object, and can be referenced in
subsequent analyses.

## Usage

``` r
occumbData(y, spec_cov = NULL, site_cov = NULL, repl_cov = NULL)
```

## Arguments

- y:

  A 3-D array or a dataframe of sequence read counts (`integer` values).
  An array's dimensions are ordered by species, site, and replicate, and
  may have a `dimnames` attribute. A dataframe's columns are ordered by
  species, site, replicate, and sequence read counts. The data for
  missing replicates are represented by zero vectors. `NA`s are not
  allowed.

- spec_cov:

  A named list of species covariates. Each covariate can be a vector of
  continuous (`numeric` or `integer`) or discrete (`logical`, `factor`,
  or `character`) variables whose length is `dim(y)[1]` (i.e., the
  number of species). The order of the species of the covariate values
  must correspond to that of the species dimension of `y`. `NA`s are not
  allowed.

- site_cov:

  A named list of site covariates. Each covariate can be a vector of
  continuous (`numeric` or `integer`) or discrete (`logical`, `factor`,
  or `character`) variables whose length is `dim(y)[1]` (i.e., the
  number of sites). The order of the sites of the covariate values must
  correspond to that of the site dimension of `y`. `NA`s are not
  allowed.

- repl_cov:

  A named list of replicate covariates. Each covariate can be a matrix
  of continuous (`numeric` or `integer`) or discrete (`logical` or
  `character`) variables with dimensions equal to `dim(y)[2:3]` (i.e.,
  number of sites \\\times\\ number of replicates). The order of the
  sites and replicates of the covariate values must correspond to that
  of the site and replicate dimensions of `y`. `NA`s are not allowed.

## Value

An S4 object of the `occumbData` class.

## Examples

``` r
# Generate the smallest random dataset (2 species * 2 sites * 2 reps)
I <- 2 # Number of species
J <- 2 # Number of sites
K <- 2 # Number of replicates
data <- occumbData(
    y = array(sample.int(I * J * K), dim = c(I, J, K)),
    spec_cov = list(cov1 = rnorm(I)),
    site_cov = list(cov2 = rnorm(J), cov3 = factor(1:J)),
    repl_cov = list(cov4 = matrix(rnorm(J * K), J, K))
)

# A case for named y (with species and site names)
y_named <- array(sample.int(I * J * K), dim = c(I, J, K))
dimnames(y_named) <- list(c("common species", "uncommon species"),
                          c("good site", "bad site"), NULL)
data_named <- occumbData(
    y = y_named,
    spec_cov = list(cov1 = rnorm(I)),
    site_cov = list(cov2 = rnorm(J), cov3 = factor(1:J)),
    repl_cov = list(cov4 = matrix(rnorm(J * K), J, K))
)
# A real data example
data(fish_raw)
fish <- occumbData(
    y = fish_raw$y,
    spec_cov = list(mismatch = fish_raw$mismatch),
    site_cov = list(riverbank = fish_raw$riverbank)
)

# Get an overview of the datasets
summary(data)
#> Sequence read counts: 
#>  Number of species, I = 2 
#>  Number of sites, J = 2 
#>  Maximum number of replicates per site, K = 2 
#>  Number of missing observations = 0 
#>  Number of replicates per site: 2 (average), 0 (sd) 
#>  Sequencing depth: 9 (average), 2.8 (sd) 
#> 
#> Species covariates: 
#>  cov1 (continuous) 
#> Site covariates: 
#>  cov2 (continuous), cov3 (categorical) 
#> Replicate covariates: 
#>  cov4 (continuous) 
#> 
#> Labels for species: 
#>  (None) 
#> Labels for sites: 
#>  (None) 
#> Labels for replicates: 
#>  (None) 
summary(data_named)
#> Sequence read counts: 
#>  Number of species, I = 2 
#>  Number of sites, J = 2 
#>  Maximum number of replicates per site, K = 2 
#>  Number of missing observations = 0 
#>  Number of replicates per site: 2 (average), 0 (sd) 
#>  Sequencing depth: 9 (average), 3.7 (sd) 
#> 
#> Species covariates: 
#>  cov1 (continuous) 
#> Site covariates: 
#>  cov2 (continuous), cov3 (categorical) 
#> Replicate covariates: 
#>  cov4 (continuous) 
#> 
#> Labels for species: 
#>  common species, uncommon species 
#> Labels for sites: 
#>  good site, bad site 
#> Labels for replicates: 
#>  (None) 
summary(fish)
#> Sequence read counts: 
#>  Number of species, I = 50 
#>  Number of sites, J = 50 
#>  Maximum number of replicates per site, K = 3 
#>  Number of missing observations = 6 
#>  Number of replicates per site: 2.88 (average), 0.33 (sd) 
#>  Sequencing depth: 77910 (average), 98034.7 (sd) 
#> 
#> Species covariates: 
#>  mismatch (continuous) 
#> Site covariates: 
#>  riverbank (categorical) 
#> Replicate covariates: 
#>  (None) 
#> 
#> Labels for species: 
#>  Abbottina rivularis, Acanthogobius lactipes, Acheilognathus macropterus, Acheilognathus rhombeus, Anguilla japonica, Biwia zezera, Carassius cuvieri, Carassius spp., Channa argus, Ctenopharyngodon idella, Cyprinus carpio, Gambusia affinis, Gnathopogon spp., Gymnogobius castaneus, Gymnogobius petschiliensis, Gymnogobius urotaenia, Hemibarbus spp., Hypomesus nipponensis, Hypophthalmichthys spp., Hyporhamphus intermedius, Ictalurus punctatus, Ischikauia steenackeri, Lepomis macrochirus macrochirus, Leucopsarion petersii, Megalobrama amblycephala, Micropterus dolomieu dolomieu, Micropterus salmoides, Misgurnus spp., Monopterus albus, Mugil cephalus cephalus, Mylopharyngodon piceus, Nipponocypris sieboldii, Nipponocypris temminckii, Opsariichthys platypus, Opsariichthys uncirostris uncirostris, Oryzias latipes, Plecoglossus altivelis altivelis, Pseudogobio spp., Pseudorasbora parva, Rhinogobius spp., Rhodeus ocellatus ocellatus, Salangichthys microdon, Sarcocheilichthys variegatus microoculus, Silurus asotus, Squalidus chankaensis biwae, Tachysurus tokiensis, Tanakia lanceolata, Tribolodon brandtii maruta, Tribolodon hakonensis, Tridentiger spp. 
#> Labels for sites: 
#>  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50 
#> Labels for replicates: 
#>  L, C, R 
```
