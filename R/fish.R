#' Fish eDNA metabarcoding dataset
#'
#' A dataset of fish eDNA metabarcoding collected in the Kasumigaura watershed, Japan.
#'
#' @format An occumbData class object containing sequence read count \code{y},
#'  \code{mismatch} as a species covariate, and \code{riverbank} as a site
#'  covariate.
#'  \code{mismatch} represents the total number of mismatched bases in the
#'  priming region of the forward and reverse primers for each species.
#'  \code{riverbank} indicates whether the riverbank of each site lacks aquatic
#'  and riparian vegetation.
#'  The sequence reads were obtained from three replicates (collected from the
#'  center of the river and near the left and right riverbanks, respectively)
#'  from 50 sites across the watershed, of which read counts from 6 samples are
#'  missing.
#'  The resulting sequence counts of the 50 freshwater fish taxa detected are
#'  recorded.
#' @source K. Fukaya, N. I. Kondo, S. S. Matsuzaki, T. Kadoya (2021) Data from: Multispecies site occupancy modeling and study design for spatially replicated environmental DNA metabarcoding. Dryad Digital Repository. \doi{10.5061/dryad.3bk3j9kkm}
"fish"
