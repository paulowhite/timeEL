#' Acute Leukemia data from Freireich et al (1963)
#'
#' In this study there were 21 pairs of subjects, and within each pair one subject
#' received 6-mercaptopurine (6-MP) and one got placebo. The data are right censored.
#' See also Gehan (1965) and Thomas & Grunkemeier (1975) who used the data as an
#' illustrative example (ignoring the pairing).
#'
#' @format A data frame with 42 rows and 3 variables:
#' \describe{
#'   \item{time}{time in remission (in weeks)}
#'   \item{status}{event status, 1 is relapse, 0 is censored}
#'   \item{group}{treatment group: 0 (placebo) or 1 (6-MP)}
#' }
#' @source Data listed in Section 5 in Thomas & Grunkemeier (1975) and Section 11 in Gehan (1965)
#' @references
#' Freireich et al (1963) Blood 21(6):699-716
#'
#' Gehan (1965) Biometrika 52:203-223
#'
#' Thomas & Grunkemeier (1975) JASA 70(352): 865-871
#' 
"Freireich"



#' Melanoma competing risks data 
#'
#' These competing risks data relate to survival of patients after operation for malignant melanoma collected at Odense University Hospital between 1962 and
#' 1977. The data are a subsample of the 'melanoma' data of the 'timereg' package (patients who had a tumor thickness of less than 5 cm).
#'
#' @format A data frame with 173 rows and 2 variables:
#' \describe{
#'   \item{time}{time to event (in years)}
#'   \item{status}{event status, 1 is death due to malignant melanoma, 2 is death due to another cause and 0 is censored}
#' }
#' @source 'timereg' package
#' @references
#' Andersen PK, Skovgaard LT (2010) Regression with linear predictors. Springer, Berlin
#'
#' Drzewiecki K, Andersen PK (1982) Survival with malignant melanoma: a regression analysis of prognostic factors. Cancer 49:2414–2419
#' 
"melanoma5"


#' Bone Marrow Transplant Registry
#'
#' The data contain observations of 408 patients treated with HLA-identical sibling
#' bone marrow transplantation for myelodysplasia. The dataset is essentially a subset of
#' the 'bmt' data of the 'timereg' package (minor changes were introduced to break the ties).
#' 
#' @format A data frame with 408 rows and 3 variables:
#' \describe{
#'   \item{time}{time to event since transplant (in months)}
#'   \item{status}{event status, 1 is dead from treatment related causes, 2 is relapse , 0 is censored.}
#'   \item{group}{platelet level: 1 if more than 100 x 10^9 per L, 0 if less}
#' }
#' @source 'timereg' package
#' @references
#' Li, J., Le-Rademacher, J., & Zhang, M. J. (2014). Weighted comparison of two cumulative incidence functions with R-CIFsmry package. Computer methods and programs in biomedicine, 116(3), 205-214.
#'     
"BMTplat"

#' Bone Marrow Transplant Registry
#'
#' The data contain observations of 408 patients treated with HLA-identical sibling
#' bone marrow transplantation for myelodysplasia. The dataset is essentially a subset of
#' the 'bmt' data of the 'timereg' package (minor changes were introduced to break the ties).
#' 
#' @format A data frame with 408 rows and 3 variables:
#' \describe{
#'   \item{time}{time to event since transplant (in months)}
#'   \item{status}{event status, 1 is dead from treatment related causes, 2 is relapse , 0 is censored.}
#'   \item{group}{presence of T-cell depletion: 1 if present, 0 otherwise}
#' }
#' @source 'timereg' package
#' @references
#' Li, J., Le-Rademacher, J., & Zhang, M. J. (2014). Weighted comparison of two cumulative incidence functions with R-CIFsmry package. Computer methods and programs in biomedicine, 116(3), 205-214.
#'     
"BMTtcell"

#' Simulated competing risks data
#'
#' The data were simulated as described in Blanche & Eriksson (2023), using  scenario A with sample size n=100.
#' 
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{time}{time to event}
#'   \item{status}{event status, 1 is main event, 2 is competing event, 0 is censored.}
#'   \item{group}{group (1 or 0)}
#' }
#' @source Simulated data
#' @references
#' Blanche & Eriksson (2023). Empirical likelihood comparison of absolute risks.
#'     
"SimA100"
