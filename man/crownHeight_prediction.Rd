% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crown_modul.R
\name{crownHeight_prediction}
\alias{crownHeight_prediction}
\title{crownHeight_prediction}
\usage{
crownHeight_prediction(
  df_fit,
  df_predict,
  site_vars = site_vars,
  species_n_threshold = 100,
  k = 10,
  eval_model_crownHeight = TRUE,
  crownHeight_model = "lm",
  BRNN_neurons = 3,
  blocked_cv = TRUE
)
}
\arguments{
\item{df_fit}{data frame with tree heights and basal areas for individual
trees}

\item{df_predict}{data frame which will be used for predictions}

\item{site_vars}{optional, character vector with names of site variables}

\item{species_n_threshold}{a positive integer defining the minimum number of
observations required to treat a species as an independent group}

\item{k}{the number of folds to be used in the k fold cross-validation}

\item{eval_model_crownHeight}{logical, should the crown height model be
evaluated and returned as the output}

\item{crownHeight_model}{character string defining the model to be used for
crown heights. Available are ANN with Bayesian regularization (brnn) or
linear regression (lm)}

\item{BRNN_neurons}{positive integer defining the number of neurons to be
used in the brnn method.}

\item{blocked_cv}{logical, should the blocked cross-validation be used in the
evaluation phase?}
}
\value{
a list with four elements:
\enumerate{
 \item $predicted_crownHeight - a data frame with imputed crown heights
 \item $eval_crownHeight - a data frame with predicted and observed crown heights, or a character string indicating that crown height model was not evaluated
 \item $model_species - the output model for crown heights (species level)
 \item $model_speciesGroups - the output model for crown heights (species group level)
}
}
\description{
Model for predicting crown height
}
\examples{
library(MLFS)
data(data_tree_heights)
data(data_v3)

# A) Example with linear model
Crown_h_predictions <- crownHeight_prediction(df_fit = data_tree_heights,
    df_predict = data_v3,
    crownHeight_model = "lm",
    site_vars = c(),
    species_n_threshold = 100,
    k = 10, blocked_cv = TRUE,
    eval_model_crownHeight = TRUE)

predicted_df <- Crown_h_predictions$predicted_crownHeight # df with imputed heights
evaluation_df <- Crown_h_predictions$eval_crownHeight # df with evaluation results

# B) Example with non-linear BRNN model
Crown_h_predictions <- crownHeight_prediction(df_fit = data_tree_heights,
    df_predict = data_v3,
    crownHeight_model = "brnn",
    BRNN_neurons = 3,
    site_vars = c(),
    species_n_threshold = 100,
    k = 10, blocked_cv = TRUE,
    eval_model_crownHeight = TRUE)

}
