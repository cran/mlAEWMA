#' mlAEWMA: Machine Learning Based Adaptive EWMA Control Charts for Dispersion
#'
#' Adaptive EWMA control charts for the variance of a process whose smoothing
#' weight is predicted by a fitted regression model. Nine models in four
#' families are provided: neural networks ([mlaewma_neural()]), support vector
#' machines ([mlaewma_svm()]), trees and boosting ([mlaewma_tree()]) and dynamic
#' quantile regression ([mlaewma_dqr()]). [mlaewma_compare()] runs all nine and
#' [calculate_arl()] estimates run lengths by simulation.
#'
#' @references
#' Noor-ul-Amin, M., Kazmi, M. W., Alkhalaf, S., Abdel-Khalek, S. and Nabi, M.
#' (2024). Machine learning based parameter-free adaptive EWMA control chart
#' to monitor process dispersion. *Scientific Reports*, 14, 31271.
#' \doi{10.1038/s41598-024-82699-8}
#'
#' Kazmi, M. W. and Noor-ul-Amin, M. (2024). Adaptive EWMA control chart by
#' using support vector regression. *Quality and Reliability Engineering
#' International*, 40(7), 3831-3843. \doi{10.1002/qre.3603}
#'
#' Castagliola, P. (2005). A new S2-EWMA control chart for monitoring the
#' process variance. *Quality and Reliability Engineering International*,
#' 21(8), 781-794. \doi{10.1002/qre.686}
#'
#' Capizzi, G. and Masarotto, G. (2003). An adaptive exponentially weighted
#' moving average control chart. *Technometrics*, 45(3), 199-207.
#' \doi{10.1198/004017003000000023}
#'
#' Atar, G., Purutçuoğlu, V. and Bursalı, A. (2026). Comparison of machine
#' learning models in a parameter-free adaptive EWMA control chart for
#' process dispersion monitoring. Accepted at the 3rd International Congress
#' on Data Science and Information Technologies (INFTEC 2026), Nicosia,
#' TRNC, 22-24 October 2026.
#'
#' @keywords internal
"_PACKAGE"
