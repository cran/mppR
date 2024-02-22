#Startup function
#this function is executed once the library is loaded
#'
#' @export

.onAttach = function(...)
  {
  packageStartupMessage("mppR has been developed with the support of")
  packageStartupMessage("Wageningen University")
  packageStartupMessage("KWS SAAT SE & Co. KGaA")
  packageStartupMessage("Swiss National Science Foundation")
  packageStartupMessage("(Grant Postdoc.Mobility P500PB_203030)")
}