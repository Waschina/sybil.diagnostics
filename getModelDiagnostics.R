#' Perform a set of sybil modelorg diagnostics tests.
#' 
#' @param mod An object of class modelorg.
#' @param top.x Maximum number of rows displayed in the "Top"-lists.
#' @return A list for the test results of \code{mod}. List item \code{rxn.info} returns results obtained for individual 
#'   reactions. Columns include besides identifiers the predcited flux and the reduced costs. List item \code{met.info}
#'   contains only identfiers at the moment, but will be extented soon.
getModelDiagnostics <- function(mod, top.x = 10) {
  require(data.table)
  require(stringr)
  require(sybil)
  
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
  
  # A - Construct metabolite analysis table
  met.info <- data.table(id = mod@met_id, name = mod@met_name)
  
  # B - Construct reaction analysis table
  rxn.info <- data.table(id = mod@react_id, name = mod@react_name,
                         lb = mod@lowbnd, ub = mod@uppbnd)
  
  #~~~~~~~~~~~~~~~~~~~~~#
  # Perform diagnostics #
  #~~~~~~~~~~~~~~~~~~~~~#
  
  # B 1 - MTF flux
  sol <- sybil::optimizeProb(mod, algorithm = "mtf")
  rxn.info$flux = sol@fluxdist@fluxes[1:mod@react_num]
  
  # B 2 - reduced costs
  rxn.info$red.costs <- getReducedCosts(mod)
  
  
  #~~~~~~~~~~~~~~~~~~~~~#
  # Report main results #
  #~~~~~~~~~~~~~~~~~~~~~#
  cat("Nutrients, whose infow reached limits (lb):\n\n")
  print(rxn.info[flux == lb & flux != 0])
  cat("\n")
  
  cat("Top", top.x, "reactions with the lowest reduced costs:\n\n")
  print(rxn.info[order(red.costs)][1:top.x])
  cat("\n")
  
  cat("Top", top.x, "exchange reactions with the lowest reduced costs:\n\n")
  print(rxn.info[grepl("^EX", id)][order(red.costs)][1:top.x])
  cat("\n")
  
  cat("Top", top.x, "exchange reactions with flux < 0 (i.e. uptake) and the lowest reduced costs:\n\n")
  print(rxn.info[grepl("^EX", id) & flux < 0][order(red.costs)][1:top.x])
  cat("\n")
  
  return(list(met.info = met.info,
              rxn.info = rxn.info))
}

getReducedCosts <- function(mod) {
  require(sybil)
  
  nrxns <- as.character(react_num(mod)-1)
  
  opt <- optimizeProb(mod,
                      solver = "cplexAPI",
                      #prCmd = list(c("getColsLowBnds", "LP_PROB", "1:77")),
                      poCmd = list(c("getDjCPLEX","LP_PROB@oobj@env","LP_PROB@oobj@lp","0", nrxns)))
  rc <- postProc(opt)
  rc <- rc@pa[[1]]
  
  return(rc)
}
