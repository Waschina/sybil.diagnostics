get.exchange.reactions <- function(mod) {
  require(sybil)
  require(data.table)
  ind <- grep("^EX_", mod@react_id)
  dt <- data.table(Exchange = mod@react_id,
                   Name     = mod@react_name,
                   lb       = mod@lowbnd)
  dt <- dt[grepl("^EX", Exchange)]
  
  dt[, is.Resource := lb < 0]
  dt <- dt[order(is.Resource, decreasing = T)]
  
  return(dt)
}

optimize.fluxes <- function(mod, exclude.unused = F) {
  sol <- sybil::optimizeProb(mod, algorithm = "mtf")
  
  dt <- data.table(rxn  = mod@react_id, 
                   name = mod@react_name,
                   flux = sol@fluxdist@fluxes[1:mod@react_num])
  
  # Add EC information if there
  if("ec" %in% colnames(mod@react_attr))
    dt$ec <- mod@react_attr$ec
  
  #dt <- cbind(dt, data.table(mod@react_attr))
  #colnames(dt)[duplicated(colnames(dt))] <- paste0(colnames(dt)[duplicated(colnames(dt))],"_2")
  
  if(exclude.unused)
    dt <- dt[flux != 0]
  
  dt$equation <- printReact.2(mod, react = dt$rxn)
  return(dt)
}

get.utilized.metabolites <- function(mod) {
  
  # get MTF solution
  #rxn.coef <- ifelse(mod@react_attr$status %in% c("good_blast",NA,"no_seq_data"),1,1)
  sol.mtf <- optimizeProb(mod, algorithm = "mtf")
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf@fluxdist@fluxes[1:mod@react_num],
                        lb = mod@lowbnd)
  dt.mtf.tmp <- copy(dt.mtf[grepl("^EX_", ex)])
  dt.mtf.tmp[mtf.flux < 0, mtf.flux := 0]
  model.tmp <- changeBounds(mod, react = dt.mtf.tmp$ex, ub = dt.mtf.tmp$mtf.flux)
  #model.tmp <- mod
  
  # get FV solution
  sol.fv <- fluxVar(model.tmp, react = mod@react_id[grep("^EX_", mod@react_id)])
  
  dt <- data.table(ex = rep(mod@react_id[grep("^EX_", mod@react_id)],2),
                   rxn.name = rep(mod@react_name[grep("^EX_", mod@react_id)],2),
                   dir = c(rep("l",length(grep("^EX_", mod@react_id))),rep("u",length(grep("^EX_", mod@react_id)))),
                   fv = sol.fv@lp_obj)
  dt <- dcast(dt, ex + rxn.name ~ dir, value.var = "fv")[(l < -1e-4 & u <= 0) | (l < -1)]
  
  
  dt <- merge(dt, dt.mtf, by = "ex")
  dt[, flux.at.limit := ifelse(mtf.flux <= lb*0.999, "*","")]
  
  return(dt[order(l)])
}

get.produced.metabolites <- function(mod) {
  
  # get MTF solution
  #rxn.coef <- ifelse(mod@react_attr$status %in% c("good_blast",NA,"no_seq_data"),1,1)
  sol.mtf <- optimizeProb(mod, algorithm = "mtf")
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf@fluxdist@fluxes[1:mod@react_num])
  dt.mtf.tmp <- copy(dt.mtf[grepl("^EX_", ex)])
  
  # this following two lines are there to prevent the case that a nutrient (e.g. L-Lactate)
  # from the environment is taken up, and thus eneables the production of D-Lactate.
  dt.mtf.tmp[mtf.flux > 0, mtf.flux := 0]
  model.tmp <- changeBounds(mod, react = dt.mtf.tmp$ex, lb = dt.mtf.tmp$mtf.flux)
  
  
  # get FV solution
  sol.fv <- fluxVar(model.tmp, react = mod@react_id[grep("^EX_", mod@react_id)])
  
  dt <- data.table(ex = rep(mod@react_id[grep("^EX_", mod@react_id)],2),
                   rxn.name = rep(mod@react_name[grep("^EX_", mod@react_id)],2),
                   dir = c(rep("l",length(grep("^EX_", mod@react_id))),rep("u",length(grep("^EX_", mod@react_id)))),
                   fv = sol.fv@lp_obj)
  dt <- dcast(dt, ex + rxn.name ~ dir, value.var = "fv")[(u>1e-4 & l >= 0) | (u > 1)]
  
  
  
  dt <- merge(dt, dt.mtf, by = "ex")
  
  return(dt[order(-u)])
}


printReact.2 <- function(mod, react, printOut = FALSE, ...) {
  check <- checkReactId(mod, react = react)
  if (is(check, "reactId")) {
    cind <- react_pos(check)
  }
  else {
    stop("check argument react")
  }
  
  mat <- mod@S[, cind, drop = FALSE]
  nnz <- apply(mat, 2, "!=", 0)
  reaction <- character(length(cind))
  
  for (j in seq(along = cind)) {
    
    met <- met_name(mod)[nnz[, j]]
    nzv <- mat[, j][nnz[, j]]
    
    ed <- nzv < 0
    pd <- nzv > 0
    
    if (sum(ed) > 0) {
      educt   <- paste(paste("(", abs(nzv[ed]), ")", sep = ""),
                       met[ed], collapse = " + ")
    }
    else {
      educt = ""
    }
    
    if (sum(pd) > 0) {
      product <- paste(paste("(", nzv[pd], ")", sep = ""),
                       met[pd], collapse = " + ")
    }
    else {
      product = ""
    }
    
    #arrow   <- ifelse(react_rev(mod)[cind[j]], " <==> ", " --> ")
    arrow   <- ifelse(lowbnd(mod)[cind[j]] < 0 & uppbnd(mod)[cind[j]] > 0, " <==> ",
                      ifelse(lowbnd(mod)[cind[j]] >= 0," --> "," <-- "))
    
    reaction[j] <- paste(educt, product, sep = arrow)
  }
  
  return(reaction)
  
}


get.fluxvar.all.reactions <- function(mod) {
  sol <- sybil::optimizeProb(mod, algorithm = "mtf")
  
  dt <- data.table(rxn      = mod@react_id, 
                   name     = mod@react_name,
                   mtf.flux = sol@fluxdist@fluxes[1:mod@react_num])
  
  # Add EC information if there
  if("ec" %in% colnames(mod@react_attr))
    dt$ec <- mod@react_attr$ec
  
  
  dt$equation <- printReact.2(mod, react = mod@react_id)
  
  # get FV solution
  sol.fv <- fluxVar(mod, react = mod@react_id)
  
  dt.fv <- data.table(rxn  = rep(mod@react_id,2),
                      dir  = c(rep("l",length(mod@react_id)),rep("u",length(mod@react_id))),
                      fv   = sol.fv@lp_obj)
  dt.fv <- dcast(dt.fv, rxn ~ dir, value.var = "fv")
  
  dt <- merge(dt, dt.fv, by = "rxn")
  
  dt[abs(mtf.flux) < 1e-10, mtf.flux := 0]
  dt[abs(l)        < 1e-10, l        := 0]
  dt[abs(u)        < 1e-10, u        := 0]
  
  return(dt)
}
