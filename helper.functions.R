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

################################################
# Function: dynamicFBA
#
# Performs a dynamic flux balance analysis
# 
# The function dynamicFBA() is inspired by the function
# dynamicFBA() contained in the COBRA Toolbox.
# The algorithm is the same.

dynamicFBA <- function (model,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,exclUptakeRxns,
                        retOptSol = TRUE,
                        fld = FALSE,verboseMode = 2, ...){
  #PARAMETERS:
  #===========
  # model                 Sybil model structure (class modelorg)
  # substrateRxns         List of exchange reaction names for substrates
  #                       initially in the media that may change (e.g. not
  #                       h2o or co2)
  # initConcentrations    Initial concentrations of substrates (in the same
  #                       structure as substrateRxns)
  # initBiomass           Initial biomass (must be non zero)
  # timeStep              Time step size
  # nSteps                Maximum number of time steps
  # fld                   indicates if all fluxes at all steps will be returned.
  # retOptSol             indicates if optsol calss will be returned or simple list
  #
  #OPTIONAL PARAMETERS
  #===================
  # exclUptakeRxns        List of uptake reactions whose substrate
  #                       concentrations do not change (Default =
  #                       {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'})
  # 
  #RETURN VALUES:
  #=============
  # concentrationMatrix   Matrix of extracellular metabolite concentrations
  # excRxnNames           Names of exchange reactions for the EC metabolites
  # timeVec               Vector of time points
  # biomassVec            Vector of biomass values
  # all_fluxes            Matrix containing the fluxes of all reactions at different steps
  
  #
  # If no initial concentration is given for a substrate that has an open
  # uptake in the model (i.e. model.lb < 0) the concentration is assumed to
  # be high enough to not be limiting. If the uptake rate for a nutrient is
  # calculated to exceed the maximum uptake rate for that nutrient specified
  # in the model and the max uptake rate specified is > 0, the maximum uptake 
  # rate specified in the model is used instead of the calculated uptake
  # rate.
  
  # reorder substrate ids and their cocentrations
  subrxnIDorder <- order(match(substrateRxns, react_id(model)))
  substrateRxns <- substrateRxns[subrxnIDorder]
  initConcentrations <- initConcentrations[subrxnIDorder]
  
  ##--------------------------------------------------------------------------##
  # check prerequisites 
  if (!is(model, "modelorg")) {
    stop("needs an object of class modelorg!")
  }
  ##--------------------------------------------------------------------------##
  
  # Uptake reactions whose substrate concentrations do not change
  if (missing(exclUptakeRxns)){
    exclUptakeRxns = c('EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)');
    if (verboseMode > 2){
      print('Default extra cellular uptake reactions will be used: ')
      print(exclUptakeRxns);
    }
  }
  
  # Find exchange reactions
  excReact = findExchReact(model);
  excReactInd=(react_id(model) %in% react_id(excReact));#excReact$exchange
  #represent extra cellular reaction with boolean vector.
  exclUptakeRxnsInd=is.element(react_id(model) ,exclUptakeRxns);
  #Exclude reactions with concentrations that will not be changed 
  excReactInd = excReactInd & !exclUptakeRxnsInd;   #excInd & ~ismember(model.rxns,exclUptakeRxns);
  #get reaction names
  excRxnNames =react_id(model)[excReactInd];                #excRxnNames = model.rxns(excInd);
  
  substrateRxnsInd=(react_id(model) %in% substrateRxns)
  # Figure out if substrate reactions are correct: all substrate reactions should be exchange reactions.
  missingSub = substrateRxnsInd & !excReactInd;
  if (sum(missingSub)!=0){
    print(sum(missingSub));
    print(react_id(model)[missingSub]);
    print('Invalid substrate uptake reaction!');
  }
  ## 	***********************************************************     ##
  # Initialize concentrations
  #substrateMatchInd = intersect(excRxnNames,substrateRxns);
  concentrations=rep(0,length(react_id(model)))#table(excRxnNames);##vector(length=length(excRxnNames),mode="numeric");
  #concentrations[1:length(concentrations)]=0;
  concentrations[substrateRxnsInd] = initConcentrations;
  
  # Deal with reactions for which there are no initial concentrations
  originalBound = -lowbnd(model);# take all to be able to directly update
  noInitConcentration = (concentrations==0)&(lowbnd(model)<0)#(concentrations == 0 & originalBound > 0);
  concentrations[noInitConcentration] = 1000;
  
  biomass = initBiomass;
  
  # Initialize bounds
  uptakeBound =  concentrations/(biomass*timeStep);
  
  # Make sure bounds are not higher than what are specified in the model
  aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
  uptakeBound[aboveOriginal] = originalBound[aboveOriginal];
  lowbnd(model)[excReactInd]  = -uptakeBound[excReactInd];
  
  concentrationMatrix = concentrations[excReactInd];
  biomassVec = biomass;
  timeVec = 0;
  ##------------------------------------- Prepare Problem object --------------------##
  # get OptObj instance
  #lpmod <- prepProbObj(model,
  #                         nCols      = react_num(model),
  #                         nRows      = met_num(model),
  #             #            alg        = "FBA",
  #                         solver     = solver,
  #                         method     = method,
  #                         lpdir      = lpdir
  #                         #solverParm = solverParm
  #             )
  
  lpmod <- sybil::sysBiolAlg(model, algorithm = "fba", ...)
  
  ##-----------------------------------------------------------------------------##
  if (verboseMode > 2) print('Step number    Biomass\n');
  # Inititialize progress bar ...');
  #if (verboseMode == 2)  progr <- .progressBar();
  all_fluxes=NULL;
  all_stat=NULL;
  
  for (stepNo in 1:nSteps){
    
    #    if (verboseMode == 2)  progr <- .progressBar(stepNo, nSteps, progr);
    
    # Run FBA
    sol = sybil::optimizeProb(lpmod);
    mu =  sol$obj;  ##objvalue sol.f
    if ( length(checkSolStat(sol$stat,solver(problem(lpmod))))!=0 ){## checkSolStat
      print('No feasible solution - nutrients exhausted\n');
      break;
    }
    all_stat = c(all_stat,sol$stat)
    uptakeFlux = sol$fluxes[excReactInd];
    biomass = biomass*exp(mu*timeStep);
    #biomass = biomass*(1+mu*timeStep);
    biomassVec = c(biomassVec,biomass);
    if(fld){
      if (stepNo == 1) {
        all_fluxes = sol$fluxes;
      }else{
        all_fluxes = cbind(all_fluxes,sol$fluxes);
      }
    }
    # Update concentrations
    concentrations[excReactInd]= concentrations[excReactInd] - uptakeFlux/mu*biomass*(1-exp(mu*timeStep));
    #concentrations = concentrations + uptakeFlux*biomass*timeStep;
    concentrations[concentrations <= 0] = 0;
    concentrationMatrix = cbind(concentrationMatrix,concentrations[excReactInd]);
    
    # Update bounds for uptake reactions
    uptakeBound[excReactInd] =  concentrations[excReactInd]/(biomass*timeStep);
    # This is to avoid any numerical issues
    uptakeBound[uptakeBound > 1000] = 1000;
    # Figure out if the computed bounds were above the original bounds
    aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
    # Revert to original bounds if the rate was too high
    uptakeBound[aboveOriginal] = originalBound[aboveOriginal];# uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
    uptakeBound=ifelse(abs(uptakeBound) < 1e-9,0,uptakeBound);
    ## Change lower bounds according to the result of last step
    #lowbnd(model)[excReactInd]  = -uptakeBound[excReactInd];  
    uppb_tmp <- getColsUppBnds(problem(lpmod), which(excReactInd));
    changeColsBnds(problem(lpmod),which(excReactInd),lb=-uptakeBound[excReactInd],ub=uppb_tmp);
    
    if (verboseMode > 2) print(paste(stepNo,sep="    ",biomass));
    #waitbar(stepNo/nSteps,h);
    timeVec = c(timeVec,stepNo*timeStep);
  }# end loop
  
  # browser();
  row.names(concentrationMatrix)=react_id(model)[excReactInd];
  ## Preparing OUTPUT
  #concentrationMatrix,excRxnNames,timeVec,biomassVec
  if (isTRUE(retOptSol)) {
    if(is.null(all_fluxes)) all_fluxes=as.matrix(NA);
    return (optsol_dynamicFBA(solver = solver(problem(lpmod)),
                              method = method(problem(lpmod)),
                              nprob  = stepNo,
                              ncols  = react_num(model),
                              nrows  = met_num(model),
                              fld    = fld,
                              all_fluxes = all_fluxes,
                              concmat=concentrationMatrix,
                              exRxn=excRxnNames,
                              tmVec=timeVec,  
                              bmVec=biomassVec
    )
    )
  }else{
    return(optsol <- list(		  nprob  = stepNo,
                              ncols  = react_num(model),
                              nrows  = met_num(model),
                              all_fluxes = all_fluxes,
                              all_stat=all_stat,
                              concentrationMatrix=concentrationMatrix,
                              excRxnNames=excRxnNames,
                              timeVec=timeVec,  
                              biomassVec=biomassVec
    ))
  }
}



# 1/6- Definition: to be added to file AllClasses.R

#------------------------------------------------------------------------------#
#                  definition of the class optsol_dynamicFBA                   #
#------------------------------------------------------------------------------#

setClass("optsol_dynamicFBA",
         representation(
           concentrationMatrix="matrix", # Matrix of extracellular metabolite concentrations
           excRxnNames="character",       # Names of exchange reactions for the EC metabolites
           timeVec="numeric",             # Vector of time points
           biomassVec="numeric"           # Vector of biomass values
           , all_fluxes="matrix" 				# Matrix of all fluxes at all steps   24/7/2015
         ),
         contains = "optsol_optimizeProb",
         package = "sybil"
)

#showClass("optsol_dynamicFBA")


#-------------------------------------------------------------------------------#

# 2/6-Constructor: to be added to file AllClasses-constructors.R
# optsol_dynamicFBAClass
optsol_dynamicFBA <- function(solver, method, nprob,
                              #lpdir,
                              ncols, nrows, 
                              #objf,
                              fld,concmat,exRxn,tmVec,bmVec,all_fluxes) {
  if (missing(solver) || 
      missing(method) ||
      missing(nprob)  ||
      #missing(lpdir)  ||
      missing(ncols)  ||
      missing(nrows)  ||
      #missing(objf)   ||
      missing(fld)    ||
      missing(bmVec) ||
      missing(tmVec)  ||
      missing(all_fluxes)
  ) {
    stop("Not enough arguments for creating an object of class optsol_dynamicFBA!")
  }
  
  if (fld == TRUE) {
    fldist <- fluxDistribution(all_fluxes, ncols, nprob)
  }
  else {
    fldist <- fluxDistribution(NA)
  }
  
  new("optsol_dynamicFBA", 
      solver       = as.character(solver),
      method       = as.character(method),
      num_of_prob  = as.integer(nprob),
      lp_num_cols  = as.integer(ncols),
      lp_num_rows  = as.integer(nrows),
      lp_obj       = numeric(nprob),
      lp_ok        = integer(nprob),
      lp_stat      = integer(nprob),
      #lp_dir       = as.character(lpdir),
      #obj_function = as.character(objf),
      fluxdist     = fldist,
      #   num_of_steps_executed=nsteps,
      concentrationMatrix=concmat,
      excRxnNames=exRxn,
      timeVec= tmVec,  
      biomassVec = bmVec,
      all_fluxes = all_fluxes
  )
}


#-----------------------------------------------------------------------------#

# 3/6-dynamicFBA: dynamicFBA.R
##              3.1 Check
##              3.2 optimizeProb -> get OptObj
##              3.3 Set Bounds
##              3.4 Call SimpleFBA
##              3.5 Store Solution
##              3.6 Adjust OUTPUT to class optsol_dynamicFBA

#-----------------------------------------------------------------------------#

# 4/6-  accessors: to be added to file optsol_dynamicFBA-accessors.R

#-----------------------------------------------------------------------------#

# 5/6-plot: to be added to file plot-methods.R

# optsol_dynamicFBAClass
setMethod("plot", signature("optsol_dynamicFBA","missing"),
          function(x,y,
                   ylim=50,
                   xlab = "",
                   ylab = "Value",
                   type = "p",
                   pch = 20,
                   col = "black",             
                   #               collower, colupper, pchupper, pchlower,
                   #               dottedline = TRUE,
                   plotRxns=NULL,
                   baseline = 0,
                   ...) {
            if(missing(plotRxns)){
              plot(x@timeVec,x@biomassVec,main='Biomass',xlab='Time',ylab=ylab);
            }
            else {
              def.par <- par(no.readonly = TRUE);
              layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE))
              #layout.show(2);
              # first plot biomass
              plot(spline(x@timeVec,x@biomassVec, n = 201, method = "natural"), col = 1
                   ,main='Cell density',xlab='Time(hrs)',ylab="X(g/l)",type="l",lwd=2);
              points(x@timeVec,x@biomassVec, col = "red",lwd=2);
              
              # plot concentrations
              ##plot(x@timeVec,2*x@biomassVec,main='Biomass',xlab='Time',ylab=ylab);
              ## define min/max ()plot(x@timeVec,
              ymin <- min(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
              ymax <- max(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
              for ( i in 1:length(plotRxns) ){
                plotInd=(x@excRxnNames %in% plotRxns[i]);
                #print( x@concentrationMatrix[plotInd]);
                if(i==1){
                  plot(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"),
                       type="l", col =i,main="Concentrations",ylab = "mmol",xlab='Time(hrs)'
                       ,ylim=c(ymin,ymax));
                }
                else{
                  lines(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"), col =i);
                }
              }
              legend(0,ymax, plotRxns, col=1:length(plotRxns), lty=1, bg = NA, cex = 0.75);
            }
            
            #if (!missing(plotRxns)){		   }
          }
)


#-----------------------------------------------------------------------------#

# 6/6-Documentation: optsol_dynamicFBA.Rd


#New version-----------------------------
