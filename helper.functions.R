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
