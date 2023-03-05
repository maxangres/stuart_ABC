call.bf.cycle <-
  function(
    software, cores, filter, bf.args
  ) { #begin function

    if (nrow(filter) > 0) {
      #parallel processing for R-internal estimations
      if (software=='lavaan') {
        if (cores>1) {
          #set up parallel processing on windows
          if (grepl('Windows',Sys.info()[1],ignore.case=TRUE)) {
            cl <- parallel::makeCluster(cores)
            #parallel::clusterExport(cl = cl, bf.args)
            bf.results <- parallel::parLapply(cl,1:nrow(filter),function(run) {
              do.call('bf.cycle',c(run,bf.args))
            })
            parallel::stopCluster(cl)
          }
          
          #run ants in parallel on unixies
          else {
            bf.results <- parallel::mclapply(1:nrow(filter), function(run) {
              do.call('bf.cycle',c(run,bf.args))},
              mc.cores=cores
            )
          }
        } else {
          bf.results <- lapply(1:nrow(filter),function(run) {
            do.call('bf.cycle',c(run,bf.args))})
        }
      }
      
      #serial processing if Mplus is used (Mplus-internal parallelization is used)
      if (software=='Mplus') {
        bf.args$filename <- filename
        bf.args$cores <- cores
        bf.results <- lapply(1:nrow(filter), function(run) {     
          do.call('bf.cycle',c(run,bf.args))})
      }
      return(bf.results)
    } else if (nrow(filter) == 0){
      bf.results <- list()
    }
  }

