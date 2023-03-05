stuart.abc <-
  function(
    short.factor.structure, short, long.equal, comparisons.equal,
    comparisons.invariance, #made on toplevel
    capacity,
    data, factor.structure, auxi, use.order,                             #simple prerequisites
    item.invariance,  
    repeated.measures, long.invariance,                                  #longitudinal relations
    mtmm, mtmm.invariance,                                               #mtmm relations
    grouping, group.invariance,                                          #grouping relations
    comparisons,
    software, cores,                                                     #Software to be used
    
    objective=NULL, ignore.errors=FALSE, opt.lavaan.limit = TRUE,       #objective function
    
    generations = 256, individuals = 32, limit = 10,                                  #algorithm specs
    ls.sel = 'wide', ls.onlooker.sel = 'narrow', selection.pressure = NULL,
    selection = "tournament", scout.selection = "tournament", cutoff = 0.5, 
    size_best_fs = 10, size_rep_t = 2, number_gen = 5, convergence.criterion = 'geno.between', 
    tolerance = NULL,
    
    schedule = 'run',
    
    suppress.model=FALSE, analysis.options=NULL,                   #Additional modeling
    seed,
    filename,
    
    ...                                                            #All the other stuff
  ) { #function begin
    
    #set random seed, if provided
    if (!is.null(seed)) {
      old.seed <- .Random.seed
      old.kind <- RNGkind()[1]
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)
    }
    
    #initialize fitness results of best solutions
    phe.ib <- 0
    phe.gb <- 0
    
    
    # bind dependent parameters together
    if (is.null(selection.pressure)) {
      selection.pressure <- selection
      selection.pressure[selection == 'tournament'] <- 2
      selection.pressure[selection == 'proportional'] <- 1
      if (is.matrix(selection.pressure)) {
        selection.pressure <- matrix(as.numeric(selection.pressure), ncol = ncol(selection.pressure))
      } else {
        selection.pressure <- as.numeric(selection.pressure)
      }
    }
    
    # decompress tolerances (for multiple convergence criteria)
    # explicit assignment to avoid check note
    if (is.null(tolerance)) tolerance <- vector('list', length(convergence.criterion))
    if (!is.list(tolerance)) tolerance <- list(tolerance)
    if (is.null(names(tolerance))) names(tolerance) <- convergence.criterion
    tolerance_va <- tolerance[['variance']]
    tolerance_md <- tolerance[['median']]
    tolerance_gw <- tolerance[['geno.within']]
    tolerance_gb <- tolerance[['geno.between']]
    
    # tolerance presets
    tols <- ls(pattern = 'tolerance_')
    pres <- c(.05, .7, .10, .005,
              .01, .8, .05, .0005)
    for (i in seq_along(tols)) {
      if (is.null(get(tols[i]))) assign(tols[i], pres[i])
    }
    
    #initialize scheduling
    scheduled <- c('generations', 'individuals', 'limit', 'selection', 'selection.pressure', 'size_best_fs', 'cutoff', 'tolerance_va', 'tolerance_md', 'tolerance_gw', 'tolerance_gb')
    
    #global assignment to avoid check note
    iteration_limit <- NA
    generations_cur <- NA 
    individuals_cur <- NA
    limits_cur <- NA 
    size_best_fs_cur <- NA
    cutoff_cur <- NA
    selection_cur <- NA 
    selection.pressure_cur <- NA 
    reproduction_cur <- NA 
    tolerance_va_cur <- NA 
    tolerance_md_cur <- NA 
    tolerance_gw_cur <- NA 
    tolerance_gb_cur <- NA 
    
    if (opt.lavaan.limit == TRUE){
      cfa_info <- "iterations" #technicals <- "iterations"
    } else {
      cfa_info <- NULL #technicals <- NULL
    }
    
    filt <- sapply(mget(scheduled),is.array)
    for (i in 1:length(scheduled[!filt])) {
      assign(paste0(scheduled[!filt][i],'_cur'),mget(scheduled[!filt][i])[[1]])
    }
    if (length(scheduled[filt])>0) {
      scheduled <- scheduled[filt]
      for (i in 1:length(scheduled)) {
        tmp <- mget(scheduled[i])[[1]]
        if (!any(c(0,1)%in%tmp[,1])) {
          stop(paste('The parameter schedule for',scheduled[i],'does not contain a value for the first generation.'),call.=FALSE)
        }
        tmp <- tmp[which.min(tmp[,1]),2]
        assign(paste0(scheduled[i],'_cur'),tmp)
      }
    } else {
      scheduled <- NULL
    }
    
    #counting
    generation <- 1
    run <- 1
    
    # genotypes
    geno <- list()
    
    # scouts
    scout_log <- list()
    
    ### Loops ###
    log <- list()
    qual.ib <- NULL
    
    #creating user feedback
    message('Running STUART with Artificial Bee Colony Algorithm.\n')
    progress <- utils::txtProgressBar(0,max(c(generations_cur,1)),style=3)
    
    # generate initial population
    full <- FALSE
    n <- individuals_cur
    combinations <- do.call('generate.combinations',mget(names(formals(generate.combinations))))
    
    # prepare input to call bf.cycle
    call.bf.arg <- list()
    call.bf.arg$software <- software
    call.bf.arg$cores <- cores
    
    # create first template for duplicate & filter
    duplicate <- combinations$duplicate
    filter <- combinations$filter #[!duplicated(duplicate), , drop = FALSE]
    combi <- combinations$combi
    tried <- matrix(NA, ncol = sum(unlist(capacity)))[-1,]
    
    repeat { #over generations
      
      output.model <- FALSE
      svalues <- FALSE
      
      ## get arguments of bf.cycle function
      bf.args <- mget(names(formals(bf.cycle))[-1])
      
      ## inititate calculation of CFA limits if needed
      if (opt.lavaan.limit == TRUE){
        bf.args$cfa_info <- "iterations"   #$technicals <- "iterations"
        if (class(iteration_limit) == "numeric"){
          bf.args$iteration_limit <- iteration_limit
        }
      }

      ## create combination matrix
      combi_mat <- do.call(cbind, combi)
      call.bf.arg$filter <- filter
      call.bf.arg$bf.args <- bf.args
      
      ## create initial foodsources - one for each bee 
      foodsources <- as.list(seq(1, individuals_cur))
      foodsources_total <- foodsources
      
      ## calculate first created solutions
      if (run == 1){
        bf.results <- do.call("call.bf.cycle", call.bf.arg)
        tmp <- vector('list', individuals_cur)
        tmp[filter[,1]] <- bf.results
        bf.results <- tmp[duplicate]
        
        combi_eb <- combi
      }
      
      ## create measure for all tried solutions
      log <- c(log, lapply(bf.results, function(x) c(run = run, x)))
      tried <- rbind(tried, combi_mat)
      #bf.results <- lapply(bf.results, function(x) c(run = run, x))
      
      #parameter schedule
      if (!is.null(scheduled)) {
        for (i in 1:length(scheduled)) {
          tmp <- mget(scheduled[i])[[1]]
          if (schedule=='run') {
            if (any(tmp[,1]==run)) {
              message(paste0('Scheduled value of ',scheduled[i],' updated to ',tmp[which(tmp[,1]==run),2],'.'))
            }
            tmp <- tmp[max(which(tmp[,1]<=run)),2]
          } 
          if (schedule=='generation') {
            if (any(tmp[,1]==generation)) {
              message(paste0('Scheduled value of ',scheduled[i],' updated to ',tmp[which(tmp[,1]==generation),2],'.'))
            }
            tmp <- tmp[max(which(tmp[,1]<=generation)),2]
          }
          assign(paste0(scheduled[i],'_cur'),tmp)
        }
      }
      
      ##### components of artificial bee colony 
      ## create list for intergenerations
      inter_gen <- vector("list", length(short.factor.structure))
      ## create solution pair once
      if (ls.sel == "narrow") {
        partners <- vector("list", length(foodsources))
        partners <- sapply(partners, function(x) sample(foodsources, 1))
      }
      ## create list for improved items
      del_items <- list()

      #### employed bee phase
      for (lm in 1:limit_cur) {
        ## end loop if all solutions are improved
        if (length(foodsources) == length(del_items)) {
          break
        }
        ## check which food sources have been improved & exclude them from further search
        if (length(del_items) > 0){
          foodsources <- foodsources[-c(unlist(del_items))]
          del_items <- list()
        }
        ## create matrix to store all selected items for bf.cycle
        new_foodsources <- matrix(ncol = capacity[[1]]) #, nrow = capacity[[1]],ncol = individuals
        for (i in 1:length(inter_gen)){
          inter_gen[[i]] <- new_foodsources
        }
        ## select a partner bee for each solution at random
        if (ls.sel == "wide"){
          partners <- vector("list", length(foodsources))
          partners <- sapply(partners, function(x) sample(foodsources_total, 1))
          matches <- data.frame(unlist(foodsources), unlist(partners))
        } else if (ls.sel == "narrow"){
          if (lm == 1){
            matches <- data.frame(unlist(foodsources), unlist(partners))
          } else if (lm > 1){
            unimproved_rows <- list()
            for (i in 1:length(foodsources)){
              unimproved_rows <- append(unimproved_rows, which(matches[,1] == foodsources[[i]])) 
            }
            matches <- matches[c(unlist(unimproved_rows)),]
          }
        }
        
        ## recombine every solution and partner solution
        for (y in 1:nrow(matches)){
          ## check that no partnered solutions are identical
          while (matches[y,1][1] == matches[y,2][1]){
            matches[y,2][1] <- sample(foodsources_total, 1)
          }
          ## get items from both solution & partner solution
          foodsource <- lapply(combi_eb, function(x) x[matches[y,1],])
          partner <- lapply(combi_eb, function(x) x[matches[y,2],])
          ns <- vector("list", length(short.factor.structure))
          ## create an item pool for each scale with all items from both solutions 
          for (i in 1:length(short.factor.structure)) {
            fs_solution <- seq_along(short.factor.structure[[i]])%in%foodsource[[i]]
            p_solution <- seq_along(short.factor.structure[[i]])%in%partner[[i]]
            item_pool <- list()#vector("list", length(short.factor.structure)) #fs_solution - p_solution
            item_pool <- append(item_pool, which(fs_solution))
            item_pool <- append(item_pool, which(p_solution))
            ## create a list for new solution
            new_food <- list()
            ## sample individual items from all scale items in item pool
            while (length(new_food) < capacity[[1]]){
              n <- sample(item_pool, 1)
              new_food <- append(new_food, n) ## new food had index [[i]] maybe doesn't need one
              new_food <- unique(new_food)
            }
            ## translate solution back into true/false vector
            tmp <- as.list(seq(1, length(short.factor.structure[[i]])))
            ns[[i]] <- match(tmp, new_food)
            ns[[i]] <- !is.na(ns[[i]])
          }
          ## fill matrix with new solution items in next_gen
          for (i in 1:length(short.factor.structure)){
            inter_gen[[i]] <- rbind(which(ns[[i]]), inter_gen[[i]])
          }
        }
        ## integrate new solution into intergeneration
        inter_gen <- lapply(inter_gen, function(x) x[1:length(foodsources), ])
        for (i in 1:length(short.factor.structure)){
          cl <- class(inter_gen[[i]])
          # transform vector into matrix if only one solution is in intergeneration
          if (cl[[1]] == "integer"){
            inter_gen[[i]] <- t(matrix(inter_gen[[i]]))
          }
        }
        
        ## check for duplication in inter_gen
        duplicate <- match(data.frame(t(do.call(cbind, inter_gen))), data.frame(t(tried)))
        filter <- data.frame(matrix(which(is.na(duplicate)),
                                    nrow=sum(is.na(duplicate)),
                                    ncol=length(short.factor.structure)))
        ## fill in filter and new solutions for bf.cycle calculations
        bf.tmp.args <- bf.args
        bf.tmp.args[[1]] <- filter
        bf.tmp.args[[2]] <- inter_gen
        ## update arg & call bf.cycle for solution
        call.bf.arg$filter <- filter
        call.bf.arg$bf.args <- bf.tmp.args
        bf.tmp.results <- do.call("call.bf.cycle", call.bf.arg)
        
        ## save new tried & calculated solutions
        log <- c(log, lapply(bf.tmp.results, function(x) c(run = run, x)))
        inter_gen_mat <- do.call(cbind, inter_gen)
        tried <- rbind(tried, inter_gen_mat)
        ## fill in results for duplicates
        tmp <- vector('list', length(foodsources))
        tmp[filter[,1]] <- bf.tmp.results
        redo <- lapply(log[stats::na.omit(duplicate)], function(x) {
          if(all(is.na(x$solution.phe[-1]))) x$solution.phe$pheromone <- 0
          else x$solution.phe$pheromone <- do.call(objective$func, x$solution.phe[-c(1,length(x$solution.phe))])
          if(is.na(x$solution.phe$pheromone)) x$solution.phe$pheromone <- 0
          return(x)})
        tmp[sapply(tmp,is.null)] <- redo
        bf.tmp.results <- tmp
        
        ## compare newly created solution with old solution & exchange if new solution is superior
        for (i in 1:length(foodsources)){  
          index <- foodsources[[i]]
          current_sol <- bf.results[[index]]
          new_sol <- bf.tmp.results[[i]]
          if (current_sol$solution.phe$pheromone < new_sol$solution.phe$pheromone){
            bf.results[[index]] <- new_sol
            del_items <- append(del_items, foodsources[[i]])
          }
        }
      }
      
      #### Scout bees employed bee phase
      ## if solutions in employed bee phase were not improved, sample new ones in their stead
      if (length(foodsources) > 0) {
        bf.results <- bf.results[-c(unlist(foodsources))]
        bf.scout.results <- list()
        scout_pheromones <- list(0)
        ## sample new solutions till all solutions are valid
        while (0 %in% scout_pheromones == TRUE) { 
          ## create random solutions
          gen_comb_arg <- mget(names(formals(generate.combinations)))
          gen_comb_arg$n <- length(foodsources)
          scout_combinations <- do.call('generate.combinations', gen_comb_arg)
          scout_combi <- scout_combinations$combi
          ## integrate random solution into matrix for all scout solutions
          for (i in 1:length(short.factor.structure)){
            scout_combi[[i]] <- scout_combi[[i]][sample(1:nrow(scout_combi[[i]]), length(foodsources)), ]
            cl <- class(scout_combi[[i]])
            if (cl[[1]] == "integer"){
              scout_combi[[i]] <- t(matrix(scout_combi[[i]]))
            }
          }
          
          ## create new duplicates & filter for scout sols
          duplicate <- match(data.frame(t(do.call(cbind, scout_combi))), data.frame(t(tried)))
          filter <- data.frame(matrix(which(is.na(duplicate)),
                                      nrow=sum(is.na(duplicate)),
                                      ncol=length(short.factor.structure)))
          
          tmp <- do.call(cbind, scout_combi)
          if (anyDuplicated(tmp)) {
            duplicate <- match(data.frame(t(tmp)), data.frame(t(tmp)))
          } else {
            duplicate <- 1:nrow(tmp)
          }
          ## get arguments for bf.cycle and integrate filter & new solutions
          bf.tmp.args <- bf.args
          bf.tmp.args[[1]] <- filter
          bf.tmp.args[[2]] <- scout_combi
          ## limit core number to the number of evaluated solutions
          if (length(foodsources) <= cores){
            bf.tmp.args$cores <- length(foodsources)
          }
          ## call bf.cycle for solution
          call.bf.arg$filter <- filter
          call.bf.arg$bf.args <- bf.tmp.args
          bf.tmp.results <- do.call("call.bf.cycle", call.bf.arg)
          
          ## save new tried scout solutions
          log <- c(log, lapply(bf.scout.results, function(x) c(run = run, x)))
          scout_combi_mat <- do.call(cbind, scout_combi)
          tried <- rbind(tried, scout_combi_mat)
          ## fill in duplicate solutions
          tmp <- vector('list', length(foodsources))
          tmp[filter[,1]] <- bf.tmp.results
          redo <- lapply(log[stats::na.omit(duplicate)], function(x) {
            if(all(is.na(x$solution.phe[-1]))) x$solution.phe$pheromone <- 0
            else x$solution.phe$pheromone <- do.call(objective$func, x$solution.phe[-c(1,length(x$solution.phe))])
            if(is.na(x$solution.phe$pheromone)) x$solution.phe$pheromone <- 0
            return(x)})
          tmp[sapply(tmp,is.null)] <- redo
          bf.tmp.results <- tmp
          ## assign generation to solutions  
          bf.tmp.results <- lapply(bf.tmp.results, function(x){replace(x, x$run[1], run)})
          ## check if all solutions are valid
          scout_pheromones <- list()
          scout_pheromones <- sapply(bf.tmp.results, function(x) x$solution.phe$pheromone)
          if (0 %in% scout_pheromones == TRUE){
            bf.tmp.results <- bf.tmp.results[-which(scout_pheromones %in% 0)]
          }
          foodsources <- foodsources[which(scout_pheromones %in% 0)]
          bf.scout.results <- append(bf.scout.results, bf.tmp.results)
        }
        ## add new solutions to already established solutions
        bf.results <- append(bf.results, bf.scout.results)
      }
      ## create new combination matrix
      combi_eb <- vector('list', length(short.factor.structure))
      selected <- lapply(bf.results, function(x) x$selected)
      for (i in 1:length(selected)){
        for (l in 1:length(short.factor.structure)){
          combi_eb[[l]] <- rbind(combi_eb[[l]], selected[[i]][[l]])
        }
      }
    
      #### onlooker bee phase
      ## create combination matrix for onlooker bee phase
      if (run == 1){
        bf.onlooker.results <- bf.results
        combi_ob <- combi_eb
      }
      ## create lists for food sources and intergenerations
      foodsources <- as.list(seq(1, individuals_cur))
      inter_gen <- vector("list", length(short.factor.structure))
      ## get fitness values of solutions
      pheromones <- sapply(bf.onlooker.results, function(x) x$solution.phe$pheromone)
  
      ## create solution pair once by either proportionate or tournament selection
      if (ls.onlooker.sel == "narrow") {
        if (selection_cur == 'proportional') {
          partners_onlooker <- as.list(sample(1:length(pheromones), length(foodsources), TRUE, pheromones^selection.pressure_cur / sum(pheromones^selection.pressure_cur)))          
        }
        if (selection_cur == 'tournament') {
          partners_onlooker <- rep(NA, length(foodsources))
          pool <- 1:length(pheromones)
          for (i in 1:length(foodsources)) {
            if (length(pool) == 0) pool <- 1:length(pheromones)
            if (length(pool) < selection.pressure_cur) {
              tmp <- pool
            } else {
              tmp <- sample(pool, selection.pressure_cur)
            }
            partners_onlooker[i] <- tmp[which.max(pheromones[tmp])]
            pool <- pool[pool!=partners_onlooker[i]]
          }
        }
      }

      ## create list for improved solutions
      del_items <- list()
      ## for loop to calculate intergenerations
      for (lm in 1:limit_cur) {
        ## end loop if all sols are improved
        if (length(del_items) > 0){
          foodsources <- foodsources[-c(unlist(del_items))]
          del_items <- list()
        }
        ## end limit iteration search if all solutionas are imrpoved
        if (length(foodsources) == length(del_items)) {
          break
        }
        
        ## sample new partner solution each limit iteration through tournament or proportionate selection
        if (ls.onlooker.sel == "wide"){
          if (selection_cur == 'proportional') {
            partners_onlooker <- as.list(sample(1:length(pheromones), length(foodsources), TRUE, pheromones^selection.pressure_cur / sum(pheromones^selection.pressure_cur)))            
          }
          if (selection_cur == 'tournament') {
            partners_onlooker <- rep(NA, length(foodsources))
            pool <- 1:length(pheromones)
            for (i in 1:length(foodsources)) {
              if (length(pool) == 0) pool <- 1:length(pheromones)
              if (length(pool) < selection.pressure_cur) {
                tmp <- pool
              } else {
                tmp <- sample(pool, selection.pressure_cur)
              }
              partners_onlooker[i] <- tmp[which.max(pheromones[tmp])]
              pool <- pool[pool!=partners_onlooker[i]]
            }
          }
          matches <- data.frame(unlist(foodsources), unlist(partners_onlooker))
        } else if (ls.onlooker.sel == "narrow"){
          # exclude improved solutions from static solution pair list
          if (lm == 1){
            matches <- data.frame(unlist(foodsources), unlist(partners_onlooker))
          } else if (lm > 1){
            unimproved_rows <- list()
            for (i in 1:length(foodsources)){
              unimproved_rows <- append(unimproved_rows, which(matches[,1] == foodsources[[i]])) 
            }
            matches <- matches[c(unlist(unimproved_rows)),]
          }
        }
  
        ## create matrix to store all selected items for bf.cycle
        new_foodsources <- matrix(ncol = capacity[[1]]) #, nrow = capacity[[1]],ncol = individuals
        for (i in 1:length(inter_gen)){
          inter_gen[[i]] <- new_foodsources
        }
       
        ## recombine every solution and partner solution
        for (y in 1:nrow(matches)){
          ## check that no partnered solutions are identical
          while (matches[y,1][1] == matches[y,2][1]){
            matches[y,2][1] <- sample(foodsources_total, 1)
          }
          ## get items from both solution & partner solution
          foodsource <- lapply(combi_ob, function(x) x[matches[y,1],])
          partner <- lapply(combi_ob, function(x) x[matches[y,2],])
          new_sol <- vector("list", length(short.factor.structure))
          
          ## create an item pool for each scale with all items from both solutions
          for (i in 1:length(short.factor.structure)) {
            fs_solution <- seq_along(short.factor.structure[[i]])%in%foodsource[[i]]
            p_solution <- seq_along(short.factor.structure[[i]])%in%partner[[i]]
            item_pool <- list() #vector("list", length(short.factor.structure)) #fs_solution - p_solution
            item_pool <- append(item_pool, which(fs_solution))
            item_pool <- append(item_pool, which(p_solution))
            ## create a list for new solution
            new_food <- list()
            ## sample individual items from all scale items in item pool
            while (length(new_food) < capacity[[1]]){
              n <- sample(item_pool, 1)
              new_food <- append(new_food, n)
              new_food <- unique(new_food)
            }
            ## translate solution back into true/false vector
            tmp <- as.list(seq(1, length(short.factor.structure[[i]])))
            ns[[i]] <- match(tmp, new_food)
            ns[[i]] <- !is.na(ns[[i]])
          }
          ## fill matrix with new solution items in next_gen
          for (i in 1:length(short.factor.structure)){
            inter_gen[[i]] <- rbind(which(ns[[i]]), inter_gen[[i]])
          }
        }
        ## integrate new solution into intergeneration
        inter_gen <- lapply(inter_gen, function(x) x[1:length(foodsources), ])
        for (i in 1:length(short.factor.structure)){
          cl <- class(inter_gen[[i]])
          # transform vector into matrix if only one solution is in intergeneration
          if (cl[[1]] == "integer"){
            inter_gen[[i]] <- t(matrix(inter_gen[[i]]))
          }
        }
        
        ## check for duplication in inter_gen
        duplicate <- match(data.frame(t(do.call(cbind, inter_gen))), data.frame(t(tried)))
        filter <- data.frame(matrix(which(is.na(duplicate)),
                                    nrow=sum(is.na(duplicate)),
                                    ncol=length(short.factor.structure)))
        ## fill in filter and new solutions for bf.cycle calculations
        bf.tmp.args <- bf.args
        bf.tmp.args[[1]] <- filter
        bf.tmp.args[[2]] <- inter_gen
        ## set cores equal to solution count if smaller than total core count
        if (length(foodsources) <= cores){
          bf.tmp.args$cores <- length(foodsources)
        }
        ## update arg & call bf.cycle for solution
        call.bf.arg$filter <- filter
        call.bf.arg$bf.args <- bf.tmp.args
        bf.tmp.results <- do.call("call.bf.cycle", call.bf.arg)
        
        ## save new tried & calculated solutions
        log <- c(log, lapply(bf.tmp.results, function(x) c(run = run, x)))
        inter_gen_mat <- do.call(cbind, inter_gen)
        tried <- rbind(tried, inter_gen_mat)
        ## fill in results for duplicates
        tmp <- vector('list', length(foodsources))
        tmp[filter[,1]] <- bf.tmp.results
        redo <- lapply(log[stats::na.omit(duplicate)], function(x) {
          if(all(is.na(x$solution.phe[-1]))) x$solution.phe$pheromone <- 0
          else x$solution.phe$pheromone <- do.call(objective$func, x$solution.phe[-c(1,length(x$solution.phe))])
          if(is.na(x$solution.phe$pheromone)) x$solution.phe$pheromone <- 0
          return(x)})
        tmp[sapply(tmp,is.null)] <- redo
        bf.tmp.results <- tmp
        
        ## compare newly created solution with old solution & exchange if new solution is superior
        for (i in 1:length(foodsources)){
          index <- foodsources[[i]]
          current_sol <- bf.onlooker.results[[index]]
          new_sol <- bf.tmp.results[[i]]
          if (current_sol$solution.phe$pheromone < new_sol$solution.phe$pheromone){
            bf.onlooker.results[[i]] <- new_sol
            del_items <- append(del_items, foodsources[[i]])
          }
        }
      }
      
      #### scout phase onlooker bee phase
      ## if solutions in employed bee phase were not improved, sample new ones in their stead
      if (length(foodsources) > 0) {
        ## get items to be repaired
        onlooker.repair <- bf.onlooker.results[c(unlist(foodsources))]
        ## create new matrix for new sols
        inter_gen <- vector("list", length(short.factor.structure))
        new_foodsources <- matrix(ncol = capacity[[1]]) 
        for (i in 1:length(short.factor.structure)){
          inter_gen[[i]] <- new_foodsources
        }
        ## calculate number of items to be deleted
        cut_items <- round(capacity[[1]]*length(short.factor.structure)*cutoff_cur, 0)
        ## create lists for reduced items and new repaired items
        new_scout_sol <- list()
        old_sol <- list()
        
        #### Repair operator Item Pool
        ## select best solutions from log to give best items
        pheromones <- sapply(log, function(x) x$solution.phe$pheromone)
        best_selected <- with(log,order(-pheromones))
        best_selected <- best_selected[1:size_best_fs]
        best_selected <- log[c(best_selected)]
        best_selected <- lapply(best_selected, function(x) x$selected) #dataframe func
        ## create list with all items in best solutions & respective frequencies
        good_items <- list()
        best_items_freq <- list()
        best_items <- list()
        ## create seperate lists for tournament & proportionate item selection
        for (l in 1:length(short.factor.structure)){
          good_items[[l]] <- unlist(lapply(best_selected, function(x) x[l]))
          best_items_freq[[l]] <- as.data.frame(table(good_items[[l]]))
        }
        
        ## check if invalid sol is in onlooker bees 
        scout_pheromones <- list()
        scout_pheromones <- sapply(onlooker.repair, function(x) x$solution.phe$pheromone)
        if (0 %in% scout_pheromones == TRUE){
          # sample solution from employed bee phase if invalid solution in onlooker bee pop
          onlooker.repair[which(scout_pheromones %in% 0)] <- sample(bf.results, length(which(scout_pheromones %in% 0)))
        }
        
        ## call lavaan function with output for modification indices
        sol.repair <- vector("list", length(short.factor.structure))
        for (i in 1:length(short.factor.structure)){
          sol.repair[[i]] <- sapply(onlooker.repair, function(x) x$selected[[i]])
          #sol.repair[[i]] <- t(sol.repair[[i]])
        }
        ## create list for to be repaired solutions
        inter_gen <- vector("list", length(short.factor.structure))
        ## get items from solutions
        for (j in 1:length(foodsources)){
          for (i in 1:length(inter_gen)){
            inter_gen[[i]] <- rbind(inter_gen[[i]], as.numeric(sol.repair[[i]][, j]))
          }
        }
        ## transform vector into matrix for singular solutions
        for (i in 1:length(short.factor.structure)){
          cl <- class(inter_gen[[i]])
          if (cl[[1]] == "numeric"){
            inter_gen[[i]] <- t(matrix(inter_gen[[i]]))
          }
        }
        ## create list for retained items & get arguments for repair function
        selected <- list()
        selected.items <- list()
        run.options <- mget(names(formals(run.lavaan)))
        run.options$iteration_limit <- NULL
        rep.op.arg <- mget(names(formals(abc.repair.operator))[-1])
        
        ## parallel processing for repair 
        if (software=='lavaan') {
          if (cores>1) {
            #set up parallel processing on windows
            if (grepl('Windows',Sys.info()[1],ignore.case=TRUE)) {
              cl <- parallel::makeCluster(cores)
              parallel::clusterExport(cl = cl, 'abc.repair.operator')
              pre_inter_gen <- parallel::parLapply(cl,1:length(foodsources),function(run) {
                do.call('abc.repair.operator',c(run,rep.op.arg))
              })
              parallel::stopCluster(cl)
            }
            
            #run ants in parallel on unixies
            else {
              pre_inter_gen <- parallel::mclapply(1:length(foodsources), function(run) {
                do.call('abc.repair.operator',c(run,rep.op.arg))},
                mc.cores=cores, mc.preschedule = FALSE
              )
            }
          } else {
            pre_inter_gen <- lapply(1:length(foodsources),function(run) {
              do.call('abc.repair.operator',c(run, rep.op.arg))})
          }
        }
        
        #serial processing if Mplus is used (Mplus-internal parallelization is used)
        if (software=='Mplus') {
          bf.args$filename <- filename
          bf.args$cores <- cores
          pre_inter_gen <- lapply(1:length(foodsources), function(run) {     
            do.call('abc.repair.operator',c(run,rep.op.arg))})
        }
        ## get the reduced solution
        c_sol <- lapply(pre_inter_gen, function(x) length(x)==length(short.factor.structure))
        pre_inter_gen <- pre_inter_gen[which(unlist(c_sol))]
        ## transform repair function output into inter_gen framework
        inter_gen <- vector("list", length(short.factor.structure))
        for (j in 1:length(pre_inter_gen)){
          for (i in 1:length(short.factor.structure)){
            inter_gen[[i]] <- rbind(inter_gen[[i]], na.omit(pre_inter_gen[[j]][[i]])) #as.numeric(pre_inter_gen[[j]][[i]][1,])
          }
        }
        ## check if inter_gen is right format if not convert to matrix
        for (i in 1:length(short.factor.structure)){
          cl <- class(inter_gen[[i]])
          if (cl[[1]] == "numeric"){ ##maybe numeric?
            inter_gen[[i]] <- t(matrix(inter_gen[[i]]))
          }
        }
        
        ## create new duplicates & filter for repaired solutions
        duplicate <- match(data.frame(t(do.call(cbind, inter_gen))), data.frame(t(tried)))
        filter <- data.frame(matrix(which(is.na(duplicate)),
                                    nrow=sum(is.na(duplicate)),
                                    ncol=length(short.factor.structure)))
        bf.tmp.args <- bf.args
        bf.tmp.args[[1]] <- filter
        bf.tmp.args[[2]] <- inter_gen
        ## update arg & call bf.cycle for solution
        call.bf.arg$filter <- filter
        call.bf.arg$bf.args <- bf.tmp.args
        bf.tmp.results <- do.call("call.bf.cycle", call.bf.arg)
        while (length(which(lengths(bf.tmp.results) == 0)) >= 1) {
          bf.tmp.results <- do.call("call.bf.cycle", call.bf.arg)
        }
        
        ## save new calculated items
        log <- c(log, lapply(bf.tmp.results, function(x) c(run = run, x)))
        inter_gen_mat <- do.call(cbind, inter_gen)
        tried <- rbind(tried, inter_gen_mat)
        ## fill in already calculated items
        tmp <- vector('list', length(bf.tmp.results))
        tmp[filter[,1]] <- bf.tmp.results
        redo <- lapply(log[stats::na.omit(duplicate)], function(x) {
          if(all(is.na(x$solution.phe[-1]))) x$solution.phe$pheromone <- 0
          else x$solution.phe$pheromone <- do.call(objective$func, x$solution.phe[-c(1,length(x$solution.phe))])
          if(is.na(x$solution.phe$pheromone)) x$solution.phe$pheromone <- 0
          return(x)})
        tmp[sapply(tmp,is.null)] <- redo
        bf.tmp.results <- tmp
        ## add generation to repaired solution
        bf.tmp.results <- lapply(bf.tmp.results, function(x){replace(x, x$run[1], run)})
        
        ## add repaired solutions to onlooker bee population
        bf.onlooker.results <- bf.onlooker.results[-c(unlist(foodsources))]
        bf.onlooker.results <- append(bf.onlooker.results, bf.tmp.results)
        ## if lavaan calculation failed add solution from employed bee phase
        while (length(bf.onlooker.results) < individuals_cur){
          bf.onlooker.results <- append(bf.onlooker.results, sample(bf.results, 1))
        }
      } 
     
      ## replace the overall worst solution with the global best
      if (exists("solution.gb") == TRUE){
        scout_pheromones <- list()
        scout_pheromones <- sapply(bf.onlooker.results, function(x) x$solution.phe$pheromone)
        bf.onlooker.results[[which.min(scout_pheromones)]] <- list(selected = selected.gb, solution.phe = list(pheromone = phe.gb))
      }
      
      ## check whether solutions with invalid design are icluded & replace them 
      for (i in 1:length(bf.onlooker.results)){
        if ("selected" %in% names(bf.onlooker.results[[i]]) == FALSE) {
          bf.onlooker.results[i] <- sample(bf.results, 1)
        } else if (length(bf.onlooker.results[[i]]$selected) != length(short.factor.structure)) {
          bf.onlooker.results[i] <- sample(bf.results, 1)
        }
      }
      
      ## get new item combinations of employed bees
      combi_ob <- vector('list', length(short.factor.structure))
      selected <- lapply(bf.onlooker.results, function(x) x$selected)
      for (i in 1:length(selected)){
        for (l in 1:length(short.factor.structure)){
          combi_ob[[l]] <- rbind(combi_ob[[l]], selected[[i]][[l]])
        }
      }

      ## calculate limit on lavaan cfa iterations
      if (opt.lavaan.limit == TRUE){
        if (run == 1){
          ## get lists for both valid and invalid solutions
          non_converged <- lapply(log, function(x) x[x$solution.phe$cfa_info$converged == FALSE]) #solution.phe$technicals$converged
          converged <- lapply(log, function(x) x[x$solution.phe$cfa_info$converged == TRUE])
          iter_non_converged <- lapply(non_converged, function(x) x$solution.phe$cfa_info$converged)
          iter_non_converged <- iter_non_converged[lapply(iter_non_converged, length)>0]
          iter_converged <- lapply(converged, function(x) x$solution.phe$cfa_info$converged)
          iter_converged <- iter_converged[lapply(iter_converged, length)>0]
          ## calculate cluster to destinct between them
          if (length(iter_non_converged) > 0){
            border <- c(iter_non_converged, iter_converged)
            ks <- kmeans(border, 2)
            proper <- which.min(ks$centers)
            bp <- boxplot.stats(unlist(border[ks$cluster==proper]))
            iteration_limit <- bp$stats[which.max(bp$stats)]+round((bp$stats[which.max(bp$stats)]*0.5), 0) 
          }
        }
      }
      
      #iteration.best memory
      all_results <- append(bf.results, bf.onlooker.results)
      individual.ib <- which.max(sapply(all_results, function(x) return(x$solution.phe$pheromone)))
      solution.ib <- list()
      for (i in seq_along(short.factor.structure)) {
        solution.ib[[i]] <- seq_along(short.factor.structure[[i]])%in%all_results[[individual.ib]]$selected[[i]] ## names of items
      }
      
      phe.ib <- all_results[[individual.ib]]$solution.phe$pheromone
      selected.ib <- all_results[[individual.ib]]$selected # number of items
      
      #global.best memory
      if (phe.ib > phe.gb | generation == 1) {
        solution.gb <- solution.ib
        phe.gb <- phe.ib
        selected.gb <- selected.ib
      }
      
      #scout log
      tmp <- list(list( run = run, number_scouts = length(foodsources), scout_bees = unlist(foodsources), best_phe = phe.gb))
      scout_log <- append(scout_log, tmp)
      
      #create log for convergence criteria
      utils::setTxtProgressBar(progress,generation)
      
      # check for convergence
      if (generation >= generations_cur) {
        end.reason <- 'Maximum number of generations exceeded.'
        break
      }
      
      conv <- FALSE
      qual.ib <- c(qual.ib, phe.ib)
      
      if ('variance' %in% convergence.criterion) {
        if (generation > max(min(c(.1*generations_cur, 10)),1) & stats::var(qual.ib/qual.ib[1]) <= tolerance_va_cur) {
          conv <- TRUE
        }
      }
      
      if ('median' %in% convergence.criterion) {
        if (generation > max(min(c(.1*generations_cur, 10)),1) & ((phe.ib - stats::median(pheromones))/phe.ib) <= tolerance_md_cur) {
          conv <- TRUE
        }
      }
      
      if ('scout' %in% convergence.criterion) {
        number_gen <- 5 ### global parameter
        if ( run >= 10){
          sel <- tail(scout_log, n = number_gen)
          which_bees <- unique(unlist(sapply(sel, function(x) x$scout_bees)))
          bee_qual <- unique(unlist(sapply(sel, function(x) x$best_phe)))
          if (length(which_bees) == individuals_cur & length(bee_qual) == 1){
            conv <- TRUE
          }
        }
      }
      
      if ('geno.within' %in% convergence.criterion) {
        geno_var <- list()
        tmp <- c(0, cumsum(unlist(capacity)))
        for (i in 1:length(short.factor.structure)) {
          all_items <- seq_along(short.factor.structure[[i]])
          sel_items <- combi_mat[, ((tmp[i]+1):tmp[i+1])]
          geno[[i]] <- t(apply(sel_items, 1, function(x) all_items %in% x))
          geno_var[[i]] <- colMeans(geno[[i]])*(1-colMeans(geno[[i]]))
        }
        if ('geno.within' %in% convergence.criterion &
            all(unlist(sapply(geno_var, function(x) x < (tolerance_gw_cur*(1-tolerance_gw_cur)))))) {
          conv <- TRUE
        }
      }
      
      if ('geno.between' %in% convergence.criterion) {
        if (run == 1) {
          geno_r1 <- vector('list', length(short.factor.structure))
          geno_r1 <- lapply(geno_r1, function(x) 0)
          geno_d1 <- geno_r1
        } else {
          geno_r1 <- geno
        }
        geno_d2 <- geno_d1
        tmp <- c(0, cumsum(unlist(capacity)))
        for (i in 1:length(short.factor.structure)) {
          all_items <- seq_along(short.factor.structure[[i]])
          sel_items <- combi_mat[, ((tmp[i]+1):tmp[i+1])]
          geno[[i]] <- colMeans(t(apply(sel_items, 1, function(x) all_items %in% x)))
          geno_d1[[i]] <- abs(geno_r1[[i]]-geno[[i]])
        }
        
        if ('geno.between' %in% convergence.criterion &
            all(sapply(geno_d1, function(x) all(x < tolerance_gb_cur))) &
            all(sapply(geno_d2, function(x) all(x < tolerance_gb_cur)))) {
          conv <- TRUE
        }
      }
      
      
      #go on to next generation
      run <- run + 1
      generation <- generation + 1
    } # end loop
    
    #feedback
    close(progress)
    message(paste('\nSearch ended.',end.reason))      
    
    # reformat log
    #generate matrix output
    mat_fil <- c('lvcor', 'lambda', 'theta', 'psi', 'alpha', 'beta', 'nu')
    mat_fil <- mat_fil[mat_fil %in% names(formals(objective$func))]
    mats <- as.list(vector('numeric', length(mat_fil)))
    names(mats) <- mat_fil
    
    for (m in seq_along(mat_fil)) {
      mats[[m]] <- sapply(log, function(x) x$solution.phe[mat_fil[m]])
      names(mats[[m]]) <- 1:length(log)
    }
    
    # apply final pheromone function retroactively (empirical objectives)
    if (inherits(objective, 'stuartEmpiricalObjective')) {
      final_pheromone <- sapply(log, function(x) {
        if (x$solution.phe$pheromone == 0) 0
        else {do.call(objective$func, x$solution.phe[-1])}
      })
    }
    
    for (i in 1:length(log)){log[[i]]$solution.phe$technicals = NULL}
    # log <- sapply(log, function(x) list(x$run, x$selected, within(x$solution.phe, rm("technicals"))), USE.NAMES = TRUE)
    # names(log) <- c("run", "selected", "solution.phe")
    n <- names(log[[1]]$solution.phe)[!names(log[[1]]$solution.phe)%in%mat_fil]
    tmp <- sapply(log, `[[`, 1)
    tmp <- unlist(lapply(table(tmp), function(x) seq(1, x)))
    log <- cbind(cumsum(tmp == 1), tmp, t(sapply(log, function(x) array(data=unlist(x$solution.phe[!names(x$solution.phe)%in%mat_fil])))))
    log <- data.frame(log)
    names(log) <- c('run', 'ind', n)
    
    #return to previous random seeds
    if (!is.null(seed)) {
      RNGkind(old.kind)
      .Random.seed <<- old.seed
    }
    
    # Preparing output
    convergence <- vector('list', length(convergence.criterion))
    names(convergence) <- convergence.criterion
    if ('variance' %in% names(convergence)) convergence[['variance']] <- stats::var(qual.ib/qual.ib[1])
    if ('median' %in% names(convergence)) convergence[['median']] <- phe.ib - stats::median(pheromones)
    if ('geno.within' %in% names(convergence)) convergence[['geno.within']] <- lapply(geno, colMeans)
    if ('geno.between' %in% names(convergence) & 'geno_d1'%in%ls()) convergence[['geno.between']] <- geno_d1
    if ('scout' %in% names(convergence)) convergence[['scout']] <- tail(scout_log, n = number_gen)
    
    tolerance <- list(variance = tolerance_va, median = tolerance_md, 
                      geno.within = tolerance_gw, geno.between = tolerance_gb)
    
    #Generating Output
    tmp <- c(0, cumsum(unlist(capacity)))
    for (i in 1:length(short.factor.structure)) {
      all_items <- seq_along(short.factor.structure[[i]])
      sel_items <- combi_mat[, ((tmp[i]+1):tmp[i+1])]
      geno[[i]] <- colMeans(t(apply(sel_items, 1, function(x) all_items %in% x)))
    }
    names(geno) <- names(short.factor.structure)
    
    for (i in seq_along(solution.gb)) names(solution.gb[[i]]) <- short.factor.structure[[i]]
    names(solution.gb) <- names(short.factor.structure)
    
    results <- mget(grep('.gb',ls(),value=TRUE))
    results$selected.items <- translate.selection(selected.gb,factor.structure,short)
    results$log <- log
    results$log_mat <- mats
    results$pheromones <- pheromones
    results$parameters <- list(generations, individuals, selection, selection.pressure, ls.sel, ls.onlooker.sel,
                               selection, scout.selection, cutoff, 
                               size_best_fs, size_rep_t, number_gen, convergence.criterion, tolerance, 
                               seed, objective, factor.structure)
    names(results$parameters) <- c('generations', 'individuals', 'selection', 'selection.pressure', "ls.sel", "ls.onlooker.sel",
                                   "selection", "scout.selection", "cutoff", 
                                   "size_best_fs", "size_rep_t", "number_gen", 'convergence.criterion', "tolerance",
                                   'seed', 'objective', 'factor.structure')
    results$convergence <- convergence
    results$genotype <- geno
    results$end.reason <- end.reason
    return(results)
    
  }
