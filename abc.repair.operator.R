abc.repair.operator <-
  function(
    run, software, cores, factor.structure, short, short.factor.structure, 
    capacity, scout.selection, cutoff_cur, run.options, onlooker.repair, 
    sol.repair, cut_items, inter_gen, good_items, best_items_freq, size_rep_t 
    #, ... #iteration_limit = NA, ...
  ) { #begin function
    
    ## get items from solutions to be repaired
    selected <- list()
    for (j in 1:length(sol.repair)) {
      selected[[j]] <- as.numeric(inter_gen[[j]][run,])
    }
    
    ## create lavaan input to get modindices for selected solution
    selected.items <- translate.selection(selected,factor.structure,short)
    run.options$selected <- selected
    run.options$selected.items <- selected.items
    run.options$output.model <- TRUE
    ## call lavaan cfa calculation
    old_sol <- do.call("run.lavaan", run.options)
    ## get items from solution
    items <- run.options$selected
    
    ## get modification indices
    if (length(short.factor.structure) > 1){
      mi <- lavaan::modindices(old_sol$model)
      mi <- mi[mi$op == '=~', ]
    } else if (length(short.factor.structure) == 1){
      mi <- lavaan::modindices(old_sol$model)
    }
  
    ## fill list with items with highest modindices until cutoff is reached
    worst_items <- list()
    while(length(worst_items) < cut_items){
      wi <- which.max(tapply(mi$mi, mi$rhs, sum))
      worst_items <- append(worst_items, names(wi))
      mi <- mi[!mi$rhs == names(wi),]
    }
    ## get index of worst items by name & delete them
    for (l in 1:length(selected.items)){
      items <- selected.items[[l]]
      wi_indices <- match(items, worst_items)
      wi_indices <- which(is.na(wi_indices))
      selected[[l]] <- selected[[l]][c(wi_indices)]
    }
    
    ## create new matrix for new sols
    inter_gen <- vector("list", length(short.factor.structure))
    new_foodsources <- matrix(ncol = capacity[[1]]) #, nrow = capacity[[1]],ncol = individuals
    for (i in 1:length(short.factor.structure)){
      inter_gen[[i]] <- new_foodsources
    }
    
    ## fill solution with random items / items from best solutions depending on item selection 
    for (l in 1:length(short.factor.structure)){
      items <- selected[[l]]
      if (scout.selection == "random"){  
        while (length(items) < capacity[[1]]){
          items <- append(items, sample(1:length(factor.structure[[l]]), 1))
          items <- unique(items)
        }    
        selected[[l]] <- items
      } else if (scout.selection == "proportional"){
        while (length(items) < capacity[[1]]){
          if (length(unique(good_items[[l]])) == capacity[[1]]){
            items <- append(items, sample(1:length(factor.structure[[l]]), 1))
            items <- unique(items)
          } else {
            items <- append(items, sample(good_items[[l]], 1))
            items <- unique(items)
          }
        }
        selected[[l]] <- items
      } else if (scout.selection == "tournament"){
        while (length(items) < capacity[[1]]){
          if (length(unique(good_items[[l]])) == capacity[[1]]){
            items <- append(items, sample(1:length(factor.structure[[l]]), 1))
            items <- unique(items)
          } else {
            pre_sel <- best_items_freq[[l]][sample(nrow(best_items_freq[[l]]), size_rep_t), ]
            pre_sel <- pre_sel[which(pre_sel[,2] == max(pre_sel[,2])),1]
            new_item <- sample(as.numeric(as.character(pre_sel)), 1) #which.max(pre_sel[,2]),1]
            items <- append(items, new_item) #as.numeric(as.character(new_item))
            items <- unique(items)
          }
        }
        selected[[l]] <- items
      }
      inter_gen[[l]] <- rbind(c(selected[[l]]), inter_gen[[l]])
    }
    
    return(inter_gen)
  }