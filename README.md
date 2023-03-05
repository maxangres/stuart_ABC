# Artificial Bee Algorithm for Psychometric Scale Shortening


This repository entails a version of the stuart R package written by Martin Schultze in collaboration with Johanna Schüller and myself. The Stuart (Subtests Using Algorithmic Rummaging Techniques) package offers a framework for metaheuristic algorithms which enable the user to construct subtests from psychometric scales. For an in-depth exploration of the package please visit the [package description](https://cran.r-project.org/web/packages/stuart/index.html) and the [manual](https://cran.r-project.org/web/packages/stuart/stuart.pdf).

The version included in this repository includes an Artificial Bee Colony Algorithm (ABC) which has not yet been merged with the main branch also accessible [here](https://bitbucket.org/martscht/stuart/src/7e21bcea11820b2924d4d3cdd8be802dddbebb6e/?at=feature%2FArtificialBees).

## Psychometric Scale Shortening

Automated item selection for scale shortening can be defined as a combinatorial problem by searching the optimal combination of items from all possibilities (Schultze, 2017). This definition enables the framing of item selection as a combinatorial optimization problem. Most psychometric scales consist of multiple facets measured by multiple scales included in a questionnaire. The framework for item selection as a combinatorial optimization problem needs to reflect this fact. Schultze (2017) found a suitable approach for the problem by representing it as an I-dimensional multiple knapsack problem with assignment restrictions (IMKAR).

In the IMKAR approach multiple knapsacks can depict multiple facets of a scale. Each facet gives a pool of items, from which items can be sampled until the weight boundary of the knapsack is violated by its components, i.e. the maximum number of items per subscale is surpassed. The weight of each item is solution unspecific, while the benefit of each item is solution specific. The benefit of each item is dependent on the other items selected in the solution due to the CFA-based approach. This makes the benefit both component and solution specific, i.e. the benefit of an item is dependent on both itself and the already selected other items in a solution. A component could therefore have a unique benefit in each knapsack, but is constrained in its selection to its specific knapsack, i.e. an item can only be selected into a solution for its original facet. This necessitates two further restrictions. First, each knapsack can only be assigned components that are eligible for assignment to it, i.e. items that are included in the facet represented by the knapsack. Second, the sets of components in two different knapsacks must be disjoint, i.e. no item can be selected into two different facets. 

Such a framework is fitting for item selection as it allows the selection of a whole scale from the pool of all possible scales, rather than selecting single items added to a scale based on individual merit alone. This enables the application for all settings occurring in item selection, like one- and multi-faceted questionnaires, multiple group-settings and longitudinal assessments, without formulating a specific strategy for each practical context (Schultze, 2017). Due to this framework metaheuristics are applicable for scale abridgement.

## Artificial Bee Colony Algorithm

The Artificial Bee Colony (ABC) algorithm is a metaheuristic optimization algorithm inspired by the foraging behavior of honeybees first developed by Karoboga (2005). The algorithm iteratively searches for the optimal solution to a given optimization problem by simulating the behavior of bees in a colony. The ABC algorithm can be used to solve I-dimensional multiple knapsack problems (Sundar et al., 2010). The main stages of the ABC algorithm for scale abridgement are:

**Initialization:** The algorithm begins by creating an initial population of solutions. Each solution represents a possible assignment of items to the knapsacks. The items are assigned randomly, and the solution is evaluated by calculating the total value of the assigned items for each knapsack and ensuring that the total weight of the assigned items does not exceed the knapsack capacity.

**Employed bee phase:** In this phase, each employed bee is randomly paired with another employed bee. Through a random recombination of the items included in their respective solutions a new solution is created. The fitness of the new solution is calculated via an objective function evaluating the solution on the basis of a combinatorial factor analysis (CFA). If the fitness of the new solution is better than the current position of the employed bee, the employed bee adopts the new solution. If the new solution is worse than the existing one another new solution is created from the recombination of the two paired employed bees. This process is repeated until either a better solution is found or a set number of maximum repetition is reached.

**Onlooker bee phase:** During this phase, onlooker bees select partner bees on the basis of their solution quality as evaluated by the objective function.  The probability of selecting a particular solution is proportional to its fitness value. The same process of recombination is repeated as in the employed bee phase.

**Scout bee phase:** If an employed bee exhausts its search without finding a better solution, it becomes a scout bee. A scout bee randomly generates a new solution in the search space and evaluates its fitness.

**Repair mechanism:** If during the onlooker bee phase a solution is not improved it will not be randomly replaced. Instead a CFA is again run on the solution to determine the worst fitting items which are subsequently eliminated from the solution and replaced by either randomly sampled items or from items from already existing solutions.

**Termination:** The algorithm terminates when a termination criterion is met, such as a maximum number of iterations or a target fitness value. The current solution is returned as the output of the algorithm.

An in-depth description of the algorithm can be found in the attached PDF from page 21 to 31.

##  ABC Function Integration in stuart

### Description
Construct subtests from a given pool of items using an artificial bee colony algorithm. Allows for multiple constructs, occasions, and groups.

### Usage

    abc(
	    data,
	    factor.structure, # data and structure
	    capacity=NULL,
	    item.weights=NULL,
	    item.invariance='congeneric',  #cross invariance
	    repeated.measures=NULL,
	    long.invariance='strict',
	    mtmm=NULL,
	    mtmm.invariance='configural',
	    grouping=NULL,
	    group.invariance='strict',
	    comparisons=NULL,
	    auxiliary=NULL,
	    use.order=FALSE,
	    software='lavaan',
	    cores=NULL,
	    objective=NULL,
	    ignore.errors=FALSE,
	    opt.lavaan.limit = TRUE,  #fitness specs
	    #algorithm specs
	    generations = 256,
	    individuals = 32,
	    limit = 10,
	    ls.sel = 'wide',
	    ls.onlooker.sel = 'narrow',
	    selection.pressure = NULL,
	    selection = "tournament",
	    scout.selection = "tournament",
	    cutoff = 0.5,
	    size_best_fs = 10,
	    size_rep_t = 2,
	    convergence.criterion = 'geno.between',
	    number_gen = 5,
	    tolerance = NULL,
	    schedule = 'run',
	    analysis.options=NULL,
	    suppress.model=FALSE,  #modeling specs
	    seed=NULL,
	    filename=NULL
    )

## Arguments

`data`  A data.frame containing all relevant data.  

`factor.structure` A list linking factors to items. The names of the list elements correspond to the factor names. Each list element must contain a character-vector of item names that are indicators of this factor.  

`capacity`  A list containing the number of items per subtest.  This must be in the same order as the  `factor.structure`  provided. If a single number, it is applied to all subtests. If  *NULL*  all items are evenly distributed among the subtests.

`item.weights`  A placeholder. Currently all weights are assumed to be one.  

`item.invariance` A character vector of length 1 or the same length as  factor.structure  containing the desired invariance levels between items pertaining to the same subtest. Currently there are five options: *congeneric*, *ess.equivalent*, *ess.parallel*,  *equivalent*, and *parallel*, the first being the default.  

`repeated.measures` A list linking factors that are repeated measures of each other. Repeated factors must be in one element of the list - other sets of factors in other elements of the list. When this is  *NULL*  (the default) a cross-sectional model is estimated.  

`long.invariance` A character vector of length 1 or the same length as  repeated.measures  containing the longitudinal invariance level of repeated items. Currently there are four options:  *configural*, *weak*, *strong*, and *strict*.  Defaults to *strict*.  When `repeated.measures`=*NULL*  this argument is ignored. 
 
`mtmm`  A list linking factors that are measurements of the same construct with different methods. Measurements of the same construct must be in one element of the list - other sets of methods in other elements of the list. When this is  *NULL*  (the  
default) a single method model is estimated.  

`mtmm.invariance` A character vector of length 1 or the same length as  mtmm  containing the invariance level of MTMM items. Currently there are five options: *none*, *configural*, *weak*, *strong*, and *strict*. Defaults to *configural*. With *none* differing items are allowed for different methods. When  `mtmm`=*NULL*  this argument is ignored.  

`grouping`  The name of the grouping variable. The grouping variable must be part of  data provided and must be a numeric variable.  

`group.invariance`  A single value describing the assumed invariance of items across groups. Currently there are four options: *configural*, *weak*, *strong*, and *strict*. Defaults to *strict*. When `grouping` = *NULL*  this argument is ignored.  

`comparisons`  A character vector containing any combination of *item*, *long*, *mtmm*, and *group* indicating which invariance should be assessed via model comparisons. The order of the vector dictates the sequence in which model comparisons are performed. Defaults to  *NULL*  meaning that no model comparisons are performed.  

`auxiliary`  The names of auxiliary variables in  data. These can be used in additional modeling steps that may be provided in  `analysis.options$model`.  

`use.order`  A logical indicating whether or not to take the selection order of the items into account. Defaults to *FALSE*.  

`software`  The name of the estimation software. Can currently be *lavaan* (the default) or  *Mplus*. Each option requires the software to be installed.  

`cores`  The number of cores to be used in parallel processing. If  *NULL*  (the default) the result of detectCores will be used. On Unix-y machines parallel processing is implemented via mclapply, on Windows machines it is realized via parLapply.  

`objective`  A function that converts the results of model estimation into a pheromone. 

`generations` Maximum number of generations to run. Defaults to 64.

`individuals` The number of individuals per generation. Defaults to 16.

`limit` A number adaptable to the state of the search that broadens or focuses the search space. Defaults to scheduled


`ls.sel` Determines whether new foodsource is selected for every limit iteration, can be *narrow* for ongoing item selection from first selected food source or *wide* for new selected solution for every iteration

`selection` Determines whether tournament or proportional selection is applied for partner selection in the onlooker bee phase

`ls.onlooker.sel` The selected method by which onlooker bees select there food sources. Can be either *proportional* for roulette random selection *tournament* (the default) for a semi-deterministic selection.

`scout.selection` Determines the kind of selection for repair operator

`cutoff` percentage of solution that will be eliminated in the repair operator

`size_best_fs` how many of the best solutions are saved

`size_rep_t size` of the tournament for repair items

`number_gen` number of generation over which scout bees need to replace whole population without improvement (maybe second parameter needed)

`selection.pressure` The pressure exerted during the selection process, depending on the \code{selection.pressure}: if *selection* = *proportional* the non-linearity coefficient of the pheromone when determining selection probability (the default is 1); if `selection` = *proportional* the number of randomly selected individuals from which to choose the best (the default is 5).

`convergence.criterion` The criterion by which convergence is determined. Can be one of four criteria *variance*, *median*, *geno.within*, and *geno.between*, *scout*, *none* . See 'details'.

`number_gen` number of generations over which a certain amount of scout iniatiations has to exceed threshold

`per_gen` number of scout initation per generation constituting the threshold

`tolerance` The tolerance for determining convergence. The default depends on the setting used for `convergence.criterion`. 

`schedule` The counter which the scheduling of parameters pertains to. Can be either *run* (the default), for a continuous schedule, *generation*, for a schedule that is restarted every time the population is reinitialized.

`analysis.options` A list additional arguments to be passed to the estimation software. The names of list elements must correspond to the arguments changed in the respective estimation software. E.g. `analysis.options$model` can contain additional modeling commands - such as regressions on auxiliary variables.

`suppress.model` A logical indicating whether to suppress the default model generation. If *TRUE* a model must be provided in `analysis.options$model`.

`seed` A random seed for the generation of random samples. 

`filename` The stem of the filenames used to save inputs, outputs, and data files when `software` = *Mplus*. This may include the file path. When *NULL* (the default) files will be saved to the temporary directory, which is deleted when the R session is ended.

### Details

The pheromone function provided via  objective  is used to assess the quality of the solutions.  These functions can contain any combination of the fit indices provided by the estimation software.  When using Mplus these fit indices are `rmsea`, `srmr`, `cfi`, `tli`, `chisq` (with `df` and ’pvalue’),  `aic`, `bic`, and `abic`. With lavaan any fit index provided by  inspect  can be used. Additionally `crel` provides an aggregate of composite reliabilites, `rel` provides a vector or a list of reliability  coefficients for the latent variables, `con` provides an aggregate consistency estimate for MTMM  analyses, and `lvcor` provides a list of the latent variable correlation matrices. For more detailed objective functions `lambda`, `theta`, `psi`, and `alpha` provide the model-implied matrices. Per  default a pheromone function using ’crel’, ’rmsea’, and ’srmr’ is used. Please be aware that the  objective  must be a function with the required fit indices as (correctly named) arguments. Using model comparisons via the  comparisons  argument compares the target model to a model  with one less degree of assumed invariance (e.g. if your target model contains strong invariance,  the comparison model contain weak invariance). Adding comparisons will change the preset for  the objective function to include model differences. With comparisons, a custom objective function  (the recommended approach) can also include all model fit indices with a preceding  delta.  To  indicate the difference in this index between the two models. If more than one type of comparison  is used, the argument of the objective function should end in the type of comparison requested (e.g. delta.cfi.group to use the difference in CFI between the model comparison of invariance across  groups).

The artificial bee algorithm offers multiple parameters with which the search process can be influenced. The three main parameters are `generations`, `individuals` and `limit`. `generations` determines after how many complete search iterations the algorithm terminates. `individuals` determines the total amount of bees employed in each bee phase. `limit` determines the number of tries a bee has to find a better solution.

The parameters `ls.sel` for the employed bee phase and ‘ls.onlooker.sel’ for the onlooker bee phase determine the pairing mechanism between bees during their intergenerational search defined by the `limit` parameter. A *wide* setting means that for each new limit iteration a new partner bee & respective solution is randomly selected. A *narrow* setting means that the same paring established during the first limit iteration is used for the sampling of new solutions during all limit iterations. The `selection` parameter determines how a partner bee is selected during the onlooker bee phase. The two possible mechanisms are `tournament` and `proportional` selection. Tournament selection lets a number of randomly created solutions compete against each other, while proportional selection selects a partner bee with a frequency directly proportional to the quality of its solution.

Three parameters `scout.selection`, `cutoff`, `size_best_fs` influence the repair operator. `scout.selection` determines the mechanism by which the worst items are eliminated from a non-improved onlooker bee solution. It can be set to *random*, *tournament* and *proportional*. `cutoff` determines the percentage of low-performing items eliminated from the solution. `size_best_fs` determines the number of high quality solution to sample items from if either tournament or proportional scout.selection is selected.

## Code examples

    # ABC selection in a simple situation  
    # requires lavaan  
    # number of cores set to 1 in all examples  
    data(fairplayer)  
    fs <- list(si = names(fairplayer)[83:92])  
    # minimal example  
    sel <- abc(fairplayer, fs, 4,  
    generations = 1, individuals = 10,  # minimal runtime, remove for application  
    seed = 55635, cores = 1)  
    summary(sel)  
    # longitudinal example  
    data(fairplayer)  
    fs <- list(si1 = names(fairplayer)[83:92],  
    si2 = names(fairplayer)[93:102],  
    si3 = names(fairplayer)[103:112])  
    repe <- list(si = c('si1',  'si2',  'si3'))  
     run to convergence  
    switching to best-last mating and 50\% mating size  
    sel <- abc(fairplayer, fs, 4,  
    repeated.measures = repe, long.invariance =  'strong',  
    seed = 55635, cores = 1)
