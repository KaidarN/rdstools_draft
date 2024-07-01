

#' Variance estimation with bootstrap chain and tree methods
#'
#' @param RESPONDENT_ID character or numeric vector; A variable indicating respondent ID
#' @param SEED_ID character or numeric vector; A variable indicating seed ID
#' @param SEED character or numeric vector; a variable indicating whether a particular respondent is seed or not
#' @param RECRUITER_ID character or numeric vector; A variable indicating recruiter ID
#' @param data.cov vector, matrix or data frame; Variable(s) of interest that will be used in the analysis
#' @param type character; One of the six types of bootstrap methods, please see the vignette for a more detailed information
#' @param n.times scalar; How many bootstrap resamples should be taken
#' @param return_data logical; By default return_data=F returns respondent ID and bootstrap iterations. If return_data TRUE, then the output will include variable(s) specified in data.cov
#'
#' @return When return_data=F returns a data frame consistiting of the following elements:
#' \item {RESPONDENT_ID}{character vector; A variable indicating respondent ID}
#' \item {boot_n}{character vector; A indicator variable for each bootstrap sample}
#'
#' @examples
#' data('RDStoydata')
#'
#' # Preprocess data with RDSdata function
#' rds_data <- RDSdata(data = RDStoydata,unique_id = "ID",
#' redeemed_coupon = "CouponR",
#' issued_coupon = c("Coupon1",
#'                  "Coupon2",
#'                  "Coupon3"),
#'                degree = "Degree",
#'                result = c('Age','Sex'))
#'
#'
#' # Run bootstrap_RDS with rds_data
#' results = bootstrap_RDS(RESPONDENT_ID = rds_data$ID, SEED_ID = rds_data$S_ID,
#' SEED = rds_data$SEED, RECRUITER_ID = rds_data$R_CP,
#' data.cov = rds_data[,c('Age', 'Sex')], type = 'boot_chain_one', n.times = 100,
#' return_data = T)
#'
#' @export
bootstrap_RDS = function(RESPONDENT_ID, SEED_ID, SEED, RECRUITER_ID,
                         data.cov, type, n.times, return_data=F){

  data.cov = as.data.frame(data.cov)
  data = as.data.frame(cbind(RESPONDENT_ID, SEED_ID, SEED, RECRUITER_ID, data.cov)) # combine data

  results_bootstrap_RDS = vector(mode='list', n.times)  # list to store all of the bootstrap iterations

  if (type == 'boot_chain_one'){
    ############################################# function for chain bootstrap 1
    boot_chain1 = function(data) {

      empty.list = list() # initialize an empty list
      df = data[which(data$SEED==1),] # only seeds are selected
      seed_id = sample(df$RESPONDENT_ID, nrow(df), replace=T) # all seeds are selected

      for (i in 1:length(seed_id)){
        empty.list[[i]] = data[data$SEED_ID==seed_id[i],]
      }

      data.df = do.call(rbind, empty.list) # data frame from do.call
      return(data.df)
    }

    for (i in seq_len(n.times)){
      results_bootstrap_RDS[[i]] = boot_chain1(data)
    }
  }

  else if (type == 'boot_chain_two'){
    ############################################# function for chain boostrap 2
    boot_chain2 = function(data) {

      bootstrap_list = list()
      seeds_df = data[data$SEED == 1, ]

      if (nrow(seeds_df) == 0) {
        stop("There are no records where SEED equals 1.")
      }

      data_df = data.frame()
      i = 1

      while (nrow(data_df) < nrow(data)) {
        selected_seed = sample(seeds_df$RESPONDENT_ID, 1, replace=TRUE)
        bootstrap_list[[i]] = data[data$SEED_ID == selected_seed, ]
        data_df = do.call(rbind, bootstrap_list)
        i = i + 1
      }

      return(data_df)
    }
    for (i in seq_len(n.times)){
      results_bootstrap_RDS[[i]] = boot_chain2(data)
    }
  }

  else if(type == 'boot_tree_one'){
    ######################################## function for tree boostrap 1
    tree.boot.one = function(data){
      all_seeds = data[which(data$SEED == 1),] # subset seeds only
      seed_ids <- sample(all_seeds$RESPONDENT_ID, size = nrow(all_seeds), replace = TRUE) # sample seeds only

      results = list()
      ### for each seed lets find their recruits
      for (i in 1:length(seed_ids)){
        seed_recruits = data[data$RECRUITER_ID %in% seed_ids[i],]
        r.recruits <- sample(seed_recruits$RESPONDENT_ID, size = nrow(seed_recruits), replace = TRUE)
        recruit_index_list <- match(r.recruits, seed_recruits$RESPONDENT_ID)
        results[[i]] <- seed_recruits[recruit_index_list, ]
      }

      results = do.call(rbind, results)


      full_recruitment_data <- results # initialize with the first wave of recruits
      current_wave_ids <- full_recruitment_data$RESPONDENT_ID # get unique respondent IDs from the first wave

      # Continue the process until there are no new recruits
      while (length(current_wave_ids) > 0) {
        recruits_df <- data[data$RECRUITER_ID %in% current_wave_ids, ] # Find recruits of the current wave

        if (nrow(recruits_df) > 0) { # Check if there are any new recruits
          r.recruits <- sample(recruits_df$RESPONDENT_ID, size = nrow(recruits_df), replace = TRUE)
          recruit_index_list <- match(r.recruits, recruits_df$RESPONDENT_ID)
          new_recruits <- recruits_df[recruit_index_list, ]

          if (nrow(new_recruits) > 0) {
            full_recruitment_data <- rbind(full_recruitment_data, new_recruits) # Append new recruits to the full data
            current_wave_ids <- new_recruits$RESPONDENT_ID # Update current_wave_ids with new recruits' IDs
          } else {
            current_wave_ids <- NULL # If no new recruits, break the loop
          }
        } else {
          current_wave_ids <- NULL # If no new recruits, break the loop
        }
      }

      seed_id = list()
      for (i in 1:length(seed_ids)){
        seed_id[[i]] = data[data$RESPONDENT_ID %in% seed_ids[i],]
      }

      seed_id = do.call(rbind, seed_id)
      results = rbind(seed_id, full_recruitment_data)
      return(results)
    }
    for (i in seq_len(n.times)){
      results_bootstrap_RDS[[i]] = tree.boot.one(data)
    }
  }



  else if (type == 'boot_tree_two'){
    ############################################ function for tree boostrap 2
    tree.boot.two = function(data){

      boot.results <- data.frame()

      while (nrow(boot.results)<nrow(data)){
        all_seeds = data[which(data$SEED == 1),] # subset seeds only
        seed_ids <- sample(all_seeds$RESPONDENT_ID, size = 1, replace = TRUE) # sample seeds only

        results = list()
        ### for each seed lets find their recruits
        for (i in 1:length(seed_ids)){
          seed_recruits = data[data$RECRUITER_ID %in% seed_ids[i],]
          r.recruits <- sample(seed_recruits$RESPONDENT_ID, size = nrow(seed_recruits), replace = TRUE)
          recruit_index_list <- match(r.recruits, seed_recruits$RESPONDENT_ID)
          results[[i]] <- seed_recruits[recruit_index_list, ]
        }

        results = do.call(rbind, results)


        full_recruitment_data <- results # initialize with the first wave of recruits
        current_wave_ids <- full_recruitment_data$RESPONDENT_ID # get unique respondent IDs from the first wave

        # Continue the process until there are no new recruits
        while (length(current_wave_ids) > 0) {
          recruits_df <- data[data$RECRUITER_ID %in% current_wave_ids, ] # Find recruits of the current wave

          if (nrow(recruits_df) > 0) { # Check if there are any new recruits
            r.recruits <- sample(recruits_df$RESPONDENT_ID, size = nrow(recruits_df), replace = TRUE)
            recruit_index_list <- match(r.recruits, recruits_df$RESPONDENT_ID)
            new_recruits <- recruits_df[recruit_index_list, ]

            if (nrow(new_recruits) > 0) {
              full_recruitment_data <- rbind(full_recruitment_data, new_recruits) # Append new recruits to the full data
              current_wave_ids <- new_recruits$RESPONDENT_ID # Update current_wave_ids with new recruits' IDs
            } else {
              current_wave_ids <- NULL # If no new recruits, break the loop
            }
          } else {
            current_wave_ids <- NULL # If no new recruits, break the loop
          }
        }
        seed_ids = data[data$RESPONDENT_ID == seed_ids,]
        boot.results = rbind(seed_ids, full_recruitment_data)
        boot.results = boot.results[complete.cases(boot.results$RESPONDENT_ID),]
      }
      return(boot.results)
    }
    for (i in seq_len(n.times)){
      results_bootstrap_RDS[[i]] = tree.boot.two(data)
    }
  }

  else if (type == 'boot_bidirectional_one'){
    bootstrap_bidirectional_one = function(data) {

      starting_ids = sample(data$RESPONDENT_ID, 10, replace = TRUE) # sample 10 IDs
      chains = vector("list", length(starting_ids))
      names(chains) = starting_ids

      for (i in seq_along(starting_ids)) {
        current_ids = starting_ids[i]     # tart with a selected ID and initialize the chain with that starting ID
        chain = c(current_ids)

        # Continue building the chain until no further connections can be added
        while (length(current_ids) > 0) {
          forward_nodes = unlist(data[data$RECRUITER_ID %in% current_ids, "RESPONDENT_ID", drop = FALSE]) # forward nodes
          backward_nodes = unlist(data[data$RESPONDENT_ID %in% current_ids, "RECRUITER_ID", drop = FALSE]) # backward nodes
          connected_nodes = c(forward_nodes, backward_nodes) # conencted nodes

          next_nodes = setdiff(connected_nodes, chain) # exclude already visited
          next_nodes = next_nodes[!is.na(next_nodes)] # sometimes used NA, exclude them
          if (length(next_nodes) == 0) break #       # If there are no new nodes to add, break the loop
          sampled_nodes = sample(next_nodes, size = length(next_nodes), replace = TRUE) # Sample with replacement from next_nodes
          chain = c(chain, sampled_nodes)  # Append the sampled nodes to the chain
          current_ids = unique(sampled_nodes) # Update current_ids for the next iteration, making sure to avoid duplications
        }
        chains[[i]] = chain     # Store the completed chain for the current starting ID
      }
      # now we need to transform chains into a dataframe, for this we create an empty data frame and the iterate over
      # each element in a list
      chains_df <- data.frame(chain_id = integer(), RESPONDENT_ID = integer(), stringsAsFactors = FALSE) # empty data frame
      for (i in seq_along(chains)) {
        single_chain_df <- data.frame(chain_id = rep(names(chains)[i],
                                                     length(chains[[i]])), RESPONDENT_ID = chains[[i]], stringsAsFactors = FALSE)
        chains_df <- rbind(chains_df, single_chain_df)
      }
      results = merge(chains_df, data, by='RESPONDENT_ID')
      return(results)
    }
    for (i in seq_len(n.times)){
      results_bootstrap_RDS[[i]] = bootstrap_bidirectional_one(data)
    }
  }

  else if (type == 'boot_bidirectional_two'){bootstrap_bidirectional_two = function(data) {

    # Initialize the data frame outside the while loop
    chains_df <- data.frame(chain_id = integer(), RESPONDENT_ID = integer(), stringsAsFactors = FALSE)

    # Loop until the number of rows in chains_df is greater than or equal to the number of rows in data
    while (nrow(chains_df) < nrow(data)) {
      starting_ids = sample(data$RESPONDENT_ID, 1, replace = TRUE) # sample IDs
      chains = vector("list", length(starting_ids))
      names(chains) = starting_ids

      for (i in seq_along(starting_ids)) {
        current_ids = starting_ids[i] # start with a selected ID and initialize the chain with that starting ID
        chain = c(current_ids)

        # Continue building the chain until no further connections can be added
        while (length(current_ids) > 0) {
          forward_nodes = unlist(data[data$RECRUITER_ID %in% current_ids, "RESPONDENT_ID", drop = FALSE])
          backward_nodes = unlist(data[data$RESPONDENT_ID %in% current_ids, "RECRUITER_ID", drop = FALSE])
          connected_nodes = c(forward_nodes, backward_nodes)

          next_nodes = setdiff(connected_nodes, chain)
          next_nodes = next_nodes[!is.na(next_nodes)]
          if (length(next_nodes) == 0) break
          sampled_nodes = sample(next_nodes, size = length(next_nodes), replace = TRUE)
          chain = c(chain, sampled_nodes)
          current_ids = unique(sampled_nodes)
        }
        chains[[i]] = chain
      }

      # Transform chains into a dataframe
      for (i in seq_along(chains)) {
        single_chain_df <- data.frame(chain_id = rep(names(chains)[i], length(chains[[i]])),
                                      RESPONDENT_ID = chains[[i]],
                                      stringsAsFactors = FALSE)
        chains_df <- rbind(chains_df, single_chain_df)
      }
    }
    results = merge(chains_df, data, by='RESPONDENT_ID')
    return(results)
  }
  for (i in seq_len(n.times)){
    results_bootstrap_RDS[[i]] = bootstrap_bidirectional_two(data)
  }
  }


  else{ warning('Specifcy one of six boostrap types')}

  results_bootstrap_RDS = lapply(seq(results_bootstrap_RDS),
                                 function(x) "[[<-"(results_bootstrap_RDS[[x]],
                                                    paste0("boot_n"), value = x))


  if (!return_data) {

    results_bootstrap_RDS = lapply(results_bootstrap_RDS,  # returns ID variables
                                   function(df) {data.frame(RESPONDENT_ID = df$RESPONDENT_ID, boot_n = df$boot_n)})

  }  else {
    results_bootstrap_RDS = do.call(rbind, results_bootstrap_RDS)
  }

  return(results_bootstrap_RDS)
}

