


#' Data processing in an RDS study
#'
#' @param data Data frame; Should contain ID numbers for the nodes in the social network, corresponding redeemed coupon number, and issued coupon number
#' @param unique_id Character vector; The column name of the column in the data that represents the ID numbers for the nodes in the social network
#' @param redeemed_coupon Character vector; The column name of the column in the data that represents coupon numbers of coupons redeemed by respondents when participating in the survey
#' @param issued_coupons Character vector; The column name of the column in the data that represents the coupon numbers of coupons issued to respondents
#' @param degree Numeric vector;  The column name of the column in the data that represents the degree of the respondents
#' @param zero_degree_method
#' @param NA_degree_method Character; This parameter is used to set the method for imputing missing values in the 'degree' variable. There are three methods to choose from: mean, median, and hotdeck. If this parameter is not set, the default method is hotdeck.
#' @param result Character; This parameter is used to set the method for imputing zero values in the 'degree' variable. Three methods are available for selection: mean, median, and hotdeck. If this parameter is not set, the default imputation method is hotdeck
#'
#' @return
#' \item{ID}{Character vector; Renamed unique_id variable}
#' \item{R_CP}{Character vector; Renamed redeemed coupon variable}
#' \item{T_CP1 - T_CPn}{Character vector(s); Renamed issued coupon variable}
#' \item{DEGREE}{Numeric vector; Renamed degree variable}
#' \item{WAVE}{Numeric vector; Indicates which round the node was introduced into the survey, and the value of Seed is 0}
#' \item{S_ID}{chacarter vector; Indicates the ID of the seed corresponding to the node. The value of the seed is itself}
#' \item{R_ID}{Character vector; Indicates the ID of the node who recruits the node joinning the survey. The value of seed is NA}
#' \item{SEED}{Numeric vector: Values are only 0 and 1, they are used to indicate whether the node is seed or not. If it is seed, the value is 1, if not, it is 0.}
#' \item{CT_T_CT}{Numeric data, indicates how many coupons have been issued to this node to invite others to join the survey}
#' \item{CT_T_CP_USED}{Numeric data, indicates how many coupons issued to this node have been used to invite others to join the survey.}
#' \item{DEGREE_MODIFY}{Numeric data, indicates the `degree` data after imputation of zero values and missing values.}
#'
#' @examples
#' # for preprocessing use RDStoydata
#'
#' data('RDStoydata')
#'
#' rds_data <- RDSdata(data = RDStoydata,unique_id = "ID",
#' redeemed_coupon = "CouponR",
#' issued_coupon = c("Coupon1",
#'                  "Coupon2",
#'                  "Coupon3"),
#'                degree = "Degree",
#'                result = c('Age','Sex'))
#' @export

RDSdata <- function(
    data,
    unique_id,
    redeemed_coupon,
    issued_coupons,
    degree,
    zero_degree_method = 'hotdeck',
    NA_degree_method = 'hotdeck',
    result
) {
  Warning_Function_ID_missing_data<- function(unique_id,data)if (any(is.na(data[,1]))) {
    stop("Function operation has been interrupted.
         Please make sure there are no missing values in ",
         unique_id," before trying again.")}

  Warning_Function_ID_missing_data(unique_id,data)

  if(unique_id == redeemed_coupon) {
    unique_id <- "ID_NEW"
    df[[unique_id]] <- df[[redeemed_coupon]]
  }

  # select relevant columns
  df<-subset(data, select = c(unique_id, redeemed_coupon, issued_coupons,degree,result))
  df <- data.frame(lapply(df, as.character))
  df[df == ''] <- NA

  # rename columns
  issued_coupons <- paste0("T_CP", 1:length(issued_coupons))
  names(df) <- c("ID", "R_CP", issued_coupons,"DEGREE",result)

  get_recruiter_ID <- function(df) {
    id <- df$ID
    rc <- replace(df$R_CP, is.na(df$R_CP), '')
    tcs<-subset(df, select = -c(ID, R_CP))

    rid <- apply(tcs, 2, function(tc) id[match(rc, tc)])
    rid[is.na(rid)] <- ''
    rid <- apply(rid, 1, paste, collapse = '')
    rid[rid == ''] <- NA
    return(rid)
  }

  # construct recruitment chains
  holder <- df
  df_recruiter <- data.frame(R_ID0 = df$ID)
  i <- 0

  while(!all(is.na(df_recruiter[[paste0("R_ID", i)]]))) {
    rid <- get_recruiter_ID(holder)
    holder <- data.frame(ID = rid)
    holder$RowOrder <- seq_along(holder$ID)
    holder <- merge(holder, df, by = "ID", all.x = TRUE)
    holder <- holder[order(holder$RowOrder),]
    holder$RowOrder <- NULL

    i <- i + 1
    df_recruiter[[paste0("R_ID", i)]] <- rid
  }

  # get wave number
  df$WAVE <- apply(df_recruiter, 2, function(c) as.numeric(!is.na(c)))
  df$WAVE <- rowSums(df$WAVE)
  df$WAVE <- df$WAVE - 1

  # get seed ID
  df$S_ID <- mapply(function(i, j) df_recruiter[i, j+1],
                    1:nrow(df_recruiter),
                    df$WAVE)

  # get recruiter ID
  df$R_ID <- df_recruiter[["R_ID1"]]

  # classify seed
  df$SEED <- ifelse(is.na(df$R_ID), 1, 0)

  # find the number of coupons issued
  df$CT_T_CP <- rowSums(apply(df[issued_coupons], 2, function(c) as.numeric(!is.na(c))))

  # find the number of coupons used
  cu <- aggregate(cbind(CT_T_CP_USED = df$R_ID) ~ R_ID, data = df, FUN = length)
  names(cu)[names(cu) == "R_ID"] <- "ID"
  cu <- na.omit(cu)

  df$RowOrder <- 1:nrow(df)
  df <- merge(df, cu, by = "ID", all.x = TRUE)
  df <- df[order(df$RowOrder), ]
  df$RowOrder <- NULL
  df$CT_T_CP_USED <- replace(df$CT_T_CP_USED, is.na(df$CT_T_CP_USED), 0)

  #degree imputation
  #For calculation trans to numeric
  df$DEGREE<-as.numeric(df$DEGREE)
  #Select the rows to provide imputation data(non zero$ non NA)
  selected_rows <- df[df$DEGREE != 0 & !is.na(df$DEGREE), ]

  #Select the rows that need to be imputation
  zero_row<-df[df$DEGREE == 0&!is.na(df$DEGREE),]
  NA_row<-df[is.na(df$DEGREE),]
  #imputation for NA
  if(nrow(NA_row)!=0){
    if(NA_degree_method=='hotdeck'){
      for (i in 1:nrow(NA_row)) {
        # Randomly select a number from the row providing the imputed data
        NA_row$DEGREE[i]<- sample(selected_rows$DEGREE, 1)
      }
    }else if(NA_degree_method=='mean'){
      NA_row$DEGREE<-mean(selected_rows$DEGREE)
    }else if(NA_degree_method=='median'){
      NA_row$DEGREE<-median(selected_rows$DEGREE)
    }else{
      stop("Invalid NA_degree_method value.")
    }
  }
  #imputation for 0
  if(nrow(zero_row)!=0){
    if(zero_degree_method=='hotdeck'){
      for (i in 1:nrow(zero_row)) {
        # Randomly select a number from the row providing the imputed data
        zero_row$DEGREE[i]<- sample(selected_rows$DEGREE, 1)
      }
    }else if(zero_degree_method=='mean'){
      zero_row$DEGREE<-mean(selected_rows$DEGREE)
    }else if(zero_degree_method=='median'){
      zero_row$DEGREE<-median(selected_rows$DEGREE)
    }else{
      stop("Invalid zero_degree_method value.")
    }
  }

  #combine modified data
  #combine the imputation result fo 0 and NA
  df_sub<-rbind(zero_row,NA_row)

  #use ID to match them
  ids_to_update <- df_sub$ID
  for(ID in ids_to_update) {
    row_to_update <- which(df$ID == ID)
    df$DEGREE[row_to_update] <- df_sub$DEGREE[df_sub$ID == ID]
  }
  #calculate the WEIGHT
  df$WEIGHT<-1/df$DEGREE
  #For storage the data trans back to character
  df$DEGREE<-as.character(df$DEGREE)
  df$WEIGHT<-as.character(df$WEIGHT)
  return(df)
}




