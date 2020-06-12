naiveRankHierarchy <- function(raw_data){
  win_lost_df <- raw_data[-which(raw_data[,"Actor"]=="Start"|raw_data[,"Actor"]=="End"),
                          c("Actor","Recipient")]
  win_lost_matrix <- get_wl_matrix(win_lost_df)
  naive_rank <- as.numeric(row.names(win_lost_matrix))[order(rowSums(win_lost_matrix), decreasing = TRUE)]
  return(naive_rank)
}