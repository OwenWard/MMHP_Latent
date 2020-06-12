expertRankHierarchy <- function(raw_data){
  win_lost_df <- raw_data[-which(raw_data[,"Actor"]=="Start"|raw_data[,"Actor"]=="End"),
                          c("Actor","Recipient")]
  win_lost_matrix <- get_wl_matrix(win_lost_df)
  isi.out <- isi98(win_lost_matrix)
  expert_rank <- as.numeric(isi.out$best_order)
  return(expert_rank)
}