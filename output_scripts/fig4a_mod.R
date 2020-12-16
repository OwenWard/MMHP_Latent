sim_data_path <- here("output", "sims_m3_update","/")

#### need to edit here to account for separate files now

m3_files <- list.files(sim_data_path, pattern = "sim_model3_fit*")
#m1_files

# load(paste0(m1_path,m1_files[1]))
sim_model3_stan_sim1_total <- list()
sim_model3_stan_sim2_total <- list()
sim_model3_stan_sim3_total <- list()

for(i in seq_along(m3_files)) {
  load(paste0(sim_data_path,m3_files[i]))
  sim_model3_stan_sim1_total[[i]] <- sim_model3_stan_sim1
  sim_model3_stan_sim2_total[[i]] <- sim_model3_stan_sim2
  sim_model3_stan_sim3_total[[i]] <- sim_model3_stan_sim3
  
}


n_sim <- 50


par(mfrow=c(3, length(object_par$f_vec_1)),
    mar=c(0,0,1,0),
    oma=c(6.5,11,3.5,1))

y.ub <- c(12,12,12,12,12)
text.y <- c(rep(8,4),8)
text.x <- c(0.1,0.2,0.4,0.7,0.9)
text.x.adj <- c(0.3,0.4,0.6,0.5,0.8)

model_colors <- col_df %>% filter(method %in% c("C-HP",
                                                "C-DCHP",
                                                "C-MMHP"))

ind_colors <- model_colors %>% pull(cols)

for(m in c(1:3)){
  for(j in c(1:length(object_par$f_vec_1))){
    plot(1,4,xlim=c(0,1),ylim=c(0,y.ub[j]),
         type="n",xlab="",ylab="",xaxt="n",yaxt="n",axes=FALSE)
    
    ## Posterior density
    # current_result <- eval(parse(text=paste("sim_model3_stan_sim",m,sep="")))
    current_result <- eval(parse(text=paste("sim_model3_stan_sim",
                                            m,
                                            "_total",
                                            sep="")))
    for(i in c(1:n_sim)){
      lines(density(current_result[[i]]$f[,j]),
            col=add.alpha(model_colors[m], alpha=0.8),cex=0.6,lwd=0.8)
    }
    
    ## plot true value
    abline(v=object_par$f_vec_1[j],lwd=1.8,lty=2,col="gray30")
    
    ## plot box
    box.wdt <- 1.5
    box.col <- "black"
    box(lty = 'solid',col=box.col,lwd=0.7)
    
    if(m==3){
      text(text.x.adj[j],text.y[j],text.x[j],cex=2.4)
    }
  }
  ## bottom
  mtext(text=expression(f[1],f[2],f[3],f[4],f[5]), side=1, line=2.4, outer=TRUE, 
        at=c(1:length(object_par$f_vec_1))/length(object_par$f_vec_1)-0.1,
        cex=2)
  mtext(text="Latent rank variables", side=1, line=5.3, outer=TRUE, cex=2.25,
        at=0.5)
  ## top
  mtext(text="(a)", side=3, line=0, outer=TRUE, cex=2.5,
        at=0.5, font=2)
  ## left
  mtext(text=rev(c("C-HP  ","C-DCHP","C-MMHP")), side=2, line=0.5, outer=TRUE, 
        at=c(0:2)/3+0.15,
        cex=1.75, las=2)
}


### histogram of posterior means instead

m1 <- sim_model3_stan_sim1_total %>%
  map( ~apply(.x$f, 2, mean)  ) %>%
  enframe() %>%
  unnest(cols = value) %>%
  mutate(node_id = rep(c("f1","f2","f3","f4","f5"), 50),
         method = "C-HP")

m2 <- sim_model3_stan_sim2_total %>%
  map( ~apply(.x$f, 2, mean)  ) %>%
  enframe() %>%
  unnest(cols = value) %>%
  mutate(node_id = rep(c("f1","f2","f3","f4","f5"), 50),
         method = "C-DCHP")

m3 <- sim_model3_stan_sim3_total %>%
  map( ~apply(.x$f, 2, mean)  ) %>%
  enframe() %>%
  unnest(cols = value) %>%
  mutate(node_id = rep(c("f1","f2","f3","f4","f5"), 50),
         method = "C-MMHP")

rank_means <- bind_rows(m1,m2,m3)

rank_means

true_values <- tibble(
  id = c("f1", "f2","f3","f4","f5"),
  values = c(0.1,0.2,0.4,0.7,0.9),
  val_names = c("   0.1","0.2","0.4","0.7","0.9")
)


rank_means %>% 
  left_join(true_values, by = c("node_id" = "id")) %>%
  ggplot(aes(value, fill = method)) + 
  geom_histogram() + 
  geom_vline(aes(xintercept = values), alpha = 0.6,
             linetype = "dashed") +
  ### add text here
  geom_text(aes(values, 16, label = val_names),
            size = 4, hjust = 1) +
  scale_fill_manual(values = ind_colors) +
  facet_grid(cols = vars(node_id), rows = vars(method)) + 
  theme_bw() + theme(
    plot.title = element_text(size = 22, hjust = 0.5, face="bold"),
    text = element_text(size = 12),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "", y = "") 

### now need to add the true values

rank_means %>% 
  left_join(true_values, by = c("node_id" = "id")) %>%
  ggplot(aes(method, value, fill = method)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = values), alpha = 0.6,
             linetype = "dashed") +
  ### add text here
  geom_text(aes(2, values, label = values), vjust = -0.75,
            hjust = -0.1,
            size = 3) +
  scale_fill_manual(values = ind_colors, name = "Method") +
  facet_grid(cols = vars(node_id), space = "free") + 
  theme_bw() + theme(
    plot.title = element_text(size = 22, hjust = 0.5, face="bold"),
    text = element_text(size = 12),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(x = "", y = "", title = "(a)")
  NULL
