# Latent Space Hawkes Processes for Ranking Interaction Data


A repo describing the final code to run and analyse the entire latent space Hawkes process procedure for [Network Hawkes Process Models for Exploring Latent Hierarchy in Social Animal Interactions](https://arxiv.org/abs/2012.09598)



## Real Data Procedure

- Use `run_scripts` folder.
- Run `c_hp.R, c_dchp.R` and `c_mmhp_dc.R` for each of the 10 Mice cohorts
- Run `i_mmhp.R` for each of the 10 Mice cohorts
- Run `predict.R` for each of the 10 Mice cohorts
- Run `overall_predict.R` once, to summarise the predictions from all cohorts
- Run `diagnostics.R` for each of the 10 Mice cohorts



## Simulation Study Procedures

- Use `run_scripts` folder.
- Run `sim_fits.R` for each of the 50 simulation studies.
- Run `sim_rank_isi.R` to do the corresponding rank simulations

## Produce the results

- `Paper_Plots.Rmd` in `output_scripts` creates the plots which are used in the
paper.

