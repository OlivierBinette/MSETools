all: real_data_analysis dga_sensitivity sparsemse_sensitivity internal_consistency empirical_trajectories theory
.PHONY: real_data_analysis dga_sensitivity sparsemse_sensitivity internal_consistency empirical_trajectories theory clean

real_data_analysis: real_data/real_data_analysis.Rmd
	Rscript -e "rmarkdown::render('real_data/real_data_analysis.Rmd', run_pandoc=FALSE)"

dga_sensitivity: sensitivity_analysis/dga/dga_sensitivity.Rmd
	Rscript -e "rmarkdown::render('sensitivity_analysis/dga/dga_sensitivity.Rmd', run_pandoc=FALSE)"

sparsemse_sensitivity: sensitivity_analysis/sparsemse/sparsemse_sensitivity.Rmd
	Rscript -e "rmarkdown::render('sensitivity_analysis/sparsemse/sparsemse_sensitivity.Rmd', run_pandoc=FALSE)"

internal_consistency: internal_consistency/internal_consistency.Rmd
	Rscript -e "rmarkdown::render('internal_consistency/internal_consistency.Rmd', run_pandoc=FALSE)"

empirical_trajectories: empirical_trajectories/empirical_trajectories.Rmd
	Rscript -e "rmarkdown::render('empirical_trajectories/empirical_trajectories.Rmd', run_pandoc=FALSE)"
	
theory: theory/beta_model.Rmd
	Rscript -e "rmarkdown::render('theory/beta_model.Rmd', run_pandoc=FALSE)"

clear:
	find . -name "*.md" -o -name "*_files" | xargs rm -rf
	
clear_cache:
	find . -name "*_cache" | xargs rm -rf

