
#' @import slurmR
#' @export
batch.estimates <- function(fits, njobs, mc.cores = 1L, job_name = "batch_MSE", tmp_path = ".", plan = "submit", partition = NULL, ...) {
  if (!is.null(partition)) opts_slurmR$set_opts(partition = partition)

  Slurm_lapply(fits, function(fit) {estimates(fit, ...)}, njobs = njobs, mc.cores = mc.cores, job_name = job_name, tmp_path = tmp_path, plan = plan)
}

#' @import dplyr tibble tidyr
#' @export
grid.estimates <- function(models, datasets, slurm=FALSE, plan="collect", tidy=TRUE, ...) {
  params = expand.grid(model=models, dataset=datasets, stringsAsFactors=FALSE)
  model_fits = apply(params, 1, function(x) x$model(x$dataset))
  if (slurm) {
    ests = batch.estimates(model_fits, plan=plan, ...)
  } else {
    ests = estimates(model_fits)
  }

  if (tidy) {
    return(params %>%
             as_tibble %>%
             mutate_all(names) %>%
             add_column(ests) %>%
             unnest_wider(ests))
  } else {
    return(ests)
  }
}
