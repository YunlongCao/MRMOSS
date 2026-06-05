library(MRMOSS)

print(list_mrmoss_examples())

show_mrmoss_example("cvd_lipid_panel")

res <- run_mrmoss_example("cvd_lipid_panel", maxiter = 500)

print(res$input_qc)
print(res$global_lrt)
print(res$outcome_lrt)
print(res$domain_lrt)
