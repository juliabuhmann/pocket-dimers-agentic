# H5 mixed models (DESIGN.md v3 H5). First attempt, not independently verified.
# Usage: Rscript scripts/h5_glmm.R <protein_residues_X.parquet> <out.tsv> <outcome: plp|am> [catcol=cat] [shuffle_seed=-1] [nAGQ=1]
# plp: glmer(y ~ cat + z(rsa_iso) + z(plddt) + (1|accession), binomial)   reference category = interface
# am : lmer(am_mean ~ cat + z(rsa_iso) + z(plddt) + (1|accession))          (Wald z p-values)
# Buried residues (and the "excluded" level of the burrow-quintile sensitivity) are excluded. Conservation covariate: not included (no cheap score available; see docs/H5_RESULTS.md).
# shuffle_seed >= 0: category labels permuted within protein (negative control).
suppressPackageStartupMessages({library(lme4); library(arrow); library(data.table)})
a <- commandArgs(TRUE)
inf <- a[1]; outf <- a[2]; outcome <- a[3]
catcol <- if (length(a) >= 4) a[4] else "cat"
seed <- if (length(a) >= 5) as.integer(a[5]) else -1L
nagq <- if (length(a) >= 6) as.integer(a[6]) else 1L
d <- as.data.table(read_parquet(inf))
d[, cat := as.character(get(catcol))]
d <- d[!(cat %in% c("buried", "excluded")) & is.finite(rsa_iso) & is.finite(plddt)]
if (outcome == "am") d <- d[is.finite(am_mean)]
if (seed >= 0) { set.seed(seed); d[, cat := if (.N > 1) sample(cat) else cat, by = accession] }
d[, cat := relevel(factor(cat), ref = "interface")]
d[, rsa_z := as.numeric(scale(rsa_iso))]; d[, plddt_z := as.numeric(scale(plddt))]
t0 <- Sys.time()
if (outcome == "plp") {
  m <- glmer(y ~ cat + rsa_z + plddt_z + (1 | accession), data = d, family = binomial, nAGQ = nagq,
             control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
} else {
  m <- lmer(am_mean ~ cat + rsa_z + plddt_z + (1 | accession), data = d, REML = FALSE,
            control = lmerControl(calc.derivs = FALSE))
}
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cf <- summary(m)$coefficients
est <- cf[, 1]; se <- cf[, 2]; z <- est / se; p <- 2 * pnorm(-abs(z))
res <- data.table(term = rownames(cf), estimate = est, se = se, z = z, p = p,
                  lo = est - 1.959964 * se, hi = est + 1.959964 * se)
if (outcome == "plp") res[, `:=`(OR = exp(estimate), OR_lo = exp(lo), OR_hi = exp(hi))]
res[, `:=`(input = basename(inf), outcome = outcome, catcol = catcol, shuffle_seed = seed, nAGQ = nagq,
           n_res = nrow(d), n_prot = uniqueN(d$accession), n_events = if (outcome == "plp") sum(d$y) else NA_integer_,
           re_sd = attr(VarCorr(m)$accession, "stddev"), converged_msgs = paste(m@optinfo$conv$lme4$messages, collapse = ";"),
           fit_seconds = round(el))]
fwrite(res, outf, sep = "\t")
print(res[, .(term, estimate, se, p)])
