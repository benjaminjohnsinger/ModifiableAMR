SHELL := /bin/bash

SCENARIO ?= main
ANALYSIS_SCENARIO ?= $(SCENARIO)
R := Rscript

.PHONY: help \
  data models burden figures tables validate \
  validate-inputs validate-metrics validate-stage-data validate-all validate-supp \
  promote-metrics \
  data-smoke models-smoke burden-smoke figures-smoke tables-smoke smoke \
  clean-derived manuscript \
  models-hic models-lmic models-nagorsen models-iqvia \
  burden-region-drug burden-region-pathogen \
  permutation-models wilcoxon-test \
  tables-supp supplementary all-scenarios

help:
	@echo ""
	@echo "AMR reproducibility workflow — $(shell date +%Y)"
	@echo ""
	@echo "  MAIN PAPER (Figures 1-4):"
	@echo "    make manuscript           Full pipeline for main-paper outputs (default SCENARIO=main)"
	@echo "    make data                 Prepare and merge data"
	@echo "    make models               Fit regression models"
	@echo "    make burden               Estimate avertable burden (bootstrapped)"
	@echo "    make figures              Generate Figures 1-4"
	@echo "    make validate-all         Run all validation checks"
	@echo ""
	@echo "  SUPPLEMENTARY ANALYSES:"
	@echo "    make supplementary        Run all supplementary analyses end-to-end"
	@echo "    make models-hic           Fit models on HIC countries only"
	@echo "    make models-lmic          Fit models on LMIC countries only"
	@echo "    make models-nagorsen      Fit models for Nagorsen/hospital scenario"
	@echo "    make models-iqvia         Fit models on raw IQVIA data"
	@echo "    make burden-region-drug   Burden aggregated by region x drug"
	@echo "    make burden-region-pathogen  Burden aggregated by region x pathogen"
	@echo "    make permutation-models   Run permuted bootstrap models (all 7 classes; ~7x model runtime)"
	@echo "    make wilcoxon-test        Run Wilcoxon test (Table 2 prerequisite)"
	@echo "    make tables-supp          Generate supplementary Tables 2 and 3 (+ validate)"
	@echo "    make validate-supp        Validate supplementary table outputs"
	@echo "    make promote-metrics      Promote passing pending metrics to canonical"
	@echo ""
	@echo "  SMOKE TESTS:"
	@echo "    make smoke                Full smoke test (all stages, main scenario)"
	@echo "    make models-smoke         Models stage smoke test"
	@echo "    make burden-smoke         Burden stage smoke test"
	@echo "    make figures-smoke        Figures stage smoke test"
	@echo "    make tables-smoke         Tables stage smoke test"
	@echo ""
	@echo "  OTHER:"
	@echo "    make all-scenarios        Run manuscript workflow for every scenario"
	@echo "    make clean-derived        Remove derived/canonical output files"
	@echo ""
	@echo "  SCENARIO variable (default: main):"
	@echo "    make <target> SCENARIO=hic|lmic|main|..."
	@echo ""

data:
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/run_prepare_data.R

models:
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/run_fit_models.R

burden:
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/run_estimate_burden.R

figures:
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/run_generate_figures.R

tables:
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/run_tables.R

validate:
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/validate_outputs.R $(ANALYSIS_SCENARIO)

validate-inputs:
	$(R) scripts/validate_data_manifest.R

validate-metrics:
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/validate_metrics.R

validate-stage-data:
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/validate_stage_data.R

validate-all: validate-inputs validate-stage-data validate validate-metrics

promote-metrics:
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/promote_metrics.R $(ANALYSIS_SCENARIO)

validate-supp:
	$(R) scripts/validate_outputs.R permutation scripts/run_tables.R
	$(R) scripts/validate_outputs.R supplementary scripts/run_tables.R
	$(R) scripts/validate_metrics.R permutation

data-smoke:
	$(MAKE) data ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO)
	$(MAKE) validate-stage-data ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO)

models-smoke:
	AMR_DEV_SMOKE=1 ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/run_fit_models.R
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/validate_outputs.R $(ANALYSIS_SCENARIO) scripts/run_fit_models.R

burden-smoke:
	AMR_DEV_SMOKE=1 ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/run_estimate_burden.R
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/validate_outputs.R $(ANALYSIS_SCENARIO) scripts/run_estimate_burden.R

figures-smoke:
	$(MAKE) figures ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO)
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/validate_outputs.R $(ANALYSIS_SCENARIO) scripts/run_generate_figures.R

tables-smoke:
	$(MAKE) tables ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO)
	ANALYSIS_SCENARIO=$(ANALYSIS_SCENARIO) $(R) scripts/validate_outputs.R $(ANALYSIS_SCENARIO) scripts/run_tables.R

smoke: data-smoke models-smoke burden-smoke figures-smoke tables-smoke

clean-derived:
	@echo "Removing declared derived intermediates and canonical artifacts..."
	@awk -F, 'NR>1 && $$2=="intermediate" {print $$3}' manifests/data_manifest.csv | while read -r p; do \
	  [ -n "$$p" ] && [ -e "$$p" ] && rm -f "$$p"; \
	done
	@awk -F, 'NR>1 && $$8=="canonical" {print $$6}' manifests/figure_manifest.csv | while read -r p; do \
	  [ -n "$$p" ] && [ -e "$$p" ] && rm -f "$$p"; \
	done
	@rm -f Figure1_lagged.meta.txt Figure2.meta.txt Figure3.meta.txt

manuscript: data models burden figures validate-all

# ---------------------------------------------------------------------------
# Supplementary model runs
# ---------------------------------------------------------------------------

models-hic:
	ANALYSIS_SCENARIO=hic $(R) scripts/run_fit_models.R

models-lmic:
	ANALYSIS_SCENARIO=lmic $(R) scripts/run_fit_models.R

models-nagorsen:
	ANALYSIS_SCENARIO=hospital_nagorsen $(R) scripts/run_fit_models.R

models-iqvia:
	ANALYSIS_SCENARIO=raw_iqvia $(R) scripts/run_fit_models.R

# ---------------------------------------------------------------------------
# Supplementary burden runs
# ---------------------------------------------------------------------------

burden-region-drug:
	ANALYSIS_SCENARIO=burden_drug_region $(R) scripts/run_estimate_burden.R

burden-region-pathogen:
	ANALYSIS_SCENARIO=burden_pathogen_region $(R) scripts/run_estimate_burden.R

# ---------------------------------------------------------------------------
# Permutation models (for Table 2 Wilcoxon test)
# Runs one permutation model per antibiotic class (~7 x full model runtime).
# Set AMR_PERMUTATION_CLASS to one class, or use this target to run all seven.
# ---------------------------------------------------------------------------

permutation-models:
	@echo "[permutation] Running permutation models for all 7 antibiotic classes..."
	@for ab in J01A J01C J01D J01E J01F J01G J01M; do \
	  echo "[permutation] Class: $$ab"; \
	  ANALYSIS_SCENARIO=permutation AMR_PERMUTATION_CLASS=$$ab $(R) scripts/run_fit_models.R || exit 1; \
	done
	@echo "[permutation] All permutation models complete."

# ---------------------------------------------------------------------------
# Wilcoxon test (Table 2 prerequisite — requires permutation model outputs)
# ---------------------------------------------------------------------------

wilcoxon-test:
	ANALYSIS_SCENARIO=main $(R) scripts/run_wilcoxon_test.R

# ---------------------------------------------------------------------------
# Supplementary tables (Table 2 + Table 3)
# Table 2 requires wilcoxon-test to have run first.
# Table 3 requires the main models stage to have run first.
# ---------------------------------------------------------------------------

tables-supp: tables validate-supp

# ---------------------------------------------------------------------------
# Convenience target: all supplementary analyses end-to-end
# ---------------------------------------------------------------------------

supplementary: models-hic models-lmic models-nagorsen models-iqvia \
               burden-region-drug burden-region-pathogen \
               permutation-models wilcoxon-test tables-supp validate-supp
	@echo ""
	@echo "Supplementary analyses complete."
	@echo ""

all-scenarios:
	@for s in main hic lmic raw_iqvia hospital_nagorsen burden_optimistic burden_pessimistic burden_lower_region burden_upper_region burden_drug_region burden_pathogen_region; do \
	  echo "Running manuscript workflow for $$s"; \
	  $(MAKE) manuscript SCENARIO=$$s || exit 1; \
	done
