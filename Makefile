SRCDIR = src
DATADIR = data
OUTDIR = out
BINDIR = /nfs/research2/beltrao/software-rh7/bin

# Parameters
TABLE_STRATEGIES = max-rows max-cols max-rows-max-cols
TABLE_STRATEGY ?= max-rows-max-cols
ASSOC_METHODS = pcor pcor-filter scor scor-filter
ASSOC_METHOD ?= scor
ASSOCNET_FILTER_METHOD ?= deconvolution
ASSOCNET_FILTER_SCALE_METHOD ?= standard
DECONVOLUTION_A ?= 1.0
DECONVOLUTION_B ?= 0.99

# External data sets
# Kinase activity predictions
GSEA_DATA = $(DATADIR)/log.wGSEA.kinase_condition.clean.Rdata
# Raw phosphosite data
PSITE_DATA = $(DATADIR)/esetNR.Rdata
# specifies kinases perturbed in each condition
KIN_COND_PAIRS = $(DATADIR)/kinase-condition-pairs.tsv
# tries to calculate which kinases are perturbed in each condition
KIN_INVIVO_CONDS = $(DATADIR)/kinase_invivoconditions.tsv 

# Generated data sets Initial kinase-activity tables, trying to
# maximize #rows, #cols or both.  Raw and imputed datasets
KINACT_DATA = $(DATADIR)/kinase-activity-$(TABLE_STRATEGY).tsv
IMP_KINACT_DATA = $(DATADIR)/kinase-activity-$(TABLE_STRATEGY)-imp.tsv
# EGF-signaling-related kinase activity
EGF_KIN_ACT_DATA = $(foreach T,perturbed unperturbed,DATADIR)/egf-kinase-activity-$(T).tsv)
IMP_EGF_KIN_ACT_DATA = $(foreach T,perturbed unperturbed,$(DATADIR)/egf-kinase-activity-$(T)-imp.tsv)
# Association table output
KINACT_ASSOC = $(OUTDIR)/kinase-activity-$(TABLE_STRATEGY)-$(ASSOC_METHOD).tsv
# Posterior probabilities
KINACT_POSTERIOR_PROB = $(OUTDIR)/kinase-activity-$(TABLE_STRATEGY)-$(ASSOC_METHOD)-posterior.tsv

# Program Options
ASSOCNET_PARAMS = --unbiased-correlation --p-method=none
ASSOCNET_FILTER_PARAMS = --method=$(ASSOCNET_FILTER_METHOD) --scale-method=$(ASSOCNET_FILTER_SCALE_METHOD) --header-in
ifeq ($(NET_FILTER_METHOD),deconvolution)
ASSOCNET_FILTER_PARAMS += --deconvolution-a=$(DECONVOLUTION_A)	\
						  --deconvolution-b=$(DECONVOLUTION_B)
endif

# Programs
PYTHON3 ?= $(BINDIR)/python3
PYTHON ?= $(PYTHON3)
RSCRIPT ?= $(BINDIR)/Rscript
ASSOCNET ?= $(BINDIR)/assocnet $(ASSOCNET_PARAMS)
ASSOCNET_FILTER ?= $(BINDIR)/assocnet-filter $(NET_FILTER_PARAMS)

# Scripts
GEN_KINACT_TBL_SCRIPT = $(SRCDIR)/gen-activity-table.r
POSTERIOR_PROB_SCRIPT = $(SRCDIR)/posterior-prob.r

# Phony targets
.PHONY: posterior
posterior: $(KINACT_POSTERIOR_PROB)

.PHONY: cor-results
assoc-results: $(KINACT_ASSOC)

.PHONY: data
data: $(KINACT_DATA) $(IMP_KINACT_DATA) $(EGF_KIN_ACT_DATA) $(IMP_EGF_KIN_ACT_DATA)

.PHONY: clean-data
clean-data:
	rm -v $(KINACT_DATA) $(IMP_KINACT_DATA) $(EGF_KIN_ACT_DATA) $(IMP_EGF_KIN_ACT_DATA)

.PHONY: clean-ossoc
clean-assoc:
	rm -v $(KINACT_ASSOC)

.PHONY: clean-posterior
clean-posterior:
	rm -v $(KINACT_POSTERIOR_PROB)

.PHONY: clean-results
clean-results: clean-assoc clean-posterior

.PHONY: clean
clean: clean-data clean-results

# Rules
$(OUTDIR):
	mkdir -p $(OUTDIR)

$(KINACT_DATA) \
$(IMP_KINACT_DATA) \
$(EGF_KIN_ACT_DATA) \
$(IMP_EGF_KIN_ACT_DATA): $(GEN_KINACT_TBL_SCRIPT) $(GSEA_DATA) $(KIN_COND_PAIRS)
	$(PYTHON) $(GEN_KINACT_TBL_SCRIPT)

$(OUTDIR)/%-pcor.tsv: $(DATADIR)/%-imp.tsv $(OUTDIR)
	$(ASSOCNET) --method=pearson $< >$@

$(OUTDIR)/%-scor.tsv: $(DATADIR)/%-imp.tsv $(OUTDIR)
	$(ASSOCNET) --method=spearman $< >$@

$(OUTDIR)/%-filter.tsv: $(OUTDIR)/%.tsv
	$(ASSOCNET_FILTER) --header-in $< >$@

$(OUTDIR)/%-posterior.tsv: $(OUTDIR)/%.tsv
	$(RSCRIPT) $(POSTERIOR_PROB_SCRIPT) $(ASSOC_METHOD) $< $@
