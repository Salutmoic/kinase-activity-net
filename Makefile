SRCDIR = src
DATADIR = data
OUTDIR = out
BINDIR = /nfs/research2/beltrao/software-rh7/bin

# Parameters
TABLE_STRATEGIES = max-rows max-cols max-rows-max-cols
TABLE_STRATEGY ?= max-rows-max-cols
ASSOC_METHODS = pcor pcor-filter scor scor-filter nfchisq
ASSOC_METHOD ?= scor
ASSOCNET_FILTER_METHOD ?= deconvolution
ASSOCNET_FILTER_SCALE_METHOD ?= standard
DISCR_METHOD ?= mclust.whole
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
# PhosphositePlus data
PHOSPHOSITE_PLUS_VERSION = 2017-09-08
FULL_KIN_SUBSTR_TABLE = $(DATADIR)/Kinase_Substrate_Dataset_$(PHOSPHOSITE_PLUS_VERSION)
FULL_REG_SITES_TABLE = $(DATADIR)/Regulatory_sites_$(PHOSPHOSITE_PLUS_VERSION)
FULL_PHOS_SITES_TABLE = $(DATADIR)/Phosphorylation_site_dataset_$(PHOSPHOSITE_PLUS_VERSION)

# Generated data sets Initial kinase-activity tables, trying to
# maximize #rows, #cols or both.  Raw and imputed datasets
KINACT_DATA = $(DATADIR)/kinase-activity-$(TABLE_STRATEGY).tsv
IMP_KINACT_DATA = $(DATADIR)/kinase-activity-$(TABLE_STRATEGY)-imp.tsv
DISCR_KINACT_DATA = $(DATADIR)/kinase-activity-$(TABLE_STRATEGY)-discr.tsv
# EGF-signaling-related kinase activity
EGF_KIN_ACT_DATA = $(foreach T,perturbed unperturbed,DATADIR)/egf-kinase-activity-$(T).tsv)
IMP_EGF_KIN_ACT_DATA = $(foreach T,perturbed unperturbed,$(DATADIR)/egf-kinase-activity-$(T)-imp.tsv)
# Association table output
KINACT_ASSOC = $(OUTDIR)/kinase-activity-$(TABLE_STRATEGY)-$(ASSOC_METHOD).tsv
# Posterior probabilities
KINACT_POSTERIOR_PROB = $(OUTDIR)/kinase-activity-$(TABLE_STRATEGY)-$(ASSOC_METHOD)-posterior.tsv
# Human amino acid frequencies
AA_FREQS = $(DATADIR)/aa-freqs.tsv
# PSSM files
KIN_KIN_SCORES = $(OUTDIR)/kinase_kinase_scores.tsv
KIN_KNOWN_PSITE_SCORES = $(OUTDIR)/known_kinase_psite-score
KIN_SCORE_DIST = $(OUTDIR)/kinase_distributions.tsv
# PhosphositePlus derived files
KIN_SUBSTR_TABLE = $(DATADIR)/reduced_kinase_table.tsv
HUMAN_KINASE_TABLE = $(DATADIR)/human_kinase_table.tsv
PHOSPHOSITES = $(DATADIR)/phosphosites_reduced.tsv

# Program Options
ASSOCNET_PARAMS = --unbiased-correlation --p-method=none
ASSOCNET_FILTER_PARAMS = --method=$(ASSOCNET_FILTER_METHOD) --scale-method=$(ASSOCNET_FILTER_SCALE_METHOD) --header-in
ifeq ($(NET_FILTER_METHOD),deconvolution)
ASSOCNET_FILTER_PARAMS += --deconvolution-a=$(DECONVOLUTION_A)	\
						  --deconvolution-b=$(DECONVOLUTION_B)
endif

# Programs
PYTHON3 ?= $(BINDIR)/python3
PYTHON2 ?= $(BINDIR)/python2
PYTHON ?= $(PYTHON3)
RSCRIPT ?= $(BINDIR)/Rscript
ASSOCNET ?= $(BINDIR)/assocnet $(ASSOCNET_PARAMS)
ASSOCNET_FILTER ?= $(BINDIR)/assocnet-filter $(NET_FILTER_PARAMS)

# Scripts
GEN_KINACT_TBL_SCRIPT = $(SRCDIR)/gen-activity-table.r
POSTERIOR_PROB_SCRIPT = $(SRCDIR)/posterior-prob.r
DISCRETIZE_SCRIPT = $(SRCDIR)/discretize.r
NFCHISQ_SCRIPT = $(SRCDIR)/nfchisq.r
PSSM_SCRIPT = $(SRCDIR)/PSSM.py
AA_FREQS_SCRIPT = $(SRCDIR)/aa-freqs.py
FILTER_PSITE_PLUS_SCRIPT = $(SRCDIR)/filter-psite-plus-tbl.r

# Precious...do not delete
.PRECIOUS: $(DISCR_KINACT_DATA)

# Phony targets
.PHONY: posterior
posterior: $(KINACT_POSTERIOR_PROB)

.PHONY: assoc-results
assoc-results: $(KINACT_ASSOC)

.PHONY: pssm
pssm: $(KIN_KIN_SCORES) $(KIN_KNOWN_PSITE_SCORES) $(KIN_SCORE_DIST)

.PHONY: data
data: $(KINACT_DATA) $(IMP_KINACT_DATA) $(EGF_KIN_ACT_DATA) $(IMP_EGF_KIN_ACT_DATA)

.PHONY: clean-data
clean-data:
	-rm -v $(KINACT_DATA) $(IMP_KINACT_DATA) $(EGF_KIN_ACT_DATA) $(IMP_EGF_KIN_ACT_DATA)

.PHONY: clean-assoc
clean-assoc:
	-rm -v $(KINACT_ASSOC)

.PHONY: clean-posterior
clean-posterior:
	-rm -v $(KINACT_POSTERIOR_PROB)

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
	$(RSCRIPT) $(GEN_KINACT_TBL_SCRIPT)

$(DATADIR)/%-discr.tsv: $(DATADIR)/%-imp.tsv $(DISCRETIZE_SCRIPT)
	$(RSCRIPT) $(DISCRETIZE_SCRIPT) $(DISCR_METHOD) $< $@

$(OUTDIR)/%-pcor.tsv: $(DATADIR)/%-imp.tsv $(OUTDIR)
	$(ASSOCNET) --method=pearson $< >$@

$(OUTDIR)/%-scor.tsv: $(DATADIR)/%-imp.tsv $(OUTDIR)
	$(ASSOCNET) --method=spearman $< >$@

$(OUTDIR)/%-nfchisq.tsv: $(DATADIR)/%-discr.tsv $(NFCHISQ_SCRIPT) $(OUTDIR)
	$(RSCRIPT) $(NFCHISQ_SCRIPT) $< $@

$(OUTDIR)/%-filter.tsv: $(OUTDIR)/%.tsv
	$(ASSOCNET_FILTER) --header-in $< >$@

$(OUTDIR)/%-posterior.tsv: $(OUTDIR)/%.tsv $(POSTERIOR_PROB_SCRIPT) \
							$(KIN_KIN_SCORES) $(KIN_SCORE_DIST)
	$(RSCRIPT) $(POSTERIOR_PROB_SCRIPT) $(ASSOC_METHOD) $< $(KIN_KIN_SCORES) \
		$(KIN_SCORE_DIST) $@

$(HUMAN_KINASE_TABLE): $(FULL_KIN_SUBSTR_TABLE)
	sed '1,3d;4s|+/-|...|' $< | awk -F"\t" '{if (NR==1 || ($$4=="human" && $$9=="human")){print}}' >$@

$(KIN_SUBSTR_TABLE): $(HUMAN_KINASE_TABLE) $(FILTER_PSITE_PLUS_SCRIPT)
	$(RSCRIPT) $(FILTER_PSITE_PLUS_SCRIPT) $< $@

$(PHOSPHOSITES): $(FULL_PHOS_SITES_TABLE) $(FILTER_PSITE_PLUS_SCRIPT)
	sed '1,3d' $< >$@.tmp
	$(RSCRIPT) $(FILTER_PSITE_PLUS_SCRIPT) $@.tmp $@
	rm $@.tmp

$(KIN_KIN_SCORES) \
$(KIN_KNOWN_PSITE_SCORES) \
$(KIN_SCORE_DIST): $(AA_FREQS) $(PSSM_SCRIPT) $(KIN_SUBSTR_TABLE) $(HUMAN_KINASE_TABLE) $(PHOSPHOSITES)
	$(PYTHON2) $(PSSM_SCRIPT)

$(AA_FREQS): $(AA_FREQS_SCRIPT)
	$(PYTHON) $(AA_FREQS_SCRIPT) >$@
