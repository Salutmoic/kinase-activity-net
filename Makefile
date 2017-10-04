#################
### Variables ###
#################

##############
## Directories

SRCDIR = src
DATADIR = data
EXT_DATADIR = $(DATADIR)/external
KEGGDIR = $(EXT_DATADIR)/kegg-relationships
OUTDIR = out
IMGDIR = img
BINDIR = /nfs/research2/beltrao/software-rh7/bin

#############
## Parameters

TABLE_STRATEGIES = max-rows max-cols max-rows-max-cols
TABLE_STRATEGY ?= max-rows-max-cols
ASSOC_METHODS = pcor pcor-filter scor scor-filter nfchisq mut_info fnn_mut_info partcor
ASSOC_METHOD ?= scor
ASSOCNET_FILTER_METHOD ?= deconvolution
ASSOCNET_FILTER_SCALE_METHOD ?= abs
DISCR_METHOD ?= mclust.whole
DECONVOLUTION_A ?= 1.0
DECONVOLUTION_B ?= 0.99

KEGG_VALSET ?= $(KEGG_ACT_VALSET)
VAL_SET = $(KEGG_VALSET) ## $(PSITE_PLUS_VALSET)
# A list of protein annotations used to group proteins so that during
# validation we don't randomly select two kinases that are known to be
# functionally related (e.g. in the same pathway) for a true negative.
PROTEIN_GROUPING = $(KEGG_PATH_REFERENCE)

#####################
## External data sets

# Kinase activity predictions
GSEA_DATA = $(EXT_DATADIR)/log.wGSEA.kinase_condition.clean.Rdata
# Ext phosphosite data
PSITE_DATA = $(EXT_DATADIR)/esetNR.Rdata
# specifies kinases perturbed in each condition
KIN_COND_PAIRS = $(EXT_DATADIR)/kinase-condition-pairs.tsv
# tries to calculate which kinases are perturbed in each condition
KIN_INVIVO_CONDS = $(EXT_DATADIR)/kinase_invivoconditions.tsv 
# PhosphositePlus data
PHOSPHOSITE_PLUS_VERSION = 2017-09-08
FULL_KIN_SUBSTR_TABLE = $(EXT_DATADIR)/Kinase_Substrate_Dataset_$(PHOSPHOSITE_PLUS_VERSION)
FULL_REG_SITES_TABLE = $(EXT_DATADIR)/Regulatory_sites_$(PHOSPHOSITE_PLUS_VERSION)
FULL_PHOS_SITES_TABLE = $(EXT_DATADIR)/Phosphorylation_site_dataset_$(PHOSPHOSITE_PLUS_VERSION)
# Kegg data
KEGG_RELATIONSHIPS = $(wildcard $(KEGGDIR)/*.ppip)
# Uniprot data
FULL_UNIPROT_ID_MAPPING = $(EXT_DATADIR)/HUMAN_9606_idmapping.dat

######################
## Generated data sets

# Initial kinase-activity tables, trying to
# maximize #rows, #cols or both.  Raw and imputed datasets
KINACT_DATA = $(DATADIR)/kinact-$(TABLE_STRATEGY).tsv
IMP_KINACT_DATA = $(DATADIR)/kinact-$(TABLE_STRATEGY)-imp.tsv
DISCR_KINACT_DATA = $(DATADIR)/kinact-$(TABLE_STRATEGY)-discr.tsv
# EGF-signaling-related kinase activity
EGF_KIN_ACT_DATA = $(foreach T,perturbed unperturbed,DATADIR)/egf-kinact-$(T).tsv)
IMP_EGF_KIN_ACT_DATA = $(foreach T,perturbed unperturbed,$(DATADIR)/egf-kinact-$(T)-imp.tsv)
# Association table output
KINACT_ASSOC = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-$(ASSOC_METHOD).tsv
# Final, merged predictor
PREDICTOR = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-$(ASSOC_METHOD)-final-predictor.tsv
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
HUMAN_REG_SITES = $(DATADIR)/human_reg_sites.tsv
REG_SITES = $(DATADIR)/reg_sites.tsv
# Kegg-derived files
KEGG_ACT = $(DATADIR)/kegg-act-rels.tsv
KEGG_PHOS_ACT = $(DATADIR)/kegg-phos-act-rels.tsv
KEGG_PHOS = $(DATADIR)/kegg-phos-rels.tsv
KEGG_PATH_REFERENCE = $(DATADIR)/kegg-path-ref.tsv
# Uniprot-derived files
UNIPROT_ID_MAPPING = $(DATADIR)/uniprot-id-map.tsv
# Validation sets
PSITE_PLUS_VALSET = $(DATADIR)/validation-set-psiteplus.tsv
KEGG_ACT_VALSET = $(DATADIR)/validation-set-kegg-act.tsv
KEGG_PHOS_ACT_VALSET = $(DATADIR)/validation-set-kegg-phos-act.tsv
KEGG_PHOS_VALSET = $(DATADIR)/validation-set-kegg-phos.tsv

#########
## Images

# Validation
ASSOC_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-$(ASSOC_METHOD)-val.pdf
PSSM_VAL_IMG = $(IMGDIR)/kinase_kinase_scores-val.pdf
PREDICTOR_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-$(ASSOC_METHOD)-final-predictor-val.pdf
VAL_IMGS = $(ASSOC_VAL_IMG) $(PSSM_VAL_IMG) $(PREDICTOR_VAL_IMG)

##################
## Program Options

ASSOCNET_PARAMS = --unbiased-correlation --p-method=none
ASSOCNET_FILTER_PARAMS = --method=$(ASSOCNET_FILTER_METHOD) --scale-method=$(ASSOCNET_FILTER_SCALE_METHOD) --header-in
ifeq ($(NET_FILTER_METHOD),deconvolution)
ASSOCNET_FILTER_PARAMS += --deconvolution-a=$(DECONVOLUTION_A)	\
						  --deconvolution-b=$(DECONVOLUTION_B)
endif

###########
## Programs

PYTHON3 ?= $(BINDIR)/python3
PYTHON2 ?= $(BINDIR)/python2
PYTHON ?= $(PYTHON3)
RSCRIPT ?= $(BINDIR)/Rscript
ASSOCNET ?= $(BINDIR)/assocnet $(ASSOCNET_PARAMS)
ASSOCNET_FILTER ?= $(BINDIR)/assocnet-filter $(NET_FILTER_PARAMS)

##########
## Scripts

GEN_KINACT_TBL_SCRIPT = $(SRCDIR)/gen-activity-table.r
FINAL_PREDICTOR_SCRIPT = $(SRCDIR)/final-predictor.r
DISCRETIZE_SCRIPT = $(SRCDIR)/discretize.r
NFCHISQ_SCRIPT = $(SRCDIR)/nfchisq.r
ASSOC_SCRIPT = $(SRCDIR)/assoc_methods.r
PSSM_SCRIPT = $(SRCDIR)/PSSM.py
AA_FREQS_SCRIPT = $(SRCDIR)/aa-freqs.py
FILTER_PSITE_PLUS_SCRIPT = $(SRCDIR)/filter-psite-plus-tbl.r
PSITE_PLUS_VAL_SCRIPT = $(SRCDIR)/make-psiteplus-valset.r
FILTER_KEGG_RELS_SCRIPT = $(SRCDIR)/filter-kegg-tbl.r
KEGG_PATH_REF_SCRIPT = $(SRCDIR)/gen-kegg-pathway-ref.py
KEGG_VAL_SCRIPT = $(SRCDIR)/make-kegg-valset.r
VAL_SCRIPT = $(SRCDIR)/validation.r

# Precious...do not delete
.PRECIOUS: $(DISCR_KINACT_DATA)

#####################
### Phony targets ###
#####################

.PHONY: final-predictor
final-predictor: $(PREDICTOR)

.PHONY: assoc-results
assoc-results: $(KINACT_ASSOC)

.PHONY: pssm
pssm: $(KIN_KIN_SCORES) $(KIN_KNOWN_PSITE_SCORES) $(KIN_SCORE_DIST)

.PHONY: data
data: $(KINACT_DATA) $(IMP_KINACT_DATA) $(EGF_KIN_ACT_DATA) $(IMP_EGF_KIN_ACT_DATA)

.PHONY: validation
validation: $(VAL_IMGS)

.PHONY: clean-data
clean-data:
	-rm -v $(KINACT_DATA) $(IMP_KINACT_DATA) $(EGF_KIN_ACT_DATA) $(IMP_EGF_KIN_ACT_DATA)

.PHONY: clean-assoc
clean-assoc:
	-rm -v $(KINACT_ASSOC)

.PHONY: clean-final-predictor
clean-final-predictor:
	-rm -v $(PREDICTOR)

.PHONY: clean-validation
clean-validation:
	-rm -v $(VAL_IMGS)

.PHONY: clean-results
clean-results: clean-assoc clean-final-predictor clean-validation

.PHONY: clean
clean: clean-data clean-results

#############
### Rules ###
#############

##############
## Directories

$(OUTDIR):
	mkdir -p $(OUTDIR)

$(IMGDIR):
	mkdir -p $(IMGDIR)

########################
## Processed data tables

# Kinase activity tables
$(KINACT_DATA) \
$(IMP_KINACT_DATA): $(GEN_KINACT_TBL_SCRIPT) $(GSEA_DATA) $(KIN_COND_PAIRS)
	$(RSCRIPT) $(GEN_KINACT_TBL_SCRIPT) $(TABLE_STRATEGY)

# Discretized kinase activity
$(DATADIR)/%-discr.tsv: $(DATADIR)/%-imp.tsv $(DISCRETIZE_SCRIPT)
	$(RSCRIPT) $(DISCRETIZE_SCRIPT) $(DISCR_METHOD) $< $@

# Grab just the human kinase-substrate relationships
$(HUMAN_KINASE_TABLE): $(FULL_KIN_SUBSTR_TABLE)
	sed '1,3d;4s|+/-|...|' $< | awk -F"\t" '{if (NR==1 || ($$4=="human" && $$9=="human")){print}}' >$@

# Filter the human kinase-substrates to just those for which we have
# kinase activities
$(KIN_SUBSTR_TABLE): $(HUMAN_KINASE_TABLE) $(FILTER_PSITE_PLUS_SCRIPT)
	$(RSCRIPT) $(FILTER_PSITE_PLUS_SCRIPT) $< $@

# Filter the PhosphoSitePlus site list to just those on kinases for
# which we have kinase activities
$(PHOSPHOSITES): $(FULL_PHOS_SITES_TABLE) $(FILTER_PSITE_PLUS_SCRIPT)
	sed '1,3d' $< >$@.tmp
	$(RSCRIPT) $(FILTER_PSITE_PLUS_SCRIPT) $@.tmp $@
	rm $@.tmp

# Filter the Regulatory sites table to have just human phosphosites
$(HUMAN_REG_SITES): $(FULL_REG_SITES_TABLE)
	sed '1,3d' $< | awk -F"\t" 'BEGIN{OFS="\t"}{if (NR==1 || $$7=="human" && $$8~/-p$$/){print}}' | sed 's/-p//' >$@

# Further filter to just include reg sites on kinases for which we
# have kinase activities
$(REG_SITES): $(HUMAN_REG_SITES) $(FILTER_PSITE_PLUS_SCRIPT)
	$(RSCRIPT) $(FILTER_PSITE_PLUS_SCRIPT) $< $@

# Combine PhosphositePlus kinase-substrate and regulatory sites tables
# to generate a set of true-positive relationships
$(PSITE_PLUS_VALSET): $(KIN_SUBSTR_TABLE) $(REG_SITES)
	$(RSCRIPT) $(PSITE_PLUS_VAL_SCRIPT) $^

# Get all activation relationships from Kegg
$(KEGG_ACT): $(KEGG_RELATIONSHIPS)
	cat $^ | sed -E -n '/activation|inhibition|repression/p' | \
		awk 'BEGIN{OFS="\t"}{print $$1, $$2}' | \
		sort | uniq >$@

# Get all activation relationships from Kegg
$(KEGG_PHOS_ACT): $(KEGG_RELATIONSHIPS)
	cat $^ | sed -E -n '/activation|inhibition|repression/p' | \
		sed -n '/phosphorylation/p' | \
		awk 'BEGIN{OFS="\t"}{print $$1, $$2}' | \
		sort | uniq >$@

# Get all phospho activation relationships from Kegg
$(KEGG_PHOS): $(KEGG_RELATIONSHIPS)
	cat $^ | sed -n '/phosphorylation/p' | \
		awk 'BEGIN{OFS="\t"}{print $$1, $$2}' | \
		sort | uniq >$@

# Create the Uniprot ID map
$(UNIPROT_ID_MAPPING): $(FULL_UNIPROT_ID_MAPPING)
	awk 'BEGIN{OFS="\t"}{if ($$2=="Gene_Name"){print $$1, $$3}}' $< >$@

# Kegg pathway reference
$(KEGG_PATH_REFERENCE): $(KEGG_RELATIONSHIPS) $(UNIPROT_ID_MAPPING)
	$(PYTHON) $(KEGG_PATH_REF_SCRIPT) | sort | uniq >$@

# Kegg-based validation sets
$(DATADIR)/validation-set-kegg-%.tsv: $(DATADIR)/kegg-%-rels.tsv $(UNIPROT_ID_MAPPING) $(KEGG_VAL_SCRIPT)
	$(RSCRIPT) $(KEGG_VAL_SCRIPT) $(wordlist 1,2,$^) $@

# Calculate human proteome amino-acid frequencies
$(AA_FREQS): $(AA_FREQS_SCRIPT)
	$(PYTHON) $(AA_FREQS_SCRIPT) >$@

###########################
## Association score tables

# Pearson's correlation
$(OUTDIR)/%-pcor.tsv: $(DATADIR)/%-imp.tsv $(OUTDIR)
	$(ASSOCNET) --method=pearson $< >$@

# Spearman's correlation
$(OUTDIR)/%-scor.tsv: $(DATADIR)/%-imp.tsv $(OUTDIR)
	$(ASSOCNET) --method=spearman $< >$@

# Normalized FunChisq
$(OUTDIR)/%-nfchisq.tsv: $(DATADIR)/%-discr.tsv $(NFCHISQ_SCRIPT) $(OUTDIR)
	$(RSCRIPT) $(NFCHISQ_SCRIPT) $< $@
# Here I create three different correlation tables with the same script assoc_methods.r

# Partial correlation
$(OUTDIR)/%-partcor.tsv: $(DATADIR)/%-imp.tsv $(ASSOC_SCRIPT) $(OUTDIR)
	$(RSCRIPT) $(ASSOC_SCRIPT) partcor $< $@	

# Mutual information
$(OUTDIR)/%-mut_info.tsv: $(DATADIR)/%-discr.tsv $(ASSOC_SCRIPT) $(OUTDIR)
	$(RSCRIPT) $(ASSOC_SCRIPT) mut_info $< $@

# Fast-nearest-neighbors (FNN) mutual information
$(OUTDIR)/%-fnn_mut_info.tsv: $(DATADIR)/%-imp.tsv $(ASSOC_SCRIPT) $(OUTDIR)
	$(RSCRIPT) $(ASSOC_SCRIPT) fnn_mut_info $< $@

# Filter out indirect associations
$(OUTDIR)/%-filter.tsv: $(OUTDIR)/%.tsv
	$(ASSOCNET_FILTER) --header-in $< >$@

###############################
## Kinase-substrate predictions

# Calculate kinase-substrate PSSM scores
$(KIN_KIN_SCORES) \
$(KIN_KNOWN_PSITE_SCORES) \
$(KIN_SCORE_DIST): $(AA_FREQS) $(PSSM_SCRIPT) $(KIN_SUBSTR_TABLE) $(HUMAN_KINASE_TABLE) $(PHOSPHOSITES)
	$(PYTHON2) $(PSSM_SCRIPT)

###################
## Merged predictor

# Calculate a merged/posterior probability given all evidence.
$(PREDICTOR): $(KINACT_ASSOC) $(FINAL_PREDICTOR_SCRIPT) \
			  $(KIN_KIN_SCORES) $(KIN_SCORE_DIST)
	$(RSCRIPT) $(FINAL_PREDICTOR_SCRIPT) $(ASSOC_METHOD) $< $(KIN_KIN_SCORES) \
		$(KIN_SCORE_DIST) $@

#############
## Validation

$(IMGDIR)/%-val.pdf: $(OUTDIR)/%.tsv $(VAL_SET) $(PROTEIN_GROUPING) $(IMGDIR) $(VAL_SCRIPT)
	$(RSCRIPT) $(VAL_SCRIPT) $(wordlist 1,3,$^)

