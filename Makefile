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
TMPDIR = tmp
ifeq ($(KSEA_USE_AUTOPHOS),FALSE)
	KSEA_TMPDIR = $(TMPDIR)/ksea
else
	KSEA_TMPDIR = $(TMPDIR)/ksea-autophos
endif
BINDIR = /nfs/research2/beltrao/software-rh7/bin

#############
## Parameters

KSEA_NUM_CONDS = 668
KSEA_USE_AUTOPHOS ?= TRUE
TABLE_STRATEGIES = max-rows max-cols balanced
TABLE_STRATEGY ?= max-rows
KSEA_MIN_SITES ?= 3
ENTROPY_FILTER ?= 0.1 			# Only keep rows/columns with at least 0.X*(max row/col entropy)
NA_THRESHOLD ?= 0.3333
ASSOC_METHODS = pcor pcor-filter scor scor-filter nfchisq mut_info	\
	fnn_mut_info partcor all fvalue paircor
ASSOC_METHOD ?= scor
DISCR_METHOD ?= trunc 	# mclust.whole mclust.by.row manual trunc
ASSOCNET_FILTER_METHOD ?= deconvolution
ASSOCNET_FILTER_PRESCALE_METHOD ?= none
ASSOCNET_FILTER_POSTSCALE_METHOD ?= standard
# Network density post-deconvolution
DECONVOLUTION_A ?= 0.5
# Eigenvalue Scaling parameter
DECONVOLUTION_B ?= 0.99
PRED_METHOD = log_reg

KEGG_VALSET ?= $(KEGG_ACT_VALSET)
VAL_SET = $(KEGG_VALSET) ## $(PSITE_PLUS_VALSET)
# A list of protein annotations used to group proteins so that during
# validation we don't randomly select two kinases that are known to be
# functionally related (e.g. in the same pathway) for a true negative.
PROTEIN_GROUPING =  $(COMBINED_GROUPING) # $(COMBINED_GROUPING) $(KEGG_PATH_REFERENCE) $(GO_CELL_LOCATION)

MERGED_PRED_SOURCES = $(KINACT_ASSOC) $(KINACT_ASSOC2)			\
	$(REG_SITE_ASSOC) $(KIN_KIN_SCORES) $(HUMAN_STRING_COEXP)	\
	$(HUMAN_STRING_COEXP2)
MERGED_PRED_SOURCES_KINOME = $(KINACT_ASSOC) $(KINACT_ASSOC2)	\
	$(KIN_KIN_SCORES_KINOME) $(HUMAN_STRING_COEXP_KINOME)		\
	$(HUMAN_STRING_COEXP2_KINOME)

#####################
## External data sets

# Ext phosphosite data
PSITE_DATA = $(EXT_DATADIR)/esetNR.Rdata
# specifies kinases perturbed in each condition
KIN_COND_PAIRS = $(EXT_DATADIR)/kinase-condition-pairs.tsv
# tries to calculate which kinases are perturbed in each condition
KIN_INVIVO_CONDS = $(EXT_DATADIR)/kinase_invivoconditions.csv 
# PhosphositePlus data
PHOSPHOSITE_PLUS_VERSION = 2017-09-08
FULL_KIN_SUBSTR_TABLE = $(EXT_DATADIR)/Kinase_Substrate_Dataset_$(PHOSPHOSITE_PLUS_VERSION)
FULL_REG_SITES_TABLE = $(EXT_DATADIR)/Regulatory_sites_$(PHOSPHOSITE_PLUS_VERSION)
FULL_PHOS_SITES_TABLE = $(EXT_DATADIR)/Phosphorylation_site_dataset_$(PHOSPHOSITE_PLUS_VERSION)
# Kegg data
KEGG_RELATIONSHIPS = $(wildcard $(KEGGDIR)/*.ppip)
# Uniprot data
FULL_UNIPROT_ID_MAPPING = $(EXT_DATADIR)/HUMAN_9606_idmapping.dat
# Gene Ontology
GO_OBO_VERSION = 2017-10-04
GO_OBO = $(EXT_DATADIR)/go-basic-$(GO_OBO_VERSION).obo
HUMAN_GO_VERSION = 2017-09-26
FULL_GO_ASSOC_TABLE = $(EXT_DATADIR)/goa_human-$(HUMAN_GO_VERSION).gaf
STRING_VERSION = 10.5
HUMAN_STRING_RAW = $(EXT_DATADIR)/9606.protein.links.detailed.v$(STRING_VERSION).txt
# Kinase beads
KINASE_BEADS_RAW = $(EXT_DATADIR)/kinase-beads-activities.tsv
# Kinome
KINOME_RAW = $(EXT_DATADIR)/pkinfam.txt
# Kinase domain phosphorylation hot-spots from Marta
HOTSPOTS_RAW = $(EXT_DATADIR)/kinase-phospho-hotspots.tsv
TYR_HOTSPOTS_RAW = $(EXT_DATADIR)/tyr-kinase-phospho-hotspots.tsv

######################
## Generated data sets

# Kinase activity predictions
KSEA_TMP_DATA = $(foreach COND_NUM,$(shell seq 1 $(KSEA_NUM_CONDS)),$(KSEA_TMPDIR)/cond-$(COND_NUM).tsv)
KSEA_DATA = $(DATADIR)/log.wKSEA.kinase_condition.clean.Rdata
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
KINACT_ASSOC2 = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-$(ASSOC_METHOD)2.tsv
# Regulatory site-based activity correlations
REG_SITE_ASSOC = $(OUTDIR)/reg-site-cor.tsv
# Final, merged predictor
PREDICTOR = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-$(ASSOC_METHOD)-final-predictor.tsv
# Human amino acid frequencies
AA_FREQS = $(DATADIR)/aa-freqs.tsv
# PSSM files
KIN_KIN_SCORES = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-pssm.tsv
KIN_KIN_SCORES_KINOME = $(OUTDIR)/human-kinome-pssm.tsv
KIN_KNOWN_PSITE_SCORES = $(OUTDIR)/known_kinase_psite-score
KIN_SCORE_DIST = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-pssm-dists.tsv
# PhosphositePlus derived files
KIN_SUBSTR_TABLE = $(DATADIR)/reduced_kinase_table.tsv
HUMAN_KINASE_TABLE = $(DATADIR)/human_kinase_table.tsv
PHOSPHOSITES = $(DATADIR)/phosphosites_reduced.tsv
HUMAN_REG_SITES = $(DATADIR)/human_reg_sites.tsv
REG_SITES = $(DATADIR)/reg_sites.tsv
KIN_SUB_OVERLAP = $(DATADIR)/kinase_substrate_overlap.tsv
# Kegg-derived files
KEGG_ACT = $(DATADIR)/kegg-act-rels.tsv
KEGG_PHOS_ACT = $(DATADIR)/kegg-phos-act-rels.tsv
KEGG_PHOS = $(DATADIR)/kegg-phos-rels.tsv
KEGG_PATH_REFERENCE = $(DATADIR)/kegg-path-ref.tsv
# Uniprot-derived files
UNIPROT_ID_MAPPING = $(DATADIR)/uniprot-id-map.tsv
ENSEMBL_ID_MAPPING = $(DATADIR)/ensembl-id-map.tsv
# GO-derived files
GO_ASSOC = $(DATADIR)/go-assoc.tsv
GO_CELL_LOCATION = $(DATADIR)/go-cell-location.tsv
# STRING (in the output directory for convenience of using it as a
# predictor)
STRING_ID_MAPPING = $(DATADIR)/string-id-map.tsv
HUMAN_STRING_COEXP = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-string-coexp.tsv
HUMAN_STRING_COEXP2 = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-string-coexp2.tsv
HUMAN_STRING_COEXP2_FULL = $(OUTDIR)/string-coexp2-full.tsv
HUMAN_STRING_COOCC = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-string-coocc.tsv
HUMAN_STRING_COOCC2 = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-string-coocc2.tsv
HUMAN_STRING_COOCC2_FULL = $(OUTDIR)/string-coocc2-full.tsv
HUMAN_STRING_EXPER = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-string-exper.tsv
HUMAN_STRING_COEXP_KINOME = $(OUTDIR)/human-kinome-string-coexp.tsv
HUMAN_STRING_COEXP2_KINOME = $(OUTDIR)/human-kinome-string-coexp2.tsv
HUMAN_STRING_EXPER_KINOME = $(OUTDIR)/human-kinome-string-exper.tsv
HUMAN_STRING_COOCC_KINOME = $(OUTDIR)/human-kinome-string-coocc.tsv
HUMAN_STRING_COOCC2_KINOME = $(OUTDIR)/human-kinome-string-coocc2.tsv
HUMAN_STRING_EXPER_KINOME = $(OUTDIR)/human-kinome-string-exper.tsv
# Protein groupings
COMBINED_GROUPING = $(DATADIR)/protein-groups.tsv
# Human kinome
HUMAN_KINOME = $(DATADIR)/human-kinome.txt
# Kinase-domain phosphorylation hot spots
HOTSPOTS = $(DATADIR)/kinase-phospho-hotspots.tsv
# Validation sets
PSITE_PLUS_VALSET = $(DATADIR)/validation-set-psiteplus.tsv
KEGG_ACT_VALSET = $(DATADIR)/validation-set-kegg-act.tsv
KEGG_PHOS_ACT_VALSET = $(DATADIR)/validation-set-kegg-phos-act.tsv
KEGG_PHOS_VALSET = $(DATADIR)/validation-set-kegg-phos.tsv
NEG_VALSET = $(DATADIR)/validation-set-negative.tsv
# Kinase beads activities
KINASE_BEADS = $(DATADIR)/kinase-beads-activities.tsv
# Merged predictor data
MERGED_PRED = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-merged.tsv
MERGED_PRED_KINOME = $(OUTDIR)/human-kinome-merged.tsv
# Final predictor
FINAL_PRED = $(OUTDIR)/kinact-$(TABLE_STRATEGY)-$(PRED_METHOD).tsv

#########
## Images

# Validation
ASSOC_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-$(ASSOC_METHOD)-val.pdf
ASSOC2_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-$(ASSOC_METHOD)2-val.pdf
PSSM_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-pssm-val.pdf
STRING_COEXP_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-string-coexp-val.pdf
STRING_COEXP2_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-string-coexp2-val.pdf
STRING_COOCC_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-string-coocc-val.pdf
STRING_COOCC2_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-string-coocc2-val.pdf
STRING_EXPER_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-string-exper-val.pdf
PREDICTOR_VAL_IMG = $(IMGDIR)/kinact-$(TABLE_STRATEGY)-$(ASSOC_METHOD)-final-predictor-val.pdf
VAL_IMGS = $(ASSOC_VAL_IMG) $(ASSOC2_VAL_IMG) $(PSSM_VAL_IMG)	\
	$(STRING_COEXP_VAL_IMG) $(STRING_COOCC_VAL_IMG)				\
	$(STRING_EXPER_VAL_IMG) # $(PREDICTOR_VAL_IMG)

##################
## Program Options

ASSOCNET_PARAMS = --unbiased-correlation --p-method=none
ASSOCNET_FILTER_PARAMS = --method=$(ASSOCNET_FILTER_METHOD) \
						--pre-scale-method=$(ASSOCNET_FILTER_PRESCALE_METHOD) \
						--post-scale-method=$(ASSOCNET_FILTER_POSTSCALE_METHOD) \
						--header-in --observed-only
ifeq ($(ASSOCNET_FILTER_METHOD),deconvolution)
ASSOCNET_FILTER_PARAMS += --deconvolution-a=$(DECONVOLUTION_A)	\
						  --deconvolution-b=$(DECONVOLUTION_B)
endif

###########
## Programs

BASH = /bin/bash
SHELL = $(BASH)
PYTHON3 ?= $(BINDIR)/python3
PYTHON2 ?= $(BINDIR)/python2
PYTHON ?= $(PYTHON3)
RSCRIPT ?= $(BINDIR)/Rscript
ASSOCNET ?= $(BINDIR)/assocnet $(ASSOCNET_PARAMS)
ASSOCNET_FILTER ?= $(BINDIR)/assocnet-filter $(ASSOCNET_FILTER_PARAMS)

##########
## Scripts

KSEA_POST_SCRIPT = $(SRCDIR)/ksea-post.r
GEN_KINACT_TBL_SCRIPT = $(SRCDIR)/gen-activity-table.r
FINAL_PREDICTOR_SCRIPT = $(SRCDIR)/final-predictor.r
DISCRETIZE_SCRIPT = $(SRCDIR)/discretize.r
NFCHISQ_SCRIPT = $(SRCDIR)/nfchisq.r
ASSOC_SCRIPT = $(SRCDIR)/assoc_methods.r
LM_PRED_SCRIPT = $(SRCDIR)/lm-pred.r
PAIR_COR_SCRIPT = $(SRCDIR)/pairwise-cor.r
PSSM_SCRIPT = $(SRCDIR)/pssm.py
PSSM_DIST_SCRIPT = $(SRCDIR)/pssm-dist.py
PSSM_LIB = $(SRCDIR)/pssms/__init__.py
AA_FREQS_SCRIPT = $(SRCDIR)/aa-freqs.py
FILTER_PSITE_PLUS_SCRIPT = $(SRCDIR)/filter-psite-plus-tbl.r
KIN_SUB_OVERLAP_SCRIPT = $(SRCDIR)/kinase-overlap.r
PSITE_PLUS_VAL_SCRIPT = $(SRCDIR)/make-psiteplus-valset.r
FILTER_KEGG_RELS_SCRIPT = $(SRCDIR)/filter-kegg-tbl.r
KEGG_PATH_REF_SCRIPT = $(SRCDIR)/gen-kegg-pathway-ref.py
KEGG_VAL_SCRIPT = $(SRCDIR)/make-kegg-valset.r
NEG_VAL_SCRIPT = $(SRCDIR)/make-negative-valset.r
REFORMAT_GO_ASSOC_SCRIPT = $(SRCDIR)/reformat-go-assoc.awk
GO_CELL_LOC_SCRIPT = $(SRCDIR)/get-go-cell-component.py
FORMAT_STRING_SCRIPT = $(SRCDIR)/format-string.py
VAL_SCRIPT = $(SRCDIR)/validation.r
MERGE_SCRIPT = $(SRCDIR)/merge.r
FILTER_BY_PROTLIST_SCRIPT = $(SRCDIR)/filter-by-protein-list.r
REG_SITE_COR_SCRIPT = $(SRCDIR)/reg-site-cor.r

# Don't delete intermediate files
.SECONDARY:

#####################
### Phony targets ###
#####################

.PHONY: final-predictor
final-predictor: $(PREDICTOR)

.PHONY: merge-predictors
merge-predictors: $(MERGED_PRED)

.PHONY: assoc-results
assoc-results: $(KINACT_ASSOC) $(KINACT_ASSOC2)

.PHONY: pssm
pssm: $(KIN_KIN_SCORES) $(KIN_KNOWN_PSITE_SCORES) $(KIN_SCORE_DIST)

.PHONY: data
data: $(KSEA_DATA) $(KINACT_DATA) $(IMP_KINACT_DATA)

.PHONY: validation
validation: $(VAL_IMGS)

.PHONY: merge-validation
merge-validation: $(IMGDIR)/validation.pdf

.PHONY: clean-data
clean-data:
	-rm -v $(KINACT_DATA) $(IMP_KINACT_DATA)

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

.PHONY: clean-logs
clean-logs:
	-rm -v log/*

.PHONY: clean
clean: clean-data clean-results

# For debugging purposes
expand-var-%:
	@echo $($*)

.PHONY: info
info:
	@printf "\033[1mTable strategy:\033[0m\t"
	@case "$(TABLE_STRATEGY)" in \
		"max-rows" ) printf "max. kinases\n" ;; \
		"max-cols" ) printf "max. conditions\n" ;; \
		"balanced" ) printf "balanced\n" ;; \
	esac
	@printf "\033[1mTable dim.:\033[0m\t"
	@printf "%s x %s\n" `sed '1d' $(KINACT_DATA) | cut -f1 | sort | uniq | wc -l` \
		`sed '1d' $(KINACT_DATA) | cut -f2 | sort | uniq | wc -l`
	@printf "\033[1m%% missing:\033[0m\t"
	@awk 'BEGIN{na=0}{if ($$3=="NA"){na=na+1}}END{printf("%f\n", na/(NR-1))}' $(KINACT_DATA)


#############
### Rules ###
#############

##############
## Directories

$(DATADIR):
	mkdir -p $(DATADIR)

$(OUTDIR):
	mkdir -p $(OUTDIR)

$(IMGDIR):
	mkdir -p $(IMGDIR)

########################
## Processed data tables

# KSEA results
# The KSEA temporary files are created outside of this Makefile using
# the do-ksea.sh script
$(KSEA_DATA): $(KSEA_TMP_DATA)
	cat $^ >$@.tmp
	$(RSCRIPT) $(KSEA_POST_SCRIPT) $@.tmp
	rm $@.tmp

# Kinase activity tables
$(KINACT_DATA) \
$(IMP_KINACT_DATA): $(GEN_KINACT_TBL_SCRIPT) $(KSEA_DATA) $(KIN_COND_PAIRS) \
					$(KIN_SUB_OVERLAP)
	$(RSCRIPT) $(GEN_KINACT_TBL_SCRIPT) $(TABLE_STRATEGY) $(KSEA_MIN_SITES) $(ENTROPY_FILTER) $(NA_THRESHOLD)

# Discretized kinase activity
$(DATADIR)/%-discr.tsv: $(DATADIR)/%-imp.tsv $(DISCRETIZE_SCRIPT)
	$(RSCRIPT) $(DISCRETIZE_SCRIPT) $(DISCR_METHOD) $< $@

# Grab just the human kinase-substrate relationships
$(HUMAN_KINASE_TABLE): $(FULL_KIN_SUBSTR_TABLE)
	sed '1,3d;4s|+/-|...|' $< | awk -F"\t" '{if (NR==1 || ($$4=="human" && $$9=="human")){print}}' >$@

# Filter the human kinase-substrates to just those for which we have
# kinase activities
$(KIN_SUBSTR_TABLE): $(HUMAN_KINASE_TABLE) $(FILTER_PSITE_PLUS_SCRIPT) $(KSEA_DATA)
	$(RSCRIPT) $(FILTER_PSITE_PLUS_SCRIPT) $< $@

# Filter the PhosphoSitePlus site list to just those on kinases for
# which we have kinase activities
$(PHOSPHOSITES): $(FULL_PHOS_SITES_TABLE) $(FILTER_PSITE_PLUS_SCRIPT) $(KSEA_DATA)
	sed '1,3d' $< | awk -F"\t" -vOFS="\t" '{if (NR == 1 || $$7 == "human"){print}}' >$@.tmp
	$(RSCRIPT) $(FILTER_PSITE_PLUS_SCRIPT) $@.tmp $@
	rm $@.tmp

# Filter the Regulatory sites table to have just human phosphosites
$(HUMAN_REG_SITES): $(FULL_REG_SITES_TABLE)
	sed '1,3d' $< | awk -F"\t" 'BEGIN{OFS="\t"}{if (NR==1 || $$7=="human" && $$8~/-p$$/){print}}' | sed 's/-p//' >$@

# Further filter to just include reg sites on kinases for which we
# have kinase activities
$(REG_SITES): $(HUMAN_REG_SITES) $(FILTER_PSITE_PLUS_SCRIPT) $(KSEA_DATA)
	$(RSCRIPT) $(FILTER_PSITE_PLUS_SCRIPT) $< $@

# Kinase-substrate overlap
$(KIN_SUB_OVERLAP): $(HUMAN_KINASE_TABLE) $(KSEA_DATA) $(KIN_SUB_OVERLAP_SCRIPT)
	$(RSCRIPT) $(KIN_SUB_OVERLAP_SCRIPT)

# Combine PhosphositePlus kinase-substrate and regulatory sites tables
# to generate a set of true-positive relationships
$(PSITE_PLUS_VALSET): $(KIN_SUBSTR_TABLE) $(REG_SITES)
	$(RSCRIPT) $(PSITE_PLUS_VAL_SCRIPT) $^

# Get all activation relationships from Kegg
$(KEGG_ACT): $(KEGG_RELATIONSHIPS)
	cat $^ | sed -E -n '/activation|inhibition|repression/p' | \
		awk 'BEGIN{OFS="\t"}{if ($$1!=$$2){print $$1, $$2}}' | \
		sort | uniq >$@

# Get all activation relationships from Kegg
$(KEGG_PHOS_ACT): $(KEGG_RELATIONSHIPS)
	cat $^ | sed -E -n '/activation|inhibition|repression/p' | \
		sed -n '/phosphorylation/p' | \
		awk 'BEGIN{OFS="\t"}{if ($$1!=$$2){print $$1, $$2}}' | \
		sort | uniq >$@

# Get all phospho activation relationships from Kegg
$(KEGG_PHOS): $(KEGG_RELATIONSHIPS)
	cat $^ | sed -n '/phosphorylation/p' | \
		awk 'BEGIN{OFS="\t"}{if ($$1!=$$2){print $$1, $$2}}' | \
		sort | uniq >$@

# Create the Uniprot ID map
$(UNIPROT_ID_MAPPING): $(FULL_UNIPROT_ID_MAPPING)
	awk 'BEGIN{OFS="\t"}{if ($$2=="Gene_Name"){print $$1, $$3}}' $< | sort >$@

$(ENSEMBL_ID_MAPPING): $(FULL_UNIPROT_ID_MAPPING) $(UNIPROT_ID_MAPPING)
	awk 'BEGIN{OFS="\t"}{if ($$2=="Ensembl_PRO" && $$1 !~ /.-[0-9]/){print $$1, $$3}}' $< | sort -d -k1 >$@.tmp
	join -t'	' $@.tmp $(UNIPROT_ID_MAPPING) | cut -f2,3 | sort >$@
	rm $@.tmp

# GO associations (formatted for GOATOOLS)
$(GO_ASSOC): $(FULL_GO_ASSOC_TABLE) $(UNIPROT_ID_MAPPING) $(REFORMAT_GO_ASSOC_SCRIPT)
	sed '/^!/d' $< | cut -f2,5 | sort | uniq | awk -f $(REFORMAT_GO_ASSOC_SCRIPT) | sort -k1 >$@.tmp
	join -t'	' $(UNIPROT_ID_MAPPING) $@.tmp | cut -f2,3 >$@
	rm $@.tmp

# Parse protein cellular locations
# GO:0005634 = Nucleus
# GO:0005737 = Cytoplasm
# GO:0005886 = Plasma Membrane
$(GO_CELL_LOCATION): $(GO_ASSOC) $(GO_OBO) $(GO_CELL_LOC_SCRIPT)
	$(PYTHON) $(GO_CELL_LOC_SCRIPT) $(GO_ASSOC) $(GO_OBO) | sed '1,2d' | sort | uniq >$@

# Kinome
$(HUMAN_KINOME): $(KINOME_RAW) $(UNIPROT_ID_MAPPING)
	sed -n '/HUMAN/{s/  */\t/g;p}' $< | cut -f3 | sed 's/(//' | \
		grep -f - $(UNIPROT_ID_MAPPING) | cut -f2 | sort >$@

# STRING
# Create the STRING ID map
$(STRING_ID_MAPPING): $(FULL_UNIPROT_ID_MAPPING) $(UNIPROT_ID_MAPPING)
	awk 'BEGIN{OFS="\t"}{if ($$2=="STRING" && $$1 !~ /.-[0-9]/){print $$1, $$3}}' $< | sort -d -k1 >$@.tmp
	join -t'	' $@.tmp $(UNIPROT_ID_MAPPING) | cut -f2,3 >$@
	rm $@.tmp

# STRING coexpression
$(HUMAN_STRING_COEXP): $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) $(KINACT_DATA) \
		$(FORMAT_STRING_SCRIPT)
	$(PYTHON) $(FORMAT_STRING_SCRIPT) $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) \
		<(cut -f1 $(KINACT_DATA) | sort | uniq) 6 >$@.tmp
	cat <(sed -n '1p' $@.tmp) <(sed '1d' $@.tmp | sort -k 1d,1 -k 2d,2) >$@
	rm $@.tmp

$(HUMAN_STRING_COEXP2_FULL): $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) $(KINACT_DATA)
	$(PYTHON) $(FORMAT_STRING_SCRIPT) $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) \
		<(cut -f1 $(KINACT_DATA) | sort | uniq) 6 --all >$@.tmp
	cat <(sed -n '1p' $@.tmp) <(sed '1d' $@.tmp | sort -k 1d,1 -k 2d,2) >$@
	rm $@.tmp

$(HUMAN_STRING_COEXP2): $(HUMAN_STRING_COEXP2_FULL) $(KINACT_DATA)
	$(ASSOCNET) --header-in --method=spearman $< | sed '1s/assoc/string.coexp.cor/' >$@.tmp
	$(RSCRIPT) $(FILTER_BY_PROTLIST_SCRIPT) $@.tmp \
		<(cut -f1 $(KINACT_DATA) | sort | uniq) $@
	rm $@.tmp

$(HUMAN_STRING_COEXP_KINOME): $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) \
		$(KINACT_DATA) $(HUMAN_KINOME) $(FORMAT_STRING_SCRIPT)
	$(PYTHON) $(FORMAT_STRING_SCRIPT) $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) \
		$(HUMAN_KINOME) 6 --all >$@.tmp
	cat <(sed -n '1p' $@.tmp) <(sed '1d' $@.tmp | sort -k 1d,1 -k 2d,2) >$@
	rm $@.tmp

$(HUMAN_STRING_COEXP2_KINOME): $(HUMAN_STRING_COEXP2_FULL) $(HUMAN_KINOME)
	$(ASSOCNET) --header-in --method=spearman $< | sed '1s/assoc/string.coexp.cor/' >$@.tmp
	$(RSCRIPT) $(FILTER_BY_PROTLIST_SCRIPT) $@.tmp $(HUMAN_KINOME) $@
	# rm $@.tmp

# STRING cooccurrence
$(HUMAN_STRING_COOCC): $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) $(KINACT_DATA) \
		$(FORMAT_STRING_SCRIPT)
	$(PYTHON) $(FORMAT_STRING_SCRIPT) $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) \
		<(cut -f1 $(KINACT_DATA) | sort | uniq) 5 >$@.tmp
	cat <(sed -n '1p' $@.tmp) <(sed '1d' $@.tmp | sort -k 1d,1 -k 2d,2) >$@
	rm $@.tmp

$(HUMAN_STRING_COOCC2_FULL): $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) $(KINACT_DATA)
	$(PYTHON) $(FORMAT_STRING_SCRIPT) $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) \
		<(cut -f1 $(KINACT_DATA) | sort | uniq) 5 --all >$@.tmp
	cat <(sed -n '1p' $@.tmp) <(sed '1d' $@.tmp | sort -k 1d,1 -k 2d,2) >$@
	rm $@.tmp

$(HUMAN_STRING_COOCC2): $(HUMAN_STRING_COOCC2_FULL) $(KINACT_DATA)
	$(ASSOCNET) --header-in --method=spearman $< | sed '1s/assoc/string.coocc.cor/' >$@.tmp
	$(RSCRIPT) $(FILTER_BY_PROTLIST_SCRIPT) $@.tmp \
		<(cut -f1 $(KINACT_DATA) | sort | uniq) $@
	# rm $@.tmp

$(HUMAN_STRING_COOCC_KINOME): $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) \
		$(KINACT_DATA) $(HUMAN_KINOME) $(FORMAT_STRING_SCRIPT)
	$(PYTHON) $(FORMAT_STRING_SCRIPT) $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) \
		$(HUMAN_KINOME) 5 --all >$@.tmp
	cat <(sed -n '1p' $@.tmp) <(sed '1d' $@.tmp | sort -k 1d,1 -k 2d,2) >$@
	rm $@.tmp

$(HUMAN_STRING_COOCC2_KINOME): $(HUMAN_STRING_COOCC2_FULL) $(HUMAN_KINOME)
	$(ASSOCNET) --header-in --method=spearman $< | sed '1s/assoc/string.coexp.cor/' >$@.tmp
	$(RSCRIPT) $(FILTER_BY_PROTLIST_SCRIPT) $@.tmp $(HUMAN_KINOME) $@
	# rm $@.tmp

# STRING experimental
$(HUMAN_STRING_EXPER): $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) $(KINACT_DATA) \
		$(FORMAT_STRING_SCRIPT)
	$(PYTHON) $(FORMAT_STRING_SCRIPT) $(HUMAN_STRING_RAW) $(STRING_ID_MAPPING) \
		<(cut -f1 $(KINACT_DATA) | sort | uniq) 7 >$@.tmp
	cat <(sed -n '1p' $@.tmp) <(sed '1d' $@.tmp | sort -k 1d,1 -k 2d,2) >$@
	rm $@.tmp

# Kegg pathway reference
$(KEGG_PATH_REFERENCE): $(KEGG_RELATIONSHIPS) $(UNIPROT_ID_MAPPING) \
		$(KEGG_PATH_REF_SCRIPT)
	$(PYTHON) $(KEGG_PATH_REF_SCRIPT) | sort | uniq >$@

# Combined protein group set
$(COMBINED_GROUPING): $(GO_CELL_LOCATION) $(KEGG_PATH_REFERENCE)
	cat $^ >$@

# Kegg-based validation sets
$(DATADIR)/validation-set-kegg-%.tsv: $(DATADIR)/kegg-%-rels.tsv \
			$(UNIPROT_ID_MAPPING) $(KEGG_VAL_SCRIPT) $(KSEA_DATA)
	$(RSCRIPT) $(KEGG_VAL_SCRIPT) $(wordlist 1,2,$^) $@

# Negative validation set
$(NEG_VALSET): $(NEG_VAL_SCRIPT) $(KINACT_ASSOC) $(VAL_SET) $(PROTEIN_GROUPING)
	$(RSCRIPT) $^ $@

# Calculate human proteome amino-acid frequencies
$(AA_FREQS): $(AA_FREQS_SCRIPT)
	$(PYTHON) $(AA_FREQS_SCRIPT) >$@

# Kinase beads activities https://doi.org/10.1101/158295
$(KINASE_BEADS): $(KINASE_BEADS_RAW)
	awk -F"\t" 'BEGIN{OFS="\t"}{for (i=1; i<=NF; i++){if ($$i==""){$$i="NA"}}; print}' $< \
		| cut -f2,`seq 4 3 61 | tr '\n' ',' | sed 's/,$$//'` \
		| sed '1{s/ log2 Fold-Change//g;s/ /./g;s/\(.*\)/\L\1/g}' >$@.tmp
	$(RSCRIPT) src/melt.r $@.tmp $@.tmp2
	rm $@.tmp
	sort -k1 $@.tmp2 >$@
	rm $@.tmp2
	$(RSCRIPT) src/impute.r $@ $(DATADIR)/kinase-beads-activities-imp.tsv

# Hotspots
$(HOTSPOTS): $(HOTSPOTS_RAW) $(TYR_HOTSPOTS_RAW) $(ENSEMBL_ID_MAPPING)
	cat $(wordlist 1,2,$^) | sed 's/, /\t/g' | cut -f2,4,6,7 | \
		awk -F"\t" -vOFS="\t" '{if ($$4==0){rel=1}else{rel=-log($$4)/log(10)/100}; print $$1, $$2, $$3, rel}' | \
		sort -k1 | join -t'	' $(ENSEMBL_ID_MAPPING) - | cut -f2-5 >$@

###########################
## Association score tables

# Pearson's correlation
$(OUTDIR)/%-pcor.tsv: $(DATADIR)/%-imp.tsv $(OUTDIR)
	$(ASSOCNET) --method=pearson $< >$@

# Spearman's correlation
$(OUTDIR)/%-scor.tsv: $(DATADIR)/%-imp.tsv $(OUTDIR)
	$(ASSOCNET) --method=spearman $< >$@

# Correlation of correlation
$(OUTDIR)/%-pcor2.tsv: $(OUTDIR)/%-pcor.tsv $(OUTDIR)
	$(ASSOCNET) --header-in --method=pearson $< >$@

$(OUTDIR)/%-scor2.tsv: $(OUTDIR)/%-scor.tsv $(OUTDIR)
	$(ASSOCNET) --header-in --method=spearman $< >$@

# Normalized FunChisq
$(OUTDIR)/%-nfchisq.tsv: $(DATADIR)/%-discr.tsv $(NFCHISQ_SCRIPT) $(OUTDIR)
	$(RSCRIPT) $(NFCHISQ_SCRIPT) $< $@
# Here I create three different correlation tables with the same script assoc_methods.r
# Funchisq, partcor and FNN_mutinfo
$(OUTDIR)/%-all.tsv: $(DATADIR)/%-imp.tsv $(ASSOC_SCRIPT) $(OUTDIR)
	$(RSCRIPT) $(ASSOC_SCRIPT) $< $@ all

# Partial correlation
$(OUTDIR)/%-partcor.tsv: $(DATADIR)/%-imp.tsv $(ASSOC_SCRIPT) $(OUTDIR)
	$(RSCRIPT) $(ASSOC_SCRIPT) $< $@ partcor

# Mutual information
$(OUTDIR)/%-mut_info.tsv: $(DATADIR)/%-discr.tsv $(ASSOC_SCRIPT) $(OUTDIR)
	$(RSCRIPT) $(ASSOC_SCRIPT) $< $@ mut_info

# Fast-nearest-neighbors (FNN) mutual information
$(OUTDIR)/%-fnn_mut_info.tsv: $(DATADIR)/%-imp.tsv $(ASSOC_SCRIPT) $(OUTDIR)
	$(RSCRIPT) $(ASSOC_SCRIPT) $< $@ fnn_mut_info

# Linear model based predictions (poor man's MRA)
$(OUTDIR)/%-fvalue.tsv: $(DATADIR)/%-imp.tsv $(LM_PRED_SCRIPT) $(KIN_INVIVO_CONDS) $(OUTDIR)
	$(RSCRIPT) $(LM_PRED_SCRIPT) $< $@

# Pairwise correlations, taking into account perturbations
$(OUTDIR)/%-paircor.tsv: $(DATADIR)/%-imp.tsv $(PAIR_COR_SCRIPT) $(KIN_INVIVO_CONDS) $(OUTDIR)
	$(RSCRIPT) $(PAIR_COR_SCRIPT) $< $@

# Filter out indirect associations
$(OUTDIR)/%-filter.tsv: $(OUTDIR)/%.tsv
	$(ASSOCNET_FILTER) --header-in $< >$@

# Regulatory site-based activity predictions
$(REG_SITE_ASSOC): $(REG_SITES) $(REG_SITE_COR_SCRIPT)
	$(RSCRIPT) $(REG_SITE_COR_SCRIPT)

###############################
## Kinase-substrate predictions

# Calculate kinase-substrate PSSM scores
$(KIN_KIN_SCORES): $(AA_FREQS) $(PSSM_SCRIPT) $(KIN_SUBSTR_TABLE)	\
					$(HUMAN_KINASE_TABLE) $(PHOSPHOSITES)			\
					$(KINACT_DATA) $(PSSM_LIB) $(HOTSPOTS)
	$(PYTHON) $(PSSM_SCRIPT) <(cut -f1 $(KINACT_DATA) | sort | uniq) $@

$(KIN_KIN_SCORES_KINOME): $(AA_FREQS) $(PSSM_SCRIPT) $(KIN_SUBSTR_TABLE)	\
					$(HUMAN_KINASE_TABLE) $(PHOSPHOSITES)			\
					$(HUMAN_KINOME) $(PSSM_LIB) $(HOTSPOTS)
	$(PYTHON) $(PSSM_SCRIPT) $(HUMAN_KINOME) $@

$(KIN_SCORE_DIST): $(AA_FREQS) $(PSSM_DIST_SCRIPT) $(KIN_SUBSTR_TABLE) \
					$(HUMAN_KINASE_TABLE) $(KINACT_DATA) $(PSSM_LIB)
	$(PYTHON) $(PSSM_DIST_SCRIPT) $(KINACT_DATA)

###################
## Merged predictor

$(MERGED_PRED): $(MERGED_PRED_SOURCES) $(MERGE_SCRIPT)
	$(RSCRIPT) $(MERGE_SCRIPT) $@ $(filter-out $(MERGE_SCRIPT),$^)

$(MERGED_PRED_KINOME): $(MERGED_PRED_SOURCES_KINOME) $(MERGE_SCRIPT)
	$(RSCRIPT) $(MERGE_SCRIPT) $@ $(filter-out $(MERGE_SCRIPT),$^)

#############
## Validation

$(IMGDIR)/%-val.pdf: $(OUTDIR)/%.tsv $(VAL_SET) $(NEG_VALSET) $(IMGDIR) $(VAL_SCRIPT)
	$(RSCRIPT) $(VAL_SCRIPT) $(wordlist 1,3,$^)

$(IMGDIR)/validation.pdf:
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@ $(IMGDIR)/*-val.pdf
