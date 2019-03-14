##' Annotate genes with GO terms and measure their activities
##' @rdname measure_activity
##' @name map_go_annot
##' @description map_go_annot() annotates genes with Gene Ontology terms including less specific terms. It automatically loads correct OrgDb package given \code{taxonomy_id}. Then it selects GOALL annotations for genes provided in \code{keys}, filter by ontology branch and evidence codes. Finally it maps GO term names using GO.db and returns annotations in list and data.table formats. This function can be used to retrieve other gene annotations from OrgDb (e.g OMIM, PFAM) but those will not be mapped to readable names.
##' @details GO consortium annotates genes with the most specific terms so there is a need to propagate annotations to less specific parent terms. For example, "apoptotic process in response to mitochondrial fragmentation" (GO:0140208) is an "apoptotic process" (GO:0006915). Conceptually, if gene functions in the first it also functions in the second but it is only directly annotated by the first term. So, using GO annotations for functional enrichment requires propagating annotations to less specific terms.
##' @return list (S3 class ParetoTI_annot) containing GO annotations in list and data.table formats, and OrgDb database link.
##' @param taxonomy_id Taxonomy id of your species. 9606 is human, 10090 is mouse. Find your species id using taxonomy search by UniProt: https://www.uniprot.org/taxonomy/.
##' @param keys gene or protein identifiers of \code{keytype}
##' @param columns columns to retrieve. Default is GOALL: GO Identifiers (includes less specific terms), \link[AnnotationDbi]{colsAndKeytypes}.
##' @param keytype type of \code{keys} identifiers. Use \code{keytypes(map_go_annot(taxonomy_id = 9606, return_org_db = T))} to find which keytypes are available.
##' @param ontology_type specify which branch of gene ontology to retrieve annotations from.
##' @param evidence_code specify which evidence should have been used to annotate genes with GO terms. Use \code{keys(map_go_annot(..., return_org_db = T), "EVIDENCE")} to find which codes exist. See http://www.geneontology.org/page/guide-go-evidence-codes for explanations and details.
##' @param localHub set localHub = FALSE for working offline. Details \link[AnnotationHub]{AnnotationHub}.
##' @param return_org_db return OrgDb database link instead of gene annotations.
##' @param ann_hub_cache If needed, set location of local AnnotationHub cache to specific directory by providing it here. To view where cache is currently located use getAnnotationHubOption("CACHE")
##' @export map_go_annot
##' @import AnnotationHub
map_go_annot = function(taxonomy_id = 9606, keys = c("TP53", "ALB"),
                        columns = c("GOALL"), keytype = "ALIAS",
                        ontology_type = c("BP", "MF", "CC"),
                        evidence_code = "all",
                        localHub = FALSE, return_org_db = FALSE,
                        ann_hub_cache = AnnotationHub::getAnnotationHubOption("CACHE")) {
  # find annotation data base for taxonomy_id
  setAnnotationHubOption("CACHE", ann_hub_cache)
  hub = AnnotationHub(localHub = localHub)
  record = subset(hub, taxonomyid == taxonomy_id & rdataclass == "OrgDb")
  if(length(record) == 0) stop(paste0("No OrgDb for taxonomy_id: ", taxonomy_id, " in the AnnotationHub(). Check that you entered correct taxonomy_id. \nElse you may need to create your own annotations: named list where names are annotation terms and each element is a vector listing genes/proteins."))
  org_db = hub[[names(record)]]
  if(isTRUE(return_org_db)) return(org_db)
  suppressMessages({
    annot = as.data.table(AnnotationDbi::select(org_db, keys = keys,
                                                columns = c(keytype, columns),
                                                keytype = keytype))
  })

  # filter by ontology branch/type:
  ont_keys = c("BP", "MF", "CC")
  if(mean(ontology_type %in% ont_keys) == 1) {
    annot = annot[ONTOLOGYALL %in% ontology_type]
  } else stop("ontology_type should be one of ", paste0(ont_keys, collapse = ", "))
  # filter by evidence type and remove EVIDENCEALL column:
  evid_keys = keys(org_db, "EVIDENCE")
  if(mean(evidence_code %in% evid_keys) == 1) {
    annot = annot[EVIDENCEALL %in% evidence_code]
  } else if(evidence_code == "all") NULL else {
    stop("evidence_code should be one of: \"all\", ", paste0(evid_keys, collapse = ", "))
  }
  annot[, EVIDENCEALL := NULL]
  annot = unique(annot)
  # add GO term names
  if(sum(c("GO", "GOALL") %in% columns[1]) >= 1){
    suppressMessages({
      annot_names = as.data.table(AnnotationDbi::select(GO.db::GO.db,
                                                        keys = annot[, unique(GOALL)],
                                                        keytype = "GOID",
                                                        columns = c("GOID", "TERM"),
                                                        verbose = F))
    })
    annot = merge(annot, annot_names, by.x = columns[1], by.y = "GOID", all = TRUE)
  }
  # generate list format and return:
  structure(list(annot_list = split(annot[, get(keytype)],
                                    annot[, get(columns[1])]),
                 annot_dt = annot, org_db = org_db),
            class = "ParetoTI_annot")
}



##' @rdname measure_activity
##' @name filter_go_annot
##' @description filter_go_annot() filters GO annotations by ontology file, GO branch, number of genes annotated
##' @param annot GO annotations, output of map_go_annot()
##' @param ontology_file path (url) to ontology file in obo format
##' @param lower lower limit on size of gene sets (inclusive)
##' @param upper upper limit on size of gene sets (exclusive)
##' @export filter_go_annot
##' @import data.table
##' @examples
##' annot = map_go_annot(taxonomy_id = 9606,
##'                      keys = c("TP53", "ALB", "GPX1"),
##'                      columns = c("GOALL"), keytype = "ALIAS",
##'                      ontology_type = c("BP", "MF", "CC"))
##' annot2 = filter_go_annot(annot, ontology_file = "http://www.geneontology.org/ontology/subsets/goslim_generic.obo",
##'            lower = 50, upper = 2000, ontology_type = "BP")
filter_go_annot = function(annot, ontology_file = NULL, lower = 1, upper = Inf,
                           ontology_type = c("BP", "MF", "CC")){
  if(!is(annot, "ParetoTI_annot")) {
    if(sum(colnames(annot$annot_dt) %in% c("GO", "GOALL")) == 0) stop("filter_go_annot(): annot should be the output of map_go_annot() containing GO annotations (GO or GOALL)")
  }
  colname = colnames(annot$annot_dt)
  go_colname = colname[colname %in% c("GO", "GOALL")]
  ont_colname = colname[colname %in% c("ONTOLOGY", "ONTOLOGYALL")]
  # load ids from ontology
  if(!is.null(ontology_file)){
    ont_ids = ontologyIndex::get_ontology(ontology_file,
                                          propagate_relationships = "is_a",
                                          extract_tags = "minimal")$id
    ont_ids = ont_ids[grepl("GO", ont_ids)]
    # filter by ontology
    annot$annot_dt = annot$annot_dt[get(go_colname) %in% ont_ids]
  }
  # filter by number of annotated genes
  annot$annot_dt[, N := .N, by = go_colname]
  annot$annot_dt = annot$annot_dt[N >= lower & N < upper]
  # filter by ontology type/branch
  ont_keys = c("BP", "MF", "CC")
  if(mean(ontology_type %in% ont_keys) == 1) {
    annot$annot_dt = annot$annot_dt[get(ont_colname) %in% ontology_type]
  }
  # filter annotations in list format
  annot$annot_list = annot$annot_list[unique(annot$annot_dt[,get(go_colname)])]
  annot
}

##' @rdname measure_activity
##' @name measure_activity
##' @description measure_activity() annotates genes with GO terms and measures their activities
##' @param expr_mat expression matrix (genes in rows, cells in columns) or one of: dgCMatrix, ExpressionSet, and SummarizedExperiment or SingleCellExperiment both of which require assay_name.
##' @param which which set activities to measure? Currently implemented "BP", "MF", "CC" gene ontology subsets. Use "gwas" for constructing gene sets with gwas_catalog v1.0.2. GWAS option works only with human hgnc identifiers. Use "promoter_TF" to look at TF activities measured by mapping TF->target gene regulation via co-occurence of multiple motifs for the same TF in promoters (courtesy of Krzysztof Polanski who generated this resource). TF data is added supplied with the package and works only for human data, only SYMBOL and ENSEMBL can be used as keytype. Use "custom" to provide your own gene sets in \code{annot_dt} table.
##' @param activity_method find activities using \link[ParetoTI]{find_set_activity_AUCell} or \link[ParetoTI]{find_set_activity_pseudoinv}?
##' @param return_as_matrix return matrix (TRUE) or data.table (FALSE)
##' @param annot_dt data.table, data.frame or matrix containing custom gene set annotations. The 1st column should contain gene set identifier or name, the 2nd column should contain gene identifiers matching \code{keys}.
##' @param assay_name name of assay in SummarizedExperiment or SingleCellExperiment, normally counts or logcounts
##' @param aucell_options additional options specific to AUCell, details: \link[ParetoTI]{find_set_activity_AUCell}
##' @export measure_activity
##' @import data.table
##' @import AnnotationHub
##' @examples
##' activ = measure_activity(expr_mat, which = "BP",
##'                          taxonomy_id = 9606, keytype = "ALIAS",
##'                          ontology_file = load_go_ontology("./data/",
##'                                                   "goslim_generic.obo"))
measure_activity = function(expr_mat, which = c("BP", "MF", "CC", "gwas", "promoter_TF", "custom")[1],
                            activity_method = c("AUCell", "pseudoinverse")[1],
                            keys = rownames(expr_mat), keytype = "ALIAS",
                            ontology_file = NULL, taxonomy_id = 10090,
                            evidence_code = "all", localHub = FALSE,
                            ann_hub_cache = AnnotationHub::getAnnotationHubOption("CACHE"),
                            lower = 1, upper = Inf,
                            return_as_matrix = FALSE,
                            annot_dt = NULL,
                            assay_name = "logcounts",
                            aucell_options = list(aucMaxRank = nrow(expr_mat) * 0.05,
                                                  binary = F, nCores = 3,
                                                  plotStats = TRUE)) {
  if(!which %in% c("BP", "MF", "CC", "gwas", "promoter_TF", "custom")) stop("which should be one of c(\"BP\", \"MF\", \"CC\", \"gwas\", \"promoter_TF\", \"custom\")")

  # Retrieve annotations ---------------------------------------------------------
  if(mean(which %in% c("BP", "MF", "CC")) == 1){
    # Map GO annotations using AnnotationHub -------------------------------------
    annot = map_go_annot(taxonomy_id = taxonomy_id, keys = keys,
                         columns = c("GOALL"), keytype = keytype,
                         ontology_type = which,
                         evidence_code = evidence_code,
                         localHub = localHub, return_org_db = FALSE,
                         ann_hub_cache = ann_hub_cache)
    # Choose terms from GO slim subset
    annot = filter_go_annot(annot, ontology_file = ontology_file,
                            lower = lower, upper = upper,
                            ontology_type = which)
    set_id_col = "GOALL"
    set_name_col = "TERM"

  } else if(which == "gwas") {

    # Map GWAS phenotype and disease associations GWAS Catalog -------------------
    annot = map_gwas_annot(taxonomy_id = taxonomy_id, keys = keys,
                           keytype = keytype, localHub = localHub,
                           ann_hub_cache = ann_hub_cache,
                           lower = lower, upper = upper)
    set_id_col = "MAPPED_TRAIT_ID"
    set_name_col = "MAPPED_TRAIT_NAME"

  } else if(which == "promoter_TF") {

    # Load TF-target map  --------------------------------------------------------
    if(!keytype %in% c("SYMBOL", "ENSEMBL")) stop("When which == \"promoter_TF\", the keytype should be SYMBOL or ENSEMBL")
    annot = map_promoter_tf(keys = keys, keytype = keytype,
                            lower = lower, upper = upper)
    set_id_col = "TF"
    set_name_col = set_id_col
    keytype = "TARGET"

  } else if(which == "custom") {

    # Use custom annotations  ----------------------------------------------------

    annot = list(annot_dt = as.data.frame(annot_dt))
    set_id_col = colnames(annot_dt)[1]
    set_name_col = colnames(annot_dt)[1]
    keytype = colnames(annot_dt)[2]

  }

  # rename rows of expression matrix to match keys
  if(length(keys) == nrow(expr_mat)) {
    # remove genes with NA key
    expr_mat = expr_mat[!is.na(keys), ]
    rownames(expr_mat) = keys[!is.na(keys)]
  }

  # Find "activity" of each functional gene group  -------------------------------
  if(activity_method == "AUCell"){

    # using AUCell  --------------------------------------------------------------
    activ = find_set_activity_AUCell(expr_mat, assay_name = assay_name,
                                     aucMaxRank = aucell_options$aucMaxRank,
                                     gene_sets = annot$annot_dt,
                                     gene_col = keytype,
                                     set_id_col = set_id_col,
                                     set_name_col = set_name_col,
                                     binary = aucell_options$binary,
                                     nCores = aucell_options$nCores,
                                     plotHist = FALSE,
                                     plotStats = aucell_options$plotStats)

  } else if (activity_method == "pseudoinverse") {

    # or pseudoinverse(annotation_matrix) * expression ---------------------------
    activ = find_set_activity_pseudoinv(expr_mat,
                                        assay_name = assay_name,
                                        gene_sets = annot$annot_dt,
                                        gene_col = keytype,
                                        set_id_col = set_id_col,
                                        set_name_col = set_name_col,
                                        noself = FALSE)

  }
  # when naming matrix dimensions ad
  setnames(activ, colnames(activ), gsub(" ", "_", colnames(activ)))
  setnames(activ, colnames(activ), gsub("-", "_", colnames(activ)))
  setnames(activ, colnames(activ), gsub(",", "_", colnames(activ)))
  setnames(activ, colnames(activ), gsub("/", "_", colnames(activ)))
  setnames(activ, colnames(activ), gsub("'", "_", colnames(activ)))
  setnames(activ, colnames(activ), gsub("\\[", "_", colnames(activ)))
  setnames(activ, colnames(activ), gsub("\\]", "_", colnames(activ)))
  setnames(activ, colnames(activ), gsub("\\(", "_", colnames(activ)))
  setnames(activ, colnames(activ), gsub("\\)", "_", colnames(activ)))
  if(return_as_matrix) activ = as.matrix(activ, rownames = "cells")
  activ
}

##' @rdname measure_activity
##' @name load_go_ontology
##' @description load_go_ontology() load GO slim ontology file and return file path
##' @param ont_dir directory where to save ontology
##' @param ont_file which ontology file to load?
##' @export load_go_ontology
load_go_ontology = function(ont_dir = "./data/", ont_file = "goslim_generic.obo"){
  file = paste0(ont_dir, ont_file)
  if(!file.exists(file)){
    download.file(paste0("http://www.geneontology.org/ontology/subsets/", ont_file),
                  file)
  }
  file
}

##' @rdname measure_activity
##' @name map_gwas_annot
##' @description map_gwas_annot() Load known GWAS loci for all phenotypes in GWAS Catalog. For each the \code{keys} find which phenotypes it is associated with using genes mapped by GWAS Catalog.
##' @param gwas_url where to download annotations? Useful for other versions of GWAS Catalog.
##' @param gwas_file name of the file where to store GWAS Catalog data locally (in ann_hub_cache directory).
##' @export map_gwas_annot
map_gwas_annot = function(taxonomy_id = 9606, keys = c("TP53", "ZZZ3"),
                          keytype = "SYMBOL", localHub = FALSE,
                          ann_hub_cache = AnnotationHub::getAnnotationHubOption("CACHE"),
                          gwas_url = "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
                          gwas_file = "gwas_catalog_v1.0.2-associations_e93_r2019-01-31.tsv",
                          lower = 1, upper = Inf) {
  if(taxonomy_id %in% c(9606)){
    # find annotation data base for taxonomy_id
    setAnnotationHubOption("CACHE", ann_hub_cache)
    hub = AnnotationHub(localHub = localHub)
    record = subset(hub, taxonomyid == taxonomy_id & rdataclass == "OrgDb")
    org_db = hub[[names(record)]]
    suppressMessages({
      annot = as.data.table(AnnotationDbi::select(org_db, keys = keys,
                                                  columns = unique(c(keytype, "SYMBOL")),
                                                  keytype = keytype))
    })
  } else stop("GWAS are annotations of human genes only (taxonomy_id = 9606)")
  gwas_file = paste0(ann_hub_cache, "/", gwas_file)
  if(!file.exists(gwas_file)){
    if(!dir.exists(ann_hub_cache)) dir.create(ann_hub_cache, recursive = TRUE)
    download.file(gwas_url,
                  destfile = gwas_file)
  }
  gwas_data = fread(gwas_file, stringsAsFactors = FALSE, header = TRUE)
  # expand mapped traits and gene associations to one trait-one gene per row
  # remove traits not mapped to ontology
  trait_id = unique(gwas_data[MAPPED_TRAIT_URI != "",
                              .(MAPPED_TRAIT_URI, MAPPED_TRAIT)])
  # produce clean trait_id
  trait_id = trait_id[, .(MAPPED_TRAIT_ID = unlist(strsplit(MAPPED_TRAIT_URI,", ")),
                          MAPPED_TRAIT_NAME = ifelse(length(unlist(strsplit(MAPPED_TRAIT,", "))) >
                                                       length(unlist(strsplit(MAPPED_TRAIT_URI,", "))),
                                                     MAPPED_TRAIT,
                                                     unlist(strsplit(MAPPED_TRAIT,", "))),
                          MAPPED_TRAIT = MAPPED_TRAIT),
                      by = MAPPED_TRAIT_URI]
  trait_id[, MAPPED_TRAIT_ID := gsub("^.+/", "",MAPPED_TRAIT_ID)]
  trait_id[, MAPPED_TRAIT_NAME := paste0(MAPPED_TRAIT_NAME, "_", MAPPED_TRAIT_ID)]
  # produce clean gene_names
  gene_names = unique(gwas_data[MAPPED_GENE != "", ][, .(MAPPED_GENE_NAME =
                                                           unlist(strsplit(MAPPED_GENE," - "))),
                                                     by = MAPPED_GENE])

  # merge MAPPED_GENE_NAMEs and MAPPED_TRAIT_IDs to all data
  gwas_data = merge(gwas_data, trait_id,
                    by = c("MAPPED_TRAIT_URI", "MAPPED_TRAIT"),
                    all.y = TRUE, as.x = FALSE,
                    allow.cartesian = TRUE)
  gwas_data = merge(gwas_data, gene_names,
                    by = c("MAPPED_GENE"),
                    all.y = TRUE, as.x = FALSE,
                    allow.cartesian = TRUE)
  # merge gene identifiers in the dataset to gene names in GWAS catalog
  gwas_data = merge(gwas_data, annot,
                    by.x = c("MAPPED_GENE_NAME"), by.y = c("SYMBOL"),
                    all.y = FALSE, as.x = FALSE,
                    allow.cartesian = TRUE)
  # filter by lower and upper
  gwas_data[, N_genes := uniqueN(MAPPED_GENE_NAME), by = .(MAPPED_TRAIT_ID)]
  gwas_data = gwas_data[N_genes >= lower & N_genes <= upper]

  list(annot_dt = unique(gwas_data))
}


##' @rdname measure_activity
##' @name map_promoter_tf
##' @description map_promoter_tf() load GO slim ontology file and return file path
##' @export map_promoter_tf
map_promoter_tf = function(keys = c("TP53", "ZZZ3"),
                           keytype = "SYMBOL",
                           lower = 1, upper = Inf) {

  data("promoter_TF_motifs", package = "ParetoTI", envir = environment())
  promoter_TF_motifs = unique(promoter_TF_motifs[, c("TF", keytype)])
  # find number of targets per TF
  n_targ = table(promoter_TF_motifs$TF)
  # keep TFs between upper and lower boundary
  tf_size_filter = names(n_targ)[n_targ >= lower & n_targ <= upper]
  promoter_TF_motifs = promoter_TF_motifs[promoter_TF_motifs[, keytype] %in% keys &
                                            promoter_TF_motifs$TF %in% tf_size_filter,]
  colnames(promoter_TF_motifs) = c("TF", "TARGET")
  # check that only letters, numbers and underscores are present - remove dots
  promoter_TF_motifs$TF = gsub("\\.", "_", promoter_TF_motifs$TF)

  annot = list(annot_dt = promoter_TF_motifs)
}
