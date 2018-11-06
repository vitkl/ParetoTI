##' Use OrgDb packages to annotate genes with GO terms (including less specific)
##' @rdname map_go_annot
##' @name map_go_annot
##' @description map_go_annot() annotates genes with Gene Ontology terms including less specific terms. It automatically loads correct OrgDb package given \code{taxonomy_id}. Then it selects GOALL annotations for genes provided in \code{keys}, filter by ontology branch and evidence codes. Finally it maps GO term names using GO.db and returns annotations in list and data.table formats. This function can be used to retrieve other gene annotations from OrgDb (e.g OMIM, PFAM) but those will not be mapped to readable names.
##' @details GO consortium annotates genes with the most specific terms so there is a need to propagate annotations to less specific parent terms. For example, "apoptotic process in response to mitochondrial fragmentation" (GO:0140208) is an "apoptotic process" (GO:0006915). Conceptually, if gene functions in the first it also functions in the second but it is only directly annotated by the first term. So, using GO annotations for functional enrichment requires propagating annotations to less specific terms.
##' @return list (S3 class ParetoTI_annot) containing GO annotations in list and data.table formats, and OrgDb database link.
##' @param taxonomy_id Taxonomy id of your species. 9606 is human, 10090 is mouse. Find your species id using taxonomy search by UniProt: https://www.uniprot.org/taxonomy/.
##' @param keys gene or protein identifiers of \code{keytype}
##' @param columns columns to retrieve. Default is GOALL: GO Identifiers (includes less specific terms), \link[AnnotationDbi]{colsAndKeytypes}.
##' @param keytype type of \code{keys} identifiers. Use \code{keytypes(map_go_annot(..., return_org_db = T))} to find which keytypes are available.
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



##' @rdname map_go_annot
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
