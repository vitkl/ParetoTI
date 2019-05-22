##' Convert gene identifiers with AnnotationHub
##' @rdname map_gene_annot
##' @name map_gene_annot
##' @description map_gene_annot() Maps gene annotations using OrgDb databases from AnnotationHub. Most useful for mapping gene identifiers. This function can pick one mapped annotation for each key.
##' @return data.frame containing mapped identifiers
##' @param taxonomy_id Taxonomy id of your species. 9606 is human, 10090 is mouse. Find your species id using taxonomy search by UniProt: https://www.uniprot.org/taxonomy/.
##' @param keys gene or protein identifiers of \code{keytype}
##' @param columns columns to retrieve. Details: \link[AnnotationDbi]{colsAndKeytypes}.
##' @param keytype type of \code{keys} identifiers. Use \code{keytypes(map_gene_annot(taxonomy_id = 9606, return_org_db = T))} to find which keytypes are available. Most frequent options: "ENSEMBL", "SYMBOL", "ALIAS".
##' @param pick_one pick one mapped annotation (from columns) for each key. Works only with one column.
##' @param localHub set localHub = FALSE for working offline. Details \link[AnnotationHub]{AnnotationHub}.
##' @param return_org_db return OrgDb database link instead of gene annotations.
##' @param ann_hub_cache If needed, set location of local AnnotationHub cache to specific directory by providing it here. To view where cache is currently located use getAnnotationHubOption("CACHE")
##' @export map_gene_annot
##' @import AnnotationHub
map_gene_annot = function(taxonomy_id = 9606, keys = "TP53", columns = c("ENSEMBL"),
                          keytype = "SYMBOL", pick_one = TRUE,
                          localHub = FALSE, return_org_db = FALSE,
                          ann_hub_cache = AnnotationHub::getAnnotationHubOption("CACHE")) {

  # find annotation data base for taxonomy_id
  setAnnotationHubOption("CACHE", ann_hub_cache)
  hub = AnnotationHub(localHub = localHub)

  record = subset(hub, taxonomyid == taxonomy_id & rdataclass ==
                    "OrgDb")
  org_db = hub[[names(record)]]

  if(isTRUE(return_org_db)) return(org_db)

  annot = as.data.table(AnnotationDbi::select(org_db, keys = keys,
                                              columns = c(columns, keytype),
                                              keytype = keytype))

  if(pick_one){
    if(length(columns) != 1) stop("pick_one works only with one column")

    # pick one mapped annotation for each key
    annot[, c(columns) := get(columns)[1], by = .(get(keytype))]
    annot = unique(annot)

    # create readable row names
    annot[is.na(get(columns)), ROWNAMES := get(keytype)]
    annot[!is.na(get(columns)), ROWNAMES := get(columns)]

  }

  annot = as.data.frame(annot)

  if(!(sum(duplicated(annot[, keytype])) > 0)) rownames(annot) = annot[, keytype]
  annot
}
