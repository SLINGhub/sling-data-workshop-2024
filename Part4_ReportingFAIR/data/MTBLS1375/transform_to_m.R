library(tidyverse)
skyl_col = cols(
  OriginalTransitionName = col_character(),
  AnalyticalID = col_character(),
  PrecursorMz = col_double(),
  ProductMz = col_double(),
  `Fragment Ion` = col_character(),
  RT = col_double(),
  Area = col_double(),
  `Product Ion Formula` = col_logical(),
  `Precursor Note` = col_character(),
  `Peptide Note` = col_character(),
  `Transition Replicate Note` = col_logical()
)
skyl <-
  read_csv("human_plasma_export.csv",
           na = c("NA", "#N/A"),
           col_types = skyl_col) %>%
  mutate("OriginalTransitionName" = str_replace(.$"OriginalTransitionName", " ([abc])$", "\\1")) %>%
  mutate(
    "OriginalTransitionName" = str_replace(
      .$"OriginalTransitionName",
      "^([A-Zderx0-9 :;\\-_\\/\\(\\)\\[\\]]+)([abc]?)$",
      "\\1"
    ),
    "isoid" = str_replace(
      .$"OriginalTransitionName",
      "^([A-Zderx0-9 :;\\-_\\/\\(\\)\\[\\]]+)([abc]?)$",
      "\\2"
    )
  ) %>%
  mutate("OriginalTransitionName" = str_replace(.$"OriginalTransitionName", "^(.*)[ ]?$", "\\1")) %>%
  select("OriginalTransitionName", "isoid", everything()) %>% group_by(OriginalTransitionName, AnalyticalID, PrecursorMz, ProductMz, `Fragment Ion`, RT) %>%
  summarise(Area = sum(Area))

bp <-
  read_csv(file = "LipidCreatorVal_QC-filtered_SAMPLES-NIST_uM_04022020v2.csv")
a <-
  read_tsv("a_MTBLS1375_LC-MS_positive_reverse-phase_metabolite_profiling.txt")

# join so that only baker panel values remain as entries
d <- a %>% right_join(bp, by = "Sample Name")

e <- d %>% select(5, 43:length(colnames(d)))

# translate names to be compatible with SwissLipids / LipidMaps
f <-
  e %>% pivot_longer(cols = 2:length(colnames(.)) , names_to = "metabolite_identification") %>%
  mutate("metabolite_identification" = str_replace(.$"metabolite_identification", "^(.*)[ ]?$", "\\1"))

fullLeftJoin <-
  f %>% left_join(
    skyl,
    by = c(
      "metabolite_identification" = "OriginalTransitionName",
      "Extract Name" = "AnalyticalID"
    )
  )

# join tables and normalize names
h <- fullLeftJoin %>%
  mutate("metabolite_identification" = str_replace(.$"metabolite_identification", " ([abc])$", "\\1")) %>%
  mutate(
    "metabolite_identification" = str_replace(
      .$"metabolite_identification",
      "^([A-Zderx0-9 :;\\-_\\/\\(\\)\\[\\]]+)([abc]?)$",
      "\\1"
    ),
    "isoid" = str_replace(
      .$"metabolite_identification",
      "^([A-Zderx0-9 :;\\-_\\/\\(\\)\\[\\]]+)([abc]?)$",
      "\\2"
    )
  ) %>%
  mutate("metabolite_identification" = str_replace(.$"metabolite_identification", "^(.*)[ ]?$", "\\1")) %>%
  select("metabolite_identification", "isoid", everything())

# sum up values
i <-
  h %>% select(
    "metabolite_identification",
    "isoid",
    "Extract Name",
    "PrecursorMz",
    "ProductMz",
    "RT",
    "Area"
  ) %>%
  group_by(
    !!sym("metabolite_identification"),
    !!sym("Extract Name"),
    !!sym("PrecursorMz"),
    !!sym("ProductMz"),
    !!sym("RT")
  ) %>%
  summarise_at(c("Area"), sum, na.rm = TRUE)

# crete abundance summary columns
jsummary <-
  i %>% ungroup() %>% group_by(
    !!sym("metabolite_identification"),
    !!sym("PrecursorMz"),
    !!sym("ProductMz"),
    !!sym("RT")
  ) %>%
  summarise(
    "smallmolecule_abundance_sub" = mean(Area, na.rm = TRUE),
    "smallmolecule_abundance_stdev_sub" = sd(Area, na.rm = TRUE),
    "smallmolecule_abundance_std_error_sub" = sd(Area, na.rm = TRUE) / sqrt(n()-sum(is.na(Area)))
  )

# pull values into columns for each Extract Name
j <-
  i %>% ungroup() %>% mutate("RunId" = as.numeric(str_replace(
    .$"Extract Name",
    "^([0-9]+)_.*$",
    "\\1"
  ))) %>%
  arrange(RunId, `metabolite_identification`) %>%
  group_by(
    !!sym("metabolite_identification"),
    !!sym("Extract Name"),
    !!sym("PrecursorMz"),
    !!sym("ProductMz"),
    !!sym("RT")
  ) %>% select(-RunId) %>%
  pivot_wider(names_from = "Extract Name", values_from = Area)
jp <-
  j %>% left_join(
    jsummary,
    by = c(
      "metabolite_identification" = "metabolite_identification",
      "PrecursorMz" = "PrecursorMz",
      "ProductMz" = "ProductMz",
      "RT" = "RT"
    )
  ) %>%
  select(
    "metabolite_identification",
    "PrecursorMz",
    "ProductMz",
    "RT",
    "smallmolecule_abundance_sub",
    "smallmolecule_abundance_stdev_sub",
    "smallmolecule_abundance_std_error_sub",
    everything()
  )

# read precursor charges and adducts / transitions
precursorChargesAndAdducts <-
  read_csv("human_plasma_transitions.csv") %>%
  select("Precursor Name", "Precursor m/z", "Precursor Charge", "PrecursorAdduct", "ProductName", "Product m/z") %>%
  mutate(
    "Precursor Name" = str_replace(.$"Precursor Name", "HexCer", "Hex1Cer"),
    "PrecursorAdduct" = str_replace(.$"PrecursorAdduct", "\\[M \\+ H\\]\\+", "[M+H]1+")
  ) %>%
  mutate(
    "Precursor Name" = str_replace(.$"Precursor Name", "^(.*)[ ]?$", "\\1"),
    "ProductName" = str_replace(.$"ProductName", "#NAME?", "")
  ) %>%
  unite(fragmentation,
             "Product m/z",
             "ProductName",
             sep = "|",
             remove = FALSE
  ) %>%
  group_by(!!sym("Precursor Name"),
           !!sym("Precursor m/z"),
           !!sym("Precursor Charge"),
           !!sym("PrecursorAdduct"),
           !!sym("ProductName"),
           !!sym("Product m/z")) %>%
  distinct()

# join for precursor and adduct information
k <-
  jp %>%
  ungroup() %>%
  mutate(
    metabolite_identification = as.character(metabolite_identification),
    PrecursorMz = as.character(round(PrecursorMz, digits = 4)),
    ProductMz = as.character(round(ProductMz, digits = 4))
  ) %>%
  left_join(
    precursorChargesAndAdducts %>% ungroup() %>% mutate(
      `Precursor Name` = as.character(`Precursor Name`),
      `Precursor m/z` = as.character(round(`Precursor m/z`, digits = 4)),
      `Product m/z` = as.character(round(`Product m/z`, digits = 4))
    ),
    by = c(
      "metabolite_identification" = "Precursor Name",
      "PrecursorMz" = "Precursor m/z",
      "ProductMz" = "Product m/z"
    )
  ) %>%
  select(
    "metabolite_identification",
    "PrecursorMz",
    "Precursor Charge",
    "PrecursorAdduct",
    fragmentation,
    "ProductMz",
    "RT",
    everything()
  )

#check that the number of metabolite identifications match up with the precursor names and that we did not miss anything
stopifnot(length(
  setdiff(
    j$metabolite_identification,
    precursorChargesAndAdducts$`Precursor Name`
  )
) == 0)

library(httr)
library(jsonlite)
library(xml2)

getSwissLipidsId <- function(originalLipidName) {
  nameTibble <- tibble(metabolite_identification = originalLipidName) %>%
    #mutate(metabolite_identification=gsub("-", "_", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub(" ?\\[NL[0-9:-]+\\]","", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("(.*)O_(.*)p(.*)", "\\1O-\\2\\3", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("(.*)P_(.*)p(.*)", "\\1P-\\2\\3", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("ChE", "CE", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("CE (.*)", "SE(27:1/\\1)", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("SM (.*)", "SM(d\\1)", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("Hex1Cer(.*)", "HexCer\\1", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("([A-Za-z0-9]+) ([0-9]+:[0-9]+);1(.*)", "\\1 m\\2\\3", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("([A-Za-z0-9]+) ([0-9]+:[0-9]+);2(.*)", "\\1 d\\2\\3", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("([A-Za-z0-9]+) ([0-9]+:[0-9]+);3(.*)", "\\1 t\\2\\3", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("TAG (.*)", "TG(\\1)", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("DAG (.*)", "DG(\\1)", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("MAG (.*)", "MG(\\1)", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("([A-Za-z0-9]+) ([OP]?-)(.*)", "\\1(\\2\\3)", metabolite_identification)) %>%
    mutate(metabolite_identification=gsub("([A-Za-z0-9]+) (.*)", "\\1(\\2)", metabolite_identification))
  lipidName <- unique(nameTibble$metabolite_identification)[[1]]
  print(paste0("Resolving ", lipidName, " against SwissLipids"))
  url <-
    URLencode(
      paste(
        #"https://www.swisslipids.org/api/advancedSearch?Name",
        "https://www.swisslipids.org/api/search?term",
        lipidName,
        sep = "="
      )
    )
  response <-
    GET(
      url = url,
      add_headers(Accept = "application/json", "Content-Type" = "application/json")
    )
  if (response$status_code != 200) {
    print(
      paste0(
        "Query for ",
        lipidName,
        " returned with status code: ",
        response$status_code
      )
    )
    tbl <-
      tibble(
        "entity_id" = character(),
        "entity_name" = character(),
        "external_id" = character(),
        "entity_type" = character(),
        "classification_level" = character(),
        "query" = character(),
        "originalLipidName" = character()
      )
    tbl[1, ] <- list("", "", "", "", "", lipidName, originalLipidName)
    return(tbl)
  }
  xmlrl <- content(response, "parsed")
  cnt <- xml_find_all(xmlrl, ".//p/text()")
  itemId <-
    jsonlite::fromJSON(txt = xml_text(cnt[[1]])) %>% as_tibble(.)
  itemId$query <- lipidName
  itemId$originalLipidName <- originalLipidName
  itemIdf <-
    itemId %>% filter(
      classification_level == "Species" |
        classification_level == "Molecular subspecies" |
        classification_level == "Structural subspecies" |
        classification_level == "Isomeric subspecies"
    ) %>%
    select(-stats, -nb_exp_annot) %>% mutate(namelen = nchar(entity_name)) %>% arrange(namelen) %>% select(-namelen)
  print(paste0("Found ", nrow(itemIdf), " hits after filtering"))
  return(itemIdf %>% slice(1))
}

resolveLipidNames <- function(lipidNames) {
  lapply(lipidNames, getSwissLipidsId)
}

normalizedNames <- k %>% ungroup() %>%
  mutate("normalizedName" = metabolite_identification)

mnames <- as.list(normalizedNames$normalizedName)

ln <- resolveLipidNames(mnames)
allNames <- bind_rows(ln)
nrow(allNames)
nentries <- nrow(k)
resolvedLipids <-
  k %>% bind_cols(allNames) %>% ungroup() %>%
  mutate(
    "database_identifier" = "",
    "chemical_formula" = "",
    "smiles" = "",
    "inchi" = "",
    "mass_to_charge" = PrecursorMz,
    "fragmentation" = fragmentation,
    "modifications" = PrecursorAdduct,
    "charge" = !!sym("Precursor Charge"),
    "retention_time" = !!sym("RT"),
    "taxid" = "",
    "species" = "",
    "database" = "SwissLipids",
    "database_version" = "2020/04/02",
    "reliability" = "MSI:2",
    "uri" = "",
    "search_engine" = "",
    "search_engine_score" = "",
    "SwissLipid_identifier" = entity_id
  ) %>%
  select(
    -"entity_id",
    -"entity_name",
    -"external_id",
    -"entity_type",
    -"classification_level",-"PrecursorMz",
    -"ProductMz",
    -"Precursor Charge",
    -"RT",
    -"PrecursorAdduct",
    -"ProductName",
    -"originalLipidName",
    -"query"
  ) %>%
  select(
    "database_identifier",
    "chemical_formula",
    "smiles",
    "inchi",
    "metabolite_identification",
    "mass_to_charge",
    "fragmentation",
    "modifications",
    "charge",
    "retention_time",
    "taxid",
    "species",
    "database",
    "database_version",
    "reliability",
    "uri",
    "search_engine",
    "search_engine_score",
    "SwissLipid_identifier",
    everything()
  )

stopifnot(nrow(resolvedLipids) == nentries)

write_tsv(
  resolvedLipids,
  "m_MTBLS1375_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv"
)
df <- k %>% ungroup()
mim <- as.matrix(df %>% select(-"metabolite_identification",-"PrecursorMz",-"Precursor Charge",-"fragmentation",-"PrecursorAdduct",-"ProductMz",-"RT",-smallmolecule_abundance_sub,-smallmolecule_abundance_stdev_sub,-smallmolecule_abundance_std_error_sub,-"ProductName")) %>% replace_na(0)
rownames(mim) <- df$metabolite_identification
library(pheatmap)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
mim_nm <- t(apply(mim, 1, cal_z_score))
pdf("plasma-heatmap.pdf", width = 14, height = 32)
p <- pheatmap(mim_nm, cutree_cols = 8, cutree_rows = 23, fontsize_row = 8, fontsize_col = 10, width = 8, height = 12)
grid::grid.draw(p)
dev.off()
