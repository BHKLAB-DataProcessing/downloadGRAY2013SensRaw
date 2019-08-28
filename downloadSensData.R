library(downloader)
library(PharmacoGx)

options(stringsAsFactors=FALSE)

options(stringsAsFactors=FALSE)
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[\\]|[.]|[_]|[ ]"

#match to cell curations

matchToIDTableCELL <- function(ids,tbl, column) {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating cell ids")
    }
    return(tbl[myx, "unique.cellid"])
  })
}


#match to drug curation 

matchToIDTableDRUG <- function(ids,tbl, column) {
  sapply(ids,function(x) {
    myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating drug ids")
    }
    return(tbl[myx, "unique.drugid"])
  })
}


cell_all <- read.csv(file = "/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
curationCell <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
curationCell <- curationCell[ , c("unique.cellid", "GRAY.cellid")]
rownames(curationCell) <- curationCell[ , "unique.cellid"]
removed_cells <- c("DU-4475","T47D_KBluc","HCC2157","MB157","HCC1500","SUM 190","HCC1007","184A1N4", "MDA-MB-435")
curationCell <- curationCell[which(!curationCell$unique.cellid %in% removed_cells),]
rownames(curationCell) <- curationCell$unique.cellid

curationTissue <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
curationTissue <- curationTissue[which(curationTissue$unique.cellid %in% curationCell$unique.cellid),]
curationTissue <- curationTissue[ , c("unique.tissueid", "GRAY.tissueid")]
rownames(curationTissue) <- curationCell[ , "unique.cellid"]

drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
curationDrug <- drug_all[which(!is.na(drug_all[ , "GRAY.drugid"])),]
curationDrug <- curationDrug[,c("unique.drugid","GRAY.drugid")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]

removed_drugs <- c("Nilotinib","CI-1040","PD 173074","GW843682X","Dichloroacetic acid","2-deoxyglucose","Celecoxib","LY-294002","Sulindac sulfide","Chloroquine","Ribavirin","TPCA-1","Cetuximab","Bromopyruvic acid","Methyl Glyoxol","Z-LL-Nva-CHO","TPT","Trastuzumab",        
                   "UO126","PS-1145","Tyrphostin AG 1024","GW5074","Ilomastat","QNZ",             
                   "TAPI-0","2C4","L779450","Chk2 Inhibitor II","Ro 32-0432","TCS JNK 5a",         
                   "PP242 hydrate")

curationDrug <- curationDrug[which(!curationDrug$unique.drugid %in% removed_drugs),]


getGRAYrawData <-
  function(result.type=c("array", "list")){
    
    gray.raw.drug.sensitivity <- read.csv(file = "/pfs/downloadGRAY2013SensRaw/gb-2013-14-10-r110-s9.txt" ,header = TRUE, sep = "\t",stringsAsFactors = FALSE)
    gray.raw.drug.sensitivity.list <- do.call(c, apply(gray.raw.drug.sensitivity, 1, list))
    
    concentrations.no <- length(grep("^c[1-9]", colnames(gray.raw.drug.sensitivity)))
    
    if(result.type == "array"){
      ## create the gray.drug.response object including information viablilities and concentrations for each cell/drug pair
      obj <- array(NA, dim=c(length(unique(gray.raw.drug.sensitivity[ , "cellline"])), length(unique(gray.raw.drug.sensitivity[ , "compound"])), 2, concentrations.no), dimnames=list(unique(gray.raw.drug.sensitivity[ , "cellline"]), unique(gray.raw.drug.sensitivity[ , "compound"]), c("concentration", "viability"), 1:concentrations.no))
    }
    fnexperiment <- 
      function(values)  {
        cellline <- values["cellline"]
        drug <- values["compound"]
        #doses <- as.numeric(unlist(strsplit(input.matrix["Doses (uM)"], split=", "))) #nature paper raw data
        doses <- as.numeric(values[grep("^c[1-9]", colnames(gray.raw.drug.sensitivity))]) * 10 ^ 6 # micro molar
        if(concentrations.no > length(doses)) {doses <- c(doses, rep(NA, concentrations.no - length(doses)))}
        
        #responses <- as.numeric(unlist(strsplit(input.matrix["Activity Data\n(raw median data)"], split=",")))  #nature paper raw data
        
        responses <- matrix(NA, ncol=concentrations.no, nrow=3)
        for( i in 1:concentrations.no)
        {
          responses[1, i] <- as.numeric(values[sprintf("od%s.%s",i,1)])/as.numeric(values[sprintf("od%s.%s",0,1)])
          responses[2, i] <- as.numeric(values[sprintf("od%s.%s",i,2)])/as.numeric(values[sprintf("od%s.%s",0,2)])
          responses[3, i] <- as.numeric(values[sprintf("od%s.%s",i,3)])/as.numeric(values[sprintf("od%s.%s",0,3)])
        }
        responses <- apply(responses, MARGIN=2 , median, na.rm=TRUE) * 100
        
        if(result.type == "array"){
          obj[cellline,drug, "concentration", 1:length(doses)] <<- doses
          obj[cellline,drug, "viability", 1:length(responses)] <<- responses
        }else{
          return(list(cell=cellline, drug=drug, doses=doses, responses=responses))#paste(doses, collapse = ","), responses=paste(responses, collapse = ",")))
        }
      }
    
    gray.raw.drug.sensitivity.res <- mapply(fnexperiment, values=gray.raw.drug.sensitivity.list)
    if(result.type == "array"){
      return(list("data"=obj, "concentrations.no"=concentrations.no))
    }else{
      return(list("data"=gray.raw.drug.sensitivity.res, "concentrations.no"=concentrations.no))
    }
  }


raw.sensitivity <- getGRAYrawData(result.type="list")

con_tested <- raw.sensitivity$concentrations.no
raw.sensitivity <- t(raw.sensitivity$data)
raw.sensitivity <- t(apply(raw.sensitivity,1, function(x){unlist(x)}))

raw.sensitivity[ ,2] <- gsub("\\s*\\([^\\)]+\\)","",raw.sensitivity[ ,2]) #REMOVES PARATHESES JUST FOR MATCHING
raw_drug_match <- matchToIDTableDRUG(raw.sensitivity[ ,2], curationDrug, "GRAY.drugid")
raw.sensitivity[ ,2] <- as.character(raw_drug_match)

curationDrug <- curationDrug[which(curationDrug$unique.drugid %in% unique(raw_drug_match)),]

raw.sensitivity[ ,1][grep("KBluc", raw.sensitivity[ ,1])] <- "T47D"
raw.sensitivity[ ,1][grep("157", raw.sensitivity[ ,1])] <- "MDAMB157"

raw_cell_match <- matchToIDTableCELL(raw.sensitivity[ ,1], curationCell, "GRAY.cellid")
raw.sensitivity[ ,1] <- as.character(raw_cell_match)

rownames(raw.sensitivity)  <- sprintf("%s_%s",as.character(raw.sensitivity[ ,2]),as.character(raw.sensitivity[ ,1]))
## handle replicates
tt <- rownames(raw.sensitivity)
for(i in 1:length(tt)) {
  xx <- which(tt == tt[i])
  if(length(xx) > 1) {
    for(j in 1:length(xx)) {
      tt[xx[j]] <- paste(tt[xx[j]], j, sep="_")
    }
  }
}
rownames(raw.sensitivity) <- tt

tt <- paste0("doses", con_tested)



sensitivity.info <- raw.sensitivity[ , c(1,2, 3, grep(tt, colnames(raw.sensitivity)))]
colnames(sensitivity.info) <- c("cellid", "drugid", "min.Dose.uM", "max.Dose.uM")
sensitivity.info <- cbind(sensitivity.info, "nbr.conc.tested"=con_tested)


raw.sensitivity <- raw.sensitivity[ ,-c(1,2)]
raw.sensitivity <- array(c(as.matrix(raw.sensitivity[ ,1:con_tested]), as.matrix(raw.sensitivity[ ,(con_tested+1):(2*con_tested)])), c(nrow(raw.sensitivity), con_tested, 2),
                         dimnames=list(rownames(raw.sensitivity), colnames(raw.sensitivity[ ,1:con_tested]), c("Dose", "Viability")))

save(raw.sensitivity, sensitivity.info, tt, con_tested, file="/pfs/out/drug_norm_post.RData")


raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/1000))

dir.create("/pfs/out/slices/")

for(i in seq_along(raw.sensitivity.x)){

  slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
  saveRDS(slce, file=paste0("/pfs/out/slices/gray2013_raw_sens_", i, ".rds"))

}
