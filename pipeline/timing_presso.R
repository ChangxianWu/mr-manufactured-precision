.libPaths("D:/Rlibs/4.6")
suppressPackageStartupMessages({library(data.table);library(MRPRESSO)})
source("D:/Horzontal_project/20B3/work/course_config.R")
O <- data.table::fread(file.path(DATA,"outcome_extract/outcome_extract_SEPSIS.tsv"), data.table=FALSE)
O$beta_out<-as.numeric(O$beta_out); O$se_out<-as.numeric(O$se_out); O$eaf_out<-as.numeric(O$eaf_out)
g <- O[O$taxon=="class.Actinobacteria.id.419",]
cat("k =", nrow(g), "\n")
d <- data.frame(by=as.numeric(g$beta_out), bx=as.numeric(g$beta_instr),
                sy=as.numeric(g$se_out), sx=as.numeric(g$se_instr))
t0<-Sys.time()
r <- MRPRESSO::mr_presso(BetaOutcome="by",BetaExposure="bx",SdOutcome="sy",SdExposure="sx",
      OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=d,NbDistribution=500,SignifThreshold=0.05)
cat("elapsed sec:", as.numeric(difftime(Sys.time(),t0,units="secs")), "\n")
print(r$`Main MR results`)
