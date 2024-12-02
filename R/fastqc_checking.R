# Alvaro Benitez, 2021

qc.dir <- system.file("test_raw", package = "fastqcr")
qc.dir
list.files(qc.dir)
qc <- qc_aggregate(qc.dir)
qc %>%
  select(sample, module, status) %>%
  filter(status %in% c("WARN", "FAIL")) %>%
  arrange(sample)
summary(qc)
qc_stats(qc)
df <- qc_stats(qc)
setwd("~/fastqc_output")
write.table(df, file = "resumen_raw.csv", quote = FALSE, sep = "\t", row.names = F, col.names = T)
