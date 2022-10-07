require(rmarkdown)
Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")
rmarkdown::render('gene_DE.Rmd')