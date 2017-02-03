all:
	Rscript -e "devtools::check()"
	Rscript -e "devtools::build()"
	Rscript -e "devtools::install()"
