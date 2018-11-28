SHELL := /bin/bash
COOL_URL=ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122958/suppl/GSE122958%5Fjurkat%5Fwt%5Fhg19%5F5k%5Fq10%2Ecool%2Egz
COOL_FILE=cool/jurkat_wt_hg19_5k_q10.cool

all: run-docker-analysis

run-docker-analysis:
	if [[ "$(shell which docker)" == "" ]]; then \
		echo "error: Docker is required to run this pipeline."; \
		exit 1; \
	fi
	docker build -t hic-analysis .
	docker run -v $(CURDIR):/analysis hic-analysis make -C /analysis analysis

clusters: clusters/global_clusters.txt

$(COOL_FILE):
	wget $(COOL_URL) -O - | gzip -d > $(COOL_FILE)

clusters/global_clusters.txt: $(COOL_FILE)
	cd clusters && python src/global_clusters_oe_ice.py ../$(COOL_FILE)

analysis: clusters/global_clusters.txt
	cd analysis && R -e 'library(rmarkdown); rmarkdown::render("maja/Master.Rmd","html_document");'
	cd analysis && R -f src/clust_analysis_multiple.R
