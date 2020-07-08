all: slides purl copy

slides: advanced-vegan.Rmd slides.css
	Rscript -e "rmarkdown::render(\"advanced-vegan.Rmd\")"

purl: advanced-vegan.Rmd
	Rscript -e "knitr::purl(\"advanced-vegan.Rmd\")"

copy: advanced-vegan.html slides.css macros.js
	cp -R -u advanced-vegan_files advanced-vegan.html macros.js slides.css libs resources ~/work/web/jekyll/blog/slides/advanced-vegan-webinar-2020/
