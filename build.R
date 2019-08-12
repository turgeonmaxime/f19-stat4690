# Build slides----

# The functions below are taken from blogdown
require_rebuild = function(html, rmd) {
  older_than(html, rmd) || length(readLines(html, n = 1)) == 0
}

older_than = function(file1, file2) {
  !file_test('-f', file1) | file_test('-ot', file1, file2)
}

is_draft = function(rmd) {
  x <- rmarkdown::yaml_front_matter(rmd)$draft
  ifelse(is.null(x), FALSE, x)
}

##############################################

