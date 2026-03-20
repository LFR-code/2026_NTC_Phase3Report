bookdown::render_book(input = "index.Rmd",
                      clean = TRUE,
                      output_format = "bookdown::pdf_document2",
                      config_file = "_bookdown.yml")

bookdown::render_book(input = "index.Rmd",
                      clean = TRUE,
                      output_format = "bookdown::word_document2",
                      config_file = "_bookdown.yml")
