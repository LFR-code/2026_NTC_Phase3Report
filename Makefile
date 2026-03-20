.PHONY: pdf
pdf:
	Rscript -e "bookdown::render_book('index.Rmd')"

.PHONY: check
check:
	@echo "Checking for common issues..."
	@echo "1. Checking for TODO markers..."
	@grep -n "TODO" 0*.Rmd || echo "  ✓ No TODO markers found"
	@echo "2. Checking line lengths (should be ≤80 chars)..."
	@awk 'length>80 {print FILENAME":"NR": Line too long ("length" chars)"; count++} END {if (count) print "  ⚠ "count" lines exceed 80 characters"; else print "  ✓ All lines ≤80 characters"}' 0*.Rmd
	@echo "3. Checking for placeholder text..."
	@grep -n "placeholder\|lorem ipsum\|TODO\|FIXME\|XXX" index.Rmd 0*.Rmd || echo "  ✓ No placeholder text found"

.PHONY: clean
clean:
	rm -rf _book/ _bookdown_files/ *.log *.aux *.out *.tex
