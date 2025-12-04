# Makefile for LaTeX Template
# Combines simplicity with essential features

# Configuration
MAIN = main
BUILD_DIR = build

# Detect OS for PDF viewer
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    OPEN = xdg-open
endif
ifeq ($(UNAME_S),Darwin)
    OPEN = open
endif
ifeq ($(OS),Windows_NT)
    OPEN = start
endif

.PHONY: all clean rebuild view help

# Default: full build with glossaries
all: $(BUILD_DIR)/$(MAIN).pdf $(MAIN).pdf
	@echo "✓ Build complete: $(MAIN).pdf"

# Main build rule with proper dependencies
$(BUILD_DIR)/$(MAIN).pdf: $(MAIN).tex
	@echo "Creating build directory..."
	@mkdir -p $(BUILD_DIR)
	@echo "Compiling LaTeX document (pass 1/3)..."
	@pdflatex -interaction=nonstopmode -output-directory=$(BUILD_DIR) $(MAIN).tex > $(BUILD_DIR)/build.log 2>&1 || (echo "❌ Error in pass 1! Check $(BUILD_DIR)/build.log for details"; exit 1)
	@echo "Processing bibliography..."
	@bibtex $(BUILD_DIR)/$(MAIN) >> $(BUILD_DIR)/build.log 2>&1 || true
	@echo "Processing glossaries and index..."
	@makeglossaries -d $(BUILD_DIR) $(MAIN) >> $(BUILD_DIR)/build.log 2>&1 || true
	@echo "Compiling LaTeX document (pass 2/3)..."
	@pdflatex -interaction=nonstopmode -output-directory=$(BUILD_DIR) $(MAIN).tex >> $(BUILD_DIR)/build.log 2>&1 || (echo "❌ Error in pass 2! Check $(BUILD_DIR)/build.log for details"; exit 1)
	@echo "Compiling LaTeX document (pass 3/3)..."
	@pdflatex -interaction=nonstopmode -output-directory=$(BUILD_DIR) $(MAIN).tex >> $(BUILD_DIR)/build.log 2>&1 || (echo "❌ Error in pass 3! Check $(BUILD_DIR)/build.log for details"; exit 1)
	@echo "Checking for LaTeX warnings..."
	@grep -i "warning" $(BUILD_DIR)/build.log | head -5 || true

# Copy PDF to main directory
$(MAIN).pdf: $(BUILD_DIR)/$(MAIN).pdf
	@cp $(BUILD_DIR)/$(MAIN).pdf $(MAIN).pdf

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(BUILD_DIR)
	@rm -f $(MAIN).pdf *.aux *.log *.out *.toc *.bbl *.blg *.synctex.gz *.fdb_latexmk *.fls

# Rebuild from scratch
rebuild: clean all

# Compile and open PDF
view: all
	@echo "Opening $(MAIN).pdf..."
	@$(OPEN) $(MAIN).pdf 2>/dev/null || echo "Please open $(MAIN).pdf manually"

# Show help
help:
	@echo "LaTeX Template Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  make          - Build the document (default)"
	@echo "  make clean    - Remove build artifacts"
	@echo "  make rebuild  - Clean and rebuild"
	@echo "  make view     - Build and open PDF"
	@echo "  make help     - Show this help"
	@echo ""
	@echo "Output: $(MAIN).pdf"
