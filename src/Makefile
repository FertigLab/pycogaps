SOURCES = $(filter-out Rpackage/src/RcppExports.cpp, $(wildcard Rpackage/src/*.cpp)) \
	$(wildcard Rpackage/src/math/*.cpp) \
	$(wildcard Rpackage/src/data_structures/*.cpp) \
	$(wildcard Rpackage/src/file_parser/*.cpp)

HEADERS = $(wildcard Rpackage/src/*.h) $(wildcard Rpackage/src/math/*.h) \
	$(wildcard Rpackage/src/data_structures/*.h) \
	$(wildcard Rpackage/src/file_parser/*.h)

# TODO get these paths from R
INCLUDES = -I/usr/share/R/include -I/usr/local/lib/R/site-library/Rcpp/include \
	-I/usr/local/lib/R/site-library/BH/include

DEFINES = -DNDEBUG -DSIMD -DNGAPS_DEBUG -DBOOST_MATH_PROMOTE_DOUBLE_POLICY=0
CFLAGS = -fpic -g -O2 -fstack-protector-strong -Wformat -Wdate-time -fopenmp=libomp \
	-msse4.2 -Werror=format-security -D_FORTIFY_SOURCE=2 -Werror -x c++

LINT_TARGETS = $(addsuffix .lint, $(SOURCES)) $(addsuffix .lint, $(HEADERS))

CLANG_TIDY_CHECKS := *,-llvm-header-guard,-android*,-cppcoreguidelines-pro-bounds-pointer-arithmetic,-google-readability-todo,-modernize-use-auto,-cppcoreguidelines-pro-type-vararg,-google-runtime-references,-modernize-loop-convert

LINT_COLOR := \033[0;34m
NO_COLOR := \033[m

## Targets for R package

install_R : build_R
	cd R_build && \
	R CMD INSTALL --configure-args="$(CONFIG_ARGS)" CoGAPS_*.tar.gz && \
	cd ..

R_configure_script :
	autoconf Rpackage/configure.ac > Rpackage/configure
	rm -rf autom4te.cache

build_R :
	mkdir -p R_build
	cd R_build && \
	rm -f CoGAPS_*.tar.gz && \
	R CMD build --no-build-vignettes ../Rpackage && \
	cd ..

build_R_with_vignettes :
	mkdir -p R_build
	cd R_build && \
	rm -f CoGAPS_*.tar.gz && \
	R CMD build ../Rpackage && \
	cd ..

vignette : build_R_with_vignettes
	cd R_build && \
	rm -rf vignette_temp && \
	mkdir vignette_temp && \
	tar -xf CoGAPS_*.tar.gz -C vignette_temp && \
	cp vignette_temp/CoGAPS/inst/doc/CoGAPS_Workflow.html . && \
	rm -rf vignette_temp && \
	cd ..

docs_R : RcppInterface
	R -e 'roxygen2::roxygenise("./Rpackage")'

RcppInterface :
	R -e 'Rcpp::compileAttributes("./Rpackage")'

test_R_package :
	R -e 'devtools::load_all("./Rpackage"); devtools::test("./Rpackage")

check_R_package : build_R_with_vignettes
	cd R_build && \
	R CMD check CoGAPS_*.tar.gz && \
	R CMD BiocCheck CoGAPS_*.tar.gz	&& \
	cd ..

example_R :
	R -e 'library(CoGAPS); data(GIST); eat <- CoGAPS(as.matrix(GIST.D), as.matrix(GIST.S), nIter=1000, nOutput=250)'

example_R_from_file :
	R -e 'library(CoGAPS); data(GIST); eat <- CoGAPSFromFile("SampleData/GIST.mtx", NULL, nIter=1000, nOutput=250)'

set_R_compiler_gcc :
	@printf "CC=ccache gcc\nCXX=ccache g++\n" > ~/.R/Makevars

set_R_compiler_clang :
	@printf "CC=ccache clang\nCXX=ccache clang++\n" > ~/.R/Makevars

install_R_debug : build_R
	@make install_R CONFIG_ARGS=--enable-debug

install_R_with_warnings : build_R
	@make install_R CONFIG_ARGS=--enable-warnings

test_R_install : build_R
	@make install_R
	@make install_R CONFIG_ARGS=--enable-warnings
	@make install_R CONFIG_ARGS=--enable-warnings,--disable-simd
	@make install_R CONFIG_ARGS=--enable-warnings,--enable-simd
	@make install_R CONFIG_ARGS=--enable-warnings,--enable-simd=sse
	@make install_R CONFIG_ARGS=--enable-warnings,--enable-simd=avx

test_R_gcc_install : set_R_compiler_gcc
	@make test_R_install --no-print-directory

test_R_clang_install : set_R_compiler_clang
	@make test_R_install --no-print-directory

lint_R :
	R -e 'eat <- sapply(list.files("R_Package/R"), function(f) lintr::lint(paste("R_Package/R/", f, sep=""))); if (is.null(warnings())) {message("Lint Free!")} else {stop("Error: Not Lint Free")}'

lint_R_Package_cpp : $(R_CPP_LINT_TARGETS)

lint_R_Package_cpp_with_errors :
	@make LINT_FLAGS=-warnings-as-errors='*' lint_R_Package_cpp

lint_R_Package : lint_R_Package_cpp_with_errors lint_R

full_R_test_suite : test check lintR lintCPP_with_errors
	@make test_R_gcc_install
	@make install_R_debug
	@make example_R
	@make test_R_clang_install
	@make install_R_debug
	@make example_R

%.cpp.lint : %.cpp
	@echo "$(LINT_COLOR)linting" $(basename $@) "$(NO_COLOR)"
	@clang-tidy $(basename $@) $(LINT_FLAGS) -checks=$(CLANG_TIDY_CHECKS) -- $(INCLUDES) \
	$(DEFINES) $(CFLAGS)

%.h.lint : %.h
	@echo "$(LINT_COLOR)linting" $(basename $@) "$(NO_COLOR)"
	@clang-tidy $(basename $@) $(LINT_FLAGS) -checks=$(CLANG_TIDY_CHECKS) -- $(INCLUDES) \
	$(DEFINES) $(CFLAGS)

help :
	@echo "make install................install CoGAPS"
	@echo "make install_debug..........install the debug version of CoGAPS "
	@echo "make install_warnings.......install CoGAPS with compiler warnings enabled"
	@echo "make test_install...........test installation success for gcc/clang"
	@echo "make test_gcc_install.......test installation success for gcc"
	@echo "make test_clang_install.....test installation success for clang"
	@echo "make build..................build the package tarball"
	@echo "make build_w_vignettes......build the package tarball with vignettes"
	@echo "make vignette...............render vignette pdf/html"
	@echo "make configure_script.......run configure.ac through autoconf"
	@echo "make set_compiler_gcc.......set the R compiler to gcc"
	@echo "make set_compiler_clang.....set the R compiler to clang"
	@echo "make full_test_suite........run all tests on CoGAPS"
	@echo "make check..................run Bioconductor package checks"
	@echo "make test...................run unit tests"
	@echo "make example................run short CoGAPS example"
	@echo "make docs...................generate documentation files"
	@echo "make RcppInterface..........generate RcppExports R and C++ file"
	@echo "make lintR..................lint R code"
	@echo "make lintCPP................lint C++ code"
	@echo "make lintCPP_with_errors....lint C++ code with failure on warnings"
