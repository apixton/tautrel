# Sage code by Aaron Pixton (apixton@math.harvard.edu)
#
# This is a work in progress consisting of tools to do calculations in
# tautological rings of moduli spaces of curves.

# initialize a few global variables
attach "initialize.sage"

# non-mathematical functions
attach "meta.sage"

# miscellaneous mathematical functions
attach "util.sage"

# some linear algebra functions (e.g. computing rank of a matrix)
attach "linear_algebra.sage"

# basic strata functionality - representing strata classes, listing strata
# classes, etc
attach "strata.sage"

# code to compute the 3-spin relations
attach "3spin.sage"

# code to multiply strata classes and compute the socle pairing
attach "gorenstein.sage"

# more strata functionality (e.g. pullback map)
attach "strata_operations.sage"

# code to find new relations inside the 3-spin relations
attach "primitive.sage"

# specialized code for powers of the universal curve
attach "curvepower.sage"

# some tests
attach "tests.sage"