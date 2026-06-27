# Plain-make build mirroring CMakeLists.txt (for environments without cmake).
# Shows a [ NN%] progress prefix per compiled target, like the cmake build.
CC      ?= gcc
CFLAGS  ?= -O3 -fomit-frame-pointer -maes -mssse3 -Wall -Wextra -Iinclude -Isrc/keygen
# Fast cut-F keygen check (friend's "truncated tail columns" idea): judge only
# the d x (d+1) leading sub-block (d+1 small d-order minors), ~4-5x faster
# keygen than the default (D-1)-order descent.  Enable with: make DP1=1
ifdef DP1
CFLAGS  += -DUSE_DP1_CHECK
endif
LDLIBS   = -lm

CORE = src/poly.c src/ntt.c src/ntt_tables.c src/sampler.c src/encoding.c \
       src/packing.c src/sign.c src/zpntt.c src/aes256ctr.c src/symmetric-shake.c src/fips202.c \
       src/keygen/kg_main.c src/keygen/kg_solver.c src/keygen/kg_zint.c \
       src/keygen/kg_ntt.c src/keygen/kg_primes.c
RNG  = src/randombytes.c

MODES = 128 256 384 512
BIN   = build/bin

# ---- progress counter -------------------------------------------------------
# All build targets (4 modes x 4 kinds = 16). A shared counter file tracks how
# many have finished so each link line can print [ NN%].
ALLTARGETS = $(foreach m,$(MODES),$(BIN)/test_rhyme$(m) $(BIN)/test_ntt$(m) $(BIN)/kat_rhyme$(m) $(BIN)/speed_rhyme$(m))
NTARGETS  := $(words $(ALLTARGETS))
COUNTER   := build/.progress

# progress: increments counter and prints [ NN%] <msg>
define progress
@mkdir -p build; \
N=$$(cat $(COUNTER) 2>/dev/null || echo 0); N=$$((N+1)); echo $$N > $(COUNTER); \
PCT=$$(( N * 100 / $(NTARGETS) )); \
printf "[%3d%%] %s\n" $$PCT "$(1)"
endef

all: reset tests kats speeds
	@echo "[100%] Built $(NTARGETS) targets."
	@rm -f $(COUNTER)

reset:
	@mkdir -p $(BIN); rm -f $(COUNTER); echo 0 > $(COUNTER)

define MODE_RULES
$(BIN)/test_rhyme$(1): test/main.c $(CORE) $(RNG) | $(BIN)
	$$(call progress,Building test_rhyme$(1))
	@$(CC) $(CFLAGS) -DRHYME_MODE=$(1) -o $$@ $$^ $(LDLIBS)
$(BIN)/test_ntt$(1): test/test_ntt.c $(CORE) $(RNG) | $(BIN)
	$$(call progress,Building test_ntt$(1))
	@$(CC) $(CFLAGS) -DRHYME_MODE=$(1) -Itest -o $$@ $$^ $(LDLIBS)
$(BIN)/kat_rhyme$(1): kat/PQCgenKAT_sign.c kat/nist_rng.c $(CORE) | $(BIN)
	$$(call progress,Building kat_rhyme$(1))
	@$(CC) $(CFLAGS) -DRHYME_MODE=$(1) -Ikat -o $$@ $$^ $(LDLIBS)
$(BIN)/speed_rhyme$(1): benchmark/speed.c benchmark/cpucycles.c $(CORE) $(RNG) | $(BIN)
	$$(call progress,Building speed_rhyme$(1))
	@$(CC) $(CFLAGS) -DRHYME_MODE=$(1) -Ibenchmark -o $$@ $$^ $(LDLIBS)
endef
$(foreach m,$(MODES),$(eval $(call MODE_RULES,$(m))))

tests:  $(foreach m,$(MODES),$(BIN)/test_rhyme$(m) $(BIN)/test_ntt$(m))
kats:   $(foreach m,$(MODES),$(BIN)/kat_rhyme$(m))
speeds: $(foreach m,$(MODES),$(BIN)/speed_rhyme$(m))

$(BIN):
	@mkdir -p $(BIN)

check: tests
	@for m in $(MODES); do $(BIN)/test_ntt$$m || exit 1; done
	@$(BIN)/test_rhyme128 && $(BIN)/test_rhyme256

clean:
	rm -rf build
.PHONY: all reset tests kats speeds check clean
