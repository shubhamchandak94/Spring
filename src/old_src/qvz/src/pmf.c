#include <stdio.h>
#include <string.h>

#include "pmf.h"

//#define log2(a) log(a)/log(2.0)
/**
 * Allocates the memory for an alphabet structure and fills the symbols
 * with a default list of 0 through size-1
 */
struct alphabet_t *alloc_alphabet(uint32_t size) {
	symbol_t i;
	struct alphabet_t *rtn = (struct alphabet_t *) calloc(1, sizeof(struct alphabet_t));
	rtn->size = size;
	rtn->symbols = (symbol_t *) calloc(size, sizeof(symbol_t));

	for (i = 0; i < size; ++i) {
		rtn->symbols[i] = i;
	}
	alphabet_compute_index(rtn);

	return rtn;
}

/**
 * Makes a copy of the given alphabet
 */
struct alphabet_t *duplicate_alphabet(const struct alphabet_t *a) {
	struct alphabet_t *rtn = (struct alphabet_t *) calloc(1, sizeof(struct alphabet_t));
	rtn->size = a->size;
	rtn->symbols = (symbol_t *) calloc(a->size, sizeof(symbol_t));

	memcpy(rtn->symbols, a->symbols, a->size*sizeof(symbol_t));
	alphabet_compute_index(rtn);

	return rtn;
}

/**
 * Allocates a PMF structure for the given alphabet, but it does not copy the alphabet
 */
struct pmf_t *alloc_pmf(const struct alphabet_t *alphabet) {
	struct pmf_t *rtn = (struct pmf_t *) calloc(1, sizeof(struct pmf_t));
	rtn->alphabet = alphabet;
	rtn->pmf = (double *) calloc(alphabet->size, sizeof(double));
	rtn->counts = (uint32_t *) calloc(alphabet->size, sizeof(uint32_t));
	return rtn;
}

/**
 * Allocates an array for tracking a list of PMFs along with the underlying PMFs
 */
struct pmf_list_t *alloc_pmf_list(uint32_t size, const struct alphabet_t *alphabet) {
	uint32_t i;
	struct pmf_list_t *rtn = (struct pmf_list_t *) calloc(1, sizeof(struct pmf_list_t));
	rtn->size = size;
	rtn->pmfs = (struct pmf_t **) calloc(size, sizeof(struct pmf_t *));
	
	for (i = 0; i < size; ++i) {
		rtn->pmfs[i] = alloc_pmf(alphabet);
	}

	return rtn;
}

/**
 * Frees an alphabet
 */
void free_alphabet(struct alphabet_t *alphabet) {
	free(alphabet->symbols);
	free(alphabet->indexes);
	free(alphabet);
}

/**
 * Frees a PMF
 */
void free_pmf(struct pmf_t *pmf) {
	free(pmf->counts);
	free(pmf);
}

/**
 * Frees a list of PMFs
 */
void free_pmf_list(struct pmf_list_t *pmfs) {
	uint32_t i;
	for (i = 0; i < pmfs->size; ++i) {
		free_pmf(pmfs->pmfs[i]);
	}
	free(pmfs->pmfs);
	free(pmfs);
}

/**
 * Determine if a pmf is valid (if it sums to 1, within some tolerance)
 */
uint32_t is_pmf_valid(struct pmf_t *pmf) {
	double sum = 0;
	uint32_t i;

	if (!pmf->pmf_ready)
		recalculate_pmf(pmf);

	for (i = 0; i < pmf->alphabet->size; ++i) {
		sum += pmf->pmf[i];
	}

	if (fabs(sum - 1.0) < 0.0001)
		return 1;
	return 0;
}

/**
 * Gets the probability for a specific location, triggering lazy re-eval if
 * necessary
 */
double get_probability(struct pmf_t *pmf, uint32_t idx) {
	if (!pmf->pmf_ready)
		recalculate_pmf(pmf);
	return pmf->pmf[idx];
}

/**
 * Gets the probability for a specific symbol, triggering lazy re-eval if
 * necessary
 */
double get_symbol_probability(struct pmf_t *pmf, symbol_t symbol) {
	uint32_t idx = get_symbol_index(pmf->alphabet, symbol);

	if (!pmf->pmf_ready)
		recalculate_pmf(pmf);
	if (idx != ALPHABET_SYMBOL_NOT_FOUND)
		return pmf->pmf[idx];
	return 0.0;
}

/**
 * Calculate the entropy of this pmf in bits
 */
double get_entropy(struct pmf_t *pmf) {
	double entropy = 0.0;
	uint32_t i = 0;

	if (!pmf->pmf_ready)
		recalculate_pmf(pmf);

	for (i = 0; i < pmf->alphabet->size; ++i) {
		if (pmf->pmf[i] > 0.0) {
			entropy -= pmf->pmf[i] * log2(pmf->pmf[i]);
		}
	}

	return entropy;
}

/**
 * Calculates the Kullbeck-Leibler Divergence between two PMFs, p and q, as D(p||q)
 */
double get_kl_divergence(struct pmf_t *p, struct pmf_t *q) {
	double d = 0.0;
	uint32_t i;
	
	if (p->alphabet != q->alphabet)
		return NAN;

	if (!p->pmf_ready)
		recalculate_pmf(p);
	if (!q->pmf_ready)
		recalculate_pmf(q);

	for (i = 0; i < p->alphabet->size; ++i) {
		if (q->pmf[i] > 0) {
			if (p->pmf[i] > 0) {
				d += p->pmf[i] * log2(p->pmf[i] / q->pmf[i]);
			}
		}
	}

	return d;
}

/**
 * Combine two PMFs with two weight parameters to scale each before adding. This
 * operates based on the probabilities, not the counts, so it is suitable for use
 * in calculating the law of total probability: p(a)p(X|Y=a) + p(b)p(X|Y=b) when
 * the empirical distributions do not contain the same number of observations
 */
struct pmf_t *combine_pmfs(struct pmf_t *a, struct pmf_t *b, double weight_a, double weight_b, struct pmf_t *result) {
	uint32_t i;

	if (a->alphabet != b->alphabet || a->alphabet != result->alphabet)
		return NULL;
	
	if (!a->pmf_ready)
		recalculate_pmf(a);
	if (!b->pmf_ready)
		recalculate_pmf(b);

	for (i = 0; i < a->alphabet->size; ++i) {
		result->pmf[i] = weight_a * a->pmf[i] + weight_b * b->pmf[i];
	}
	result->pmf_ready = 1;
	return result;
}

/**
 * When counting symbols, this handles incrementing everything for the given
 * index
 */
void pmf_increment(struct pmf_t *pmf, uint32_t index) {
	pmf->counts[index] += 1;
	pmf->total += 1;
}

/**
 * Recalculates the PMF as a series of doubles from the empirical counts and total
 */
void recalculate_pmf(struct pmf_t *pmf) {
	uint32_t i;
	double total = (double) pmf->total;

	pmf->pmf_ready = 1;
	if (pmf->total == 0)
		return;
	
	for (i = 0; i < pmf->alphabet->size; ++i) {
		pmf->pmf[i] = ((double) pmf->counts[i]) / total;
	}
}

/**
 * Renormalizes a PMF if it is nonzero
 */
void renormalize_pmf(struct pmf_t *pmf) {
	double total = 0;
	uint32_t i;

	// PMFs still in counts form never need renormalization
	if (!pmf->pmf_ready)
		return;

	// Find total
	for (i = 0; i < pmf->alphabet->size; ++i) {
		total += pmf->pmf[i];
	}

	// If nonzero, scale every entry to ensure we sum to 1
	if (total > 0) {
		for (i = 0; i < pmf->alphabet->size; ++i) {
			pmf->pmf[i] = pmf->pmf[i] / total;
		}
	}
}

/**
 * Converts a PMF that is stored as a series of doubles back to the counts representation,
 * or alternatively this can be viewed as quantizing it into a fixed point representation in
 * 0.m format
 */
void pmf_to_counts(struct pmf_t *pmf, uint32_t m) {
	uint32_t i;
	double scale = ((1 << m) - 1);

	pmf->total = 0;
	for (i = 0; i < pmf->alphabet->size; ++i) {
		pmf->counts[i] = (uint32_t) (pmf->counts[i] * scale);
		pmf->total += pmf->counts[i];
	}
}

/**
 * Zeros out the counts and probabilities for a PMF to let us reuse the same memory allocation
 */
void clear_pmf(struct pmf_t *pmf) {
	memset(pmf->counts, 0, pmf->alphabet->size * sizeof(uint32_t));
	memset(pmf->pmf, 0, pmf->alphabet->size * sizeof(double));
	pmf->pmf_ready = 0;
	pmf->total = 0;
}

/**
 * Zeros out every pmf in the given list, so we can reuse the entire pmf list without
 * deallocating/reallocating memory
 */
void clear_pmf_list(struct pmf_list_t *list) {
	uint32_t i;
	for (i = 0; i < list->size; ++i) {
		clear_pmf(list->pmfs[i]);
	}
}

/**
 * Determines if the given alphabet contains the given symbol
 */
uint32_t alphabet_contains(const struct alphabet_t *alphabet, symbol_t symbol) {
	return alphabet->indexes[symbol] != ALPHABET_SYMBOL_NOT_FOUND ? 1 : 0;
}

/**
 * Looks up the index of a symbol in the given alphabet, which may be useful
 * if the alphabet doesn't start at zero, has gaps, etc.
 */
uint32_t get_symbol_index(const struct alphabet_t *alphabet, symbol_t symbol) {
	return alphabet->indexes[symbol];
}

/**
 * Finds the unique set of symbols across both input alphabets and creates an
 * output alphabet
 */
void alphabet_union(const struct alphabet_t *restrict a, const struct alphabet_t *restrict b, struct alphabet_t *result) {
	symbol_t *sym = (symbol_t *) _alloca((a->size+b->size)*sizeof(symbol_t));
	uint32_t i = 0;
	uint32_t j = 0;
	uint32_t k = 0;

	// Combine with a merge algorithm since alphabets are required to be sorted
	while (i < a->size && j < b->size) {
		if (a->symbols[i] < b->symbols[j]) {
			sym[k] = a->symbols[i];
			i += 1;
		}
		else if (a->symbols[i] == b->symbols[j]) {
			sym[k] = a->symbols[i];
			i += 1;
			j += 1;
		}
		else {
			sym[k] = b->symbols[j];
			j += 1;
		}
		k += 1;
	}

	// Tail of the merge
	while (i < a->size) {
		sym[k] = a->symbols[i];
		k += 1;
		i += 1;
	}
	while (j < b->size) {
		sym[k] = b->symbols[j];
		k += 1;
		j += 1;
	}

	// If we already have an output array, replace it with a new one
	if (result->symbols)
		free(result->symbols);
	result->symbols = (symbol_t *) calloc(k, sizeof(symbol_t));

	// Copy over temporary data
	memcpy(result->symbols, sym, k*sizeof(symbol_t));
	result->size = k;
	alphabet_compute_index(result);
}

/**
 * Computes the index table (reverse mapping of symbols in the alphabet)
 * that is used to speed up searches for symbols). This isn't a proper
 * hash table and it will consume exponential memory if symbol_t changes
 * size, so be careful
 */
void alphabet_compute_index(struct alphabet_t *A) {
	uint32_t i;

	if (A->indexes)
		free(A->indexes);
	
	// Cheating but whatever
	A->indexes = (uint32_t *) calloc(ALPHABET_INDEX_SIZE_HINT, sizeof(uint32_t));

	// Fill gaps in the table with an appropriate index so we can use this for search too
	for (i = 0; i < ALPHABET_INDEX_SIZE_HINT; ++i) {
		A->indexes[i] = ALPHABET_SYMBOL_NOT_FOUND;
	}

	for (i = 0; i < A->size; ++i) {
		A->indexes[A->symbols[i]] = i;
	}
}

/**
 * Displays an alphabet as "(index): 'character' <number>" one per line
 */
void print_alphabet(const struct alphabet_t *alphabet) {
	uint32_t i;
	for (i = 0; i < alphabet->size; ++i) {
		printf("(%d): '%c' <%d>\n", i, alphabet->symbols[i], alphabet->symbols[i]);
	}
}

/**
 * Displays a PMF
 */
void print_pmf(struct pmf_t *pmf) {
	uint32_t i;

	if (!pmf->pmf_ready)
		recalculate_pmf(pmf);

	for (i = 0; i < pmf->alphabet->size; ++i) {
		printf("<%d>: %.5f (%d/%d)\n", pmf->alphabet->symbols[i], pmf->pmf[i], pmf->counts[i], pmf->total);
	}
}
