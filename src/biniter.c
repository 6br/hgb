
/*
 * invariant conditions for bin composition:
 * 1. the first bin is always placed at offset zero of the `bins` array.
 * 2. the length (range) of the first bin is always 2^63 bp, so any genomic
 * position (any value of `genomic_pos`) is contained in the range of
 * the first bin.
 * 3. inferred from the second condition, the depth 0 (the depth with
 * the longest span) always consists of a single bin.
 */

#include <stdlib.h>
#include <stdint.h>
#include <x86intrin.h>

typedef struct {
	size_t start, end;			/* [start, end) */
} genomic_range_t;

typedef struct {
	size_t ofs;
	size_t len;
} bin_header_t;

typedef struct {
	bin_header_t const *ptr;
	size_t cnt;
} bin_slice_t;

typedef struct {
	genomic_range_t range;		/* (start pos of the first bin, and what?) */
	size_t bin_size;
} bin_slice_info_t;

typedef struct {
	/*
	 * bin_count_mask describes depth and bin count mappings.
	 * 
	 * description by example, showing least significant 16bits:
	 *
	 *   three bits are populated in the `bin_count_mask`, at 0, 5, and 7.
	 *     bin_count_mask: 0b...0000 1   0   1   0   0   0   0   1
	 *
	 *       the number of bins in each depth is inferred from the position of the bits.
	 *       the first depth, which corresponds to the least significant set bit (the rightmost bit),
	 *       has only one (2^0) bin, because it is located at bit 0.
	 *       the second depth, which corresponds to the second set bit at 5th,
	 *       indicates the second depth has 2^5 bins, and so on (2^7) for the third (last) depth.
	 *
	 *
	 *   these three bits corresponds to the first three elements of the `bin_pitch_indices` array.
	 *     bin_pitch_indices in reverse order:
	 *                  [0, ..., 0, 15,     18,                 62]
	 *                    (note: unset bits in the `bin_count_mask` do not correspond to any element
	 *                     in the `bin_pitch_indices` array. the zero elements are just unused and left zero.)
	 *
	 *       the first element 62 tells bin(s) in the first depth has a span of 2^(62+1) bp,
	 *       because bin size is twice as large as bin pitch for every layer.
	 *       the second element 18 tells all the bins in the second depth have a span of 2^(18+1),
	 *       and the third element 15 tells a span of 2^(15+1).
	 */
	uint64_t bin_count_mask;
	uint8_t bin_pitch_indices[64];
	bin_header_t bins[];
} bins_t;

typedef struct {
	genomic_range_t range;	/* immutable (just storing input), packed coordinate */
	uint64_t finished;		/* bitmask, updated from `bin_count_mask` and previous `finished` */
} bin_iterator_t;

static inline
bin_iterator_t bin_iterator_init(bins_t const *bins, genomic_range_t range)
{
	return((bin_iterator_t){
		.range    = range,
		.finished = 0	/* zero because nothing is done */
	});
}

typedef struct {
	bin_slice_t slice;
	bin_slice_info_t info;
} bin_iterator_next_t;

static inline
bin_iterator_next_t bin_iterator_next(bins_t const *bins, bin_iterator_t *it)
{
	/* load constant and previous state */
	uint64_t const finished   = it->finished;
	uint64_t const everything = bins->bin_count_mask;

	/* mask finished bits out to compute remaining bits */
	uint64_t const remaining  = everything & ~finished;

	/* if there is no remaining bits, iteration is done */
	if(remaining == 0) {
		return((bin_iterator_next_t){
			.slice = { .ptr   = NULL,     .cnt      = 0 },
			.info  = { .range = { 0, 0 }, .bin_size = 0 }
		});
	}

	/* compute bin offset; offset base is the sum of the number of finished bins */
	size_t const bin_ofs_base = (size_t)(everything & finished);

	/* determine bin pitch and size; note that bin span (size) is twice as large as bin pitch, because bins are half-overlapping */
	size_t const iterations_done = _mm_popcnt_u64(bin_ofs_base);
	size_t const bin_pitch_index = (size_t)bins->bin_pitch_indices[iterations_done];
	size_t const bin_pitch = (size_t)(0x01ULL<<bin_pitch_index);
	size_t const bin_size  = 2 * bin_pitch;

	/* compute where to slice */
	size_t const bin_ofs_disp_start = ({
		size_t t = it->range.start>>bin_pitch_index;
		if(t > 0) { t--; }						/* make the range margined so that any read overlaps with the input genomic range is collected */
		t;
	});
	size_t const bin_ofs_disp_end = ({
		size_t t = it->range.end>>bin_pitch_index;
		if(t > finished) { t = finished; }		/* previous finished (not new one) is <bin count for this depth> - 1 */
		t + 1;									/* make the range margined (tail margin) */
	});

	/* compute bin range */
	size_t const bin_range_start =  bin_ofs_disp_start   <<bin_pitch_index;
	size_t const bin_range_end   = (bin_ofs_disp_end + 1)<<bin_pitch_index;

	/* update finished mask; find least significant set bit to locate the next unfinished bit, then xor to squish further remaining bits out */
	it->finished = remaining ^ (remaining - 1);

	/* done; compute bin pointer */
	return((bin_iterator_next_t){
		.slice = {
			.ptr = &bins->bins[bin_ofs_base + bin_ofs_disp_start],
			.cnt = bin_ofs_disp_end - bin_ofs_disp_start
		},
		.info = {
			.range = {
				.start = bin_range_start,
				.end   = bin_range_end
			},
			.bin_size = bin_size
		}
	});
}


