#include <assert.h>
#include <stdio.h>
#include "qv_compressor.h"

Arithmetic_code initialize_arithmetic_encoder(uint32_t m) {
    Arithmetic_code a_code;
    
    a_code = (Arithmetic_code) calloc(1, sizeof(struct Arithmetic_code_t));
    
    a_code->m = m;
	a_code->r = 1 << (m - 3);
    a_code->l = 0;
	a_code->u = (1 << m) - 1;
    
    return a_code;
}

/**
 * E1/E2 check for the MSB of the lower and upper regions being the same, indicating that a bit has
 * been determined and must be sent to the output stream
 * E3 checks for upper being 10xxxx... and lower being 01xxxx... indicating that after rescaling the
 * range we are still in the indetermined central region
*/
void arithmetic_encoder_step(Arithmetic_code a, stream_stats_ptr_t stats, int32_t x, osStream os) {
    uint64_t range = 0;
    uint8_t msbU = 0, msbL = 0, E1_E2 = 0, E3 = 0, smsbL = 0, smsbU = 0;
    uint32_t cumCountX, cumCountX_1;
    int32_t i;

	// These are actually constants, need to lift a->m out of the struct because it is compile-time constant
	uint32_t msb_shift = a->m - 1;
	uint32_t smsb_shift = a->m - 2;
	uint32_t msb_clear_mask = (1 << msb_shift) - 1;
    
    range = a->u - a->l + 1;

	assert(x < stats->alphabetCard);
    
	cumCountX_1 = 0;
	for (i = 0; i < x; ++i) {
		cumCountX_1 += stats->counts[i];
	}
	cumCountX = cumCountX_1 + stats->counts[x];

	assert(cumCountX_1 < cumCountX);
    
    a->u = a->l + (uint32_t)((range * cumCountX) / stats->n) - 1;
    a->l = a->l + (uint32_t)((range * cumCountX_1) / stats->n);
    
	assert(a->l <= a->u);
    
    // Check the rescaling conditions
    msbL = a->l >> msb_shift;
    msbU = a->u >> msb_shift;
    E1_E2 = (msbL == msbU);
	E3 = 0;
    
    if (!E1_E2) {
		smsbL = a->l >> smsb_shift;
		smsbU = a->u >> smsb_shift;
		E3 = (smsbL == 0x01 && smsbU == 0x02);
    }
    
	// While the bounds need rescaling
    while (E1_E2 || E3) {
        if (E1_E2) {
			// We are in one half of the integer range so the next bit is fixed as the MSB
			stream_write_bit(os, msbL);

			// Clear the msb from both bounds and rescale them
			a->l = (a->l & msb_clear_mask) << 1;
			a->u = ((a->u & msb_clear_mask) << 1) + 1;
            
			// Write any extra bits based on the number of rescalings without an output before now
            while (a->scale3 > 0) {
				stream_write_bit(os, !msbL);
                a->scale3 -= 1;
            }
        }
		else { // E3 is true
            a->scale3 += 1;
			a->u = (((a->u << 1) & msb_clear_mask) | (1 << msb_shift)) + 1;
			a->l = (a->l << 1) & msb_clear_mask;
        }

        msbL = a->l >> msb_shift;
        msbU = a->u >> msb_shift;
        E1_E2 = (msbL == msbU);
		E3 = 0;

        if (!E1_E2) {
			smsbL = a->l >> smsb_shift;
			smsbU = a->u >> smsb_shift;
			E3 = (smsbL == 0x01 && smsbU == 0x02);
        }
    }
}

int encoder_last_step(Arithmetic_code a, osStream os) {
    uint8_t msbL = a->l >> (a->m - 1);

    // Write the msb of the tag (l)
	stream_write_bit(os, msbL);
    
    // write as many !msbL as scale3 left
    while (a->scale3 > 0) {
		stream_write_bit(os, !msbL);
        a->scale3 -= 1;
    }
    
    // write the rest of the tag (l)
	stream_write_bits(os, a->l, a->m - 1);
	stream_finish_byte(os);
    
    return os->written;
}

uint32_t arithmetic_decoder_step(Arithmetic_code a, stream_stats_ptr_t stats, osStream is) {
    uint64_t range = 0, tagGap = 0;
    int32_t k = 0, x = -1, i;
    uint32_t subRange = 0, cumCountX = 0, cumCountX_1 = 0, cumCount = 0;
    
    uint8_t msbU = 0, msbL = 0, E1_E2 = 0, E3 = 0, smsbL = 0, smsbU = 0;
    
	// Again, these are actually constants
	uint32_t msb_shift = a->m - 1;
	uint32_t smsb_shift = a->m - 2;
	uint32_t msb_clear_mask = (1 << msb_shift) - 1;

    range = a->u - a->l + 1;
    tagGap = a->t - a->l + 1;
    
	// @todo figure this out
    subRange = (uint32_t)((tagGap * stats->n - 1) / range);
    while (subRange >= cumCount)
        cumCount += stats->counts[k++];
    x = --k;
  
	cumCountX_1 = 0;
	for (i = 0; i < x; ++i) {
		cumCountX_1 += stats->counts[i];
	}
	cumCountX = cumCountX_1 + stats->counts[x];
    
    a->u = a->l + (uint32_t)((range * cumCountX) / stats->n) - 1;
    a->l = a->l + (uint32_t)((range * cumCountX_1) / stats->n);
    
    // Check the rescaling conditions.
    msbL = a->l >> msb_shift;
    msbU = a->u >> msb_shift;
    
    E1_E2 = (msbL == msbU);
	E3 = 0;
    
    // If E1 or E2 doen't hold, check E3
    if (!E1_E2) {
		smsbL = a->l >> smsb_shift;
		smsbU = a->u >> smsb_shift;
		E3 = (smsbL == 0x01 && smsbU == 0x02);
    }
    
    // While any of E conditions hold
    while (E1_E2 || E3) {
        if (E1_E2) {
			a->l = (a->l & msb_clear_mask) << 1;
			a->u = ((a->u & msb_clear_mask) << 1) + 1;
			a->t = ((a->t & msb_clear_mask) << 1) + stream_read_bit(is);
        }
        else { // E3 is true
			a->l = (a->l << 1) & msb_clear_mask;
			a->u = (((a->u << 1) & msb_clear_mask) | (1 << msb_shift)) + 1;
			a->t = (((a->t & msb_clear_mask) << 1) ^ (1 << msb_shift)) + stream_read_bit(is);
        }
        
		msbL = a->l >> msb_shift;
        msbU = a->u >> msb_shift;
        E1_E2 = (msbL == msbU);
		E3 = 0;
            
        if (!E1_E2) {
			smsbL = a->l >> smsb_shift;
			smsbU = a->u >> smsb_shift;
			E3 = (smsbL == 0x01 && smsbU == 0x02);
        }
    }
    
    return x;
}

uint32_t decoder_last_step(Arithmetic_code a, stream_stats_ptr_t stats) {
    uint64_t range, tagGap, subRange;
    uint32_t k = 0, cumCount = 0, x;
    
    range = a->u - a->l + 1;
    tagGap = a->t - a->l + 1;
    
    subRange = (tagGap * stats->n -1) / range;
    
    while (subRange >= cumCount)
        cumCount += stats->counts[k++];
    
    x = --k;
    
    return x;
}

