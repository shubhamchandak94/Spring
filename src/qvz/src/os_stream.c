#include "qv_compressor.h"

/**
 * Allocates a file stream wrapper for the arithmetic encoder, with a given
 * already opened file handle
 */
struct os_stream_t *alloc_os_stream(FILE *fp, uint8_t in) {
	struct os_stream_t *rtn = (struct os_stream_t *) calloc(1, sizeof(struct os_stream_t));

	rtn->fp = fp;
	rtn->buf = (uint8_t *) calloc(OS_STREAM_BUF_LEN, sizeof(uint8_t));

	if (in) {
		fread(rtn->buf, sizeof(uint8_t), OS_STREAM_BUF_LEN, fp);
	}
	rtn->bufPos = 0;
	rtn->bitPos = 0;
	rtn->written = 0;

	return rtn;
}

/**
 * Deallocate the output stream. Note that this doesn't close the file because
 * this stream doesn't own it
 */
void free_os_stream(struct os_stream_t *os) {
	free(os->buf);
	free(os);
}

/**
 * Reads a single bit from the stream
 */
uint8_t stream_read_bit(struct os_stream_t *os) {
	uint8_t rtn = os->buf[os->bufPos] >> 7;

	os->buf[os->bufPos] = os->buf[os->bufPos] << 1;
	os->bitPos += 1;

	if (os->bitPos == 8) {
		os->bitPos = 0;
		os->bufPos += 1;
		if (os->bufPos == OS_STREAM_BUF_LEN) {
			fread(os->buf, sizeof(uint8_t), OS_STREAM_BUF_LEN, os->fp);
			os->bufPos = 0;
		}
	}

	return rtn;
}

/**
 * Reads a grouping of bits to be interpreted as a single integer, regardless of length
 * Bits are implicitly written in bit endian order ONE AT A TIME elsewhere in the code,
 * so this must read that way too
 */
uint32_t stream_read_bits(struct os_stream_t *os, uint8_t len) {
	uint32_t rtn = 0;
	int8_t bit;

	for (bit = len-1; bit >= 0; --bit) {
		rtn |= stream_read_bit(os) << bit;
	}

	return rtn;
}

/**
 * Writes a single bit to the stream
 */
void stream_write_bit(struct os_stream_t *os, uint8_t bit) {
	bit = (bit & 1);
	os->buf[os->bufPos] |= bit;

	os->bitPos += 1;

	if (os->bitPos == 8) {
		os->bitPos = 0;
		os->bufPos += 1;
		if (os->bufPos == OS_STREAM_BUF_LEN) {
			stream_write_buffer(os);
		}
	}
	else {
		os->buf[os->bufPos] <<= 1;
	}
}

/**
 * Writes a grouping of bits to be interpreted as a single integer and read back the
 * same way. Bits need to be written msb first
 */
void stream_write_bits(struct os_stream_t *os, uint32_t dw, uint8_t len) {
	int8_t bit;

	for (bit = len-1; bit >= 0; --bit) {
		stream_write_bit(os, (uint8_t)(dw >> bit));
	}
}

/**
 * Finishes the current byte in progress and writes the buffer out
 */
void stream_finish_byte(struct os_stream_t *os) {
	os->buf[os->bufPos] <<= (7 - os->bitPos);
	os->bitPos = 0;
	os->bufPos += 1;
	stream_write_buffer(os);
}

/**
 * Writes out the current stream buffer regardless of fill amount
 */
void stream_write_buffer(struct os_stream_t *os) {
	fwrite(os->buf, sizeof(uint8_t), os->bufPos, os->fp);
	memset(os->buf, 0, sizeof(uint8_t)*os->bufPos);
	os->written += os->bufPos;
	os->bufPos = 0;
}
