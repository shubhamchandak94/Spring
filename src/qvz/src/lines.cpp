/**
 * Utility functions for manipulating the data from files, like reading it into
 * memory
 * and converting between the different formats we use
 */

#include "algorithms/SPRING/qvz/include/util.h"

#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>

#include "algorithms/SPRING/qvz/include/lines.h"

namespace spring {
namespace qvz {

/**
 * This reads data from the given file pointer into memory, breaking it into
 * segments
 * of the given number of lines, to ease memory management issues at the cost of
 * some
 * overhead. This assumes that the file consists entirely of quality scores with
 * no
 * other lines in between
 * @param path Path of the file to read
 * @param info Information structure to store in, this must be a valid pointer
 * already
 * @param max_lines Maximum number of lines to read, will override the actual
 * number in the file if >0
 * @todo @xxx This assumes we have only newlines in the file despite some vague
 * attempts to handle \r\n as well
 * @todo @xxx It WILL break the mapping if given a file with \r characters
 * @todo Implement windows analog to mmap to provide the same facility
 */
/*
uint32_t load_file(const char *path, struct quality_file_t *info, uint64_t
max_lines) {
        uint32_t status, block_idx, line_idx;
        char line[READ_LINEBUF_LENGTH];
        FILE *fp;
        int fd;
        struct _stat finfo;
        void *file_mmap;

        // Load metadata into the info structure
        info->path = strdup(path);
        fp = fopen(path, "rt");
        fd = open(path, O_RDONLY);
        if (!fp || fd == -1) {
                return LF_ERROR_NOT_FOUND;
        }

        // Use the first line to figure out how long the file is
        fgets(line, READ_LINEBUF_LENGTH, fp);
        info->columns = strlen(line) - 1;
        if (info->columns > MAX_READS_PER_LINE) {
                fclose(fp);
                return LF_ERROR_TOO_LONG;
        }
        fclose(fp);

        // Figure out how many lines we'll need depending on whether we were
limited or not
        _stat(path, &finfo);
        info->lines = finfo.st_size / ((uint64_t) (info->columns+1));
        if (max_lines > 0 && info->lines > max_lines) {
                info->lines = max_lines;
        }

        status = alloc_blocks(info);
        if (status != LF_ERROR_NONE)
                return status;

        // mmap the file to set up constant pointers indexing it
        file_mmap = mmap(NULL, finfo.st_size, PROT_READ, MAP_SHARED, fd, 0);

        // Process the file
        block_idx = 0;
        line_idx = 0;
        while ((block_idx * MAX_LINES_PER_BLOCK + line_idx) < info->lines) {
                // Setting up mmap indexing assumes we have only one line
ending!
                info->blocks[block_idx].lines[line_idx].m_data = file_mmap +
((uint64_t) (block_idx * MAX_LINES_PER_BLOCK + line_idx)) * (info->columns+1);

                // Increment line/block pointers as necesary
                line_idx += 1;
                if (line_idx == info->blocks[block_idx].count) {
                        line_idx = 0;
                        block_idx += 1;
                }
        }

        return LF_ERROR_NONE;
}
*/
/**
 * Allocate an array of line block pointers and the memory within each block, so
 * that we can
 * use it to store the results of reading the file
 */
// uint32_t alloc_blocks(struct quality_file_t *info) {
//	uint64_t lines_left = info->lines;
//	struct line_block_t *cblock;
//
//	// Figure out how many blocks we'll need to store this file
//	info->block_count = (uint32_t) (info->lines /
//(uint64_t)MAX_LINES_PER_BLOCK);
//	if (info->block_count * MAX_LINES_PER_BLOCK != info->lines) {
//		info->block_count += 1;
//	}
//
//	info->blocks = (struct line_block_t *) calloc(info->block_count,
// sizeof(struct line_block_t));
//	if (!info->blocks) {
//		return LF_ERROR_NO_MEMORY;
//	}
//	cblock = info->blocks;
//
//	while (lines_left > 0) {
//		// Figure out how many lines we'll have in this block
//		if (lines_left > MAX_LINES_PER_BLOCK) {
//			lines_left -= MAX_LINES_PER_BLOCK;
//			cblock->count = MAX_LINES_PER_BLOCK;
//		}
//		else {
//			cblock->count = (uint32_t) lines_left;
//			lines_left = 0;
//		}
//
//		// Allocate array of line info structs for the block
//		cblock->lines = (struct line_t *) calloc(cblock->count,
// sizeof(struct line_t));
//		if (!cblock->lines) {
//			return LF_ERROR_NO_MEMORY;
//		}
//
//		// Advance to the next line block
//		cblock += 1;
//	}
//
//	return LF_ERROR_NONE;
//}

/**
 * Deallocates the memory used to store file information in blocks
 */
void free_blocks(struct quality_file_t *info) {
  // Array of block pointers is a single allocation
  // For each block, array of lines is a single allocations
  //	uint32_t i;

  //	for (i = 0; i < info->block_count; ++i) {
  //		free(info->blocks[i].lines);
  //	}
  free(info->blocks);
}

} // namespace qvz
} // namespace spring
