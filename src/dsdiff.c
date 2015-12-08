/* DSDIFF file support
 *
 * Copyright (c) 2015 Mans Rullgard <mans@mansr.com>
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* File format specification available at
 * http://dsd-guide.com/sites/default/files/white-papers/DSDIFF_1.5_Spec.pdf
 */

#include "sox_i.h"

struct dsdiff {
	uint32_t sample_rate;
	uint16_t num_channels;
	uint64_t data_size;
	uint8_t *buf;
};

#define ID(a, b, c, d) ((a) << 24 | (b) << 16 | (c) << 8 | (d))

static int dff_startread(sox_format_t *ft)
{
	struct dsdiff *dff = ft->priv;
	uint32_t ckid;
	uint32_t cktype;
	uint64_t cksize;
	uint64_t f8size;
	uint32_t fver;
	uint64_t spos, epos;

	if (lsx_readdw(ft, &ckid) || ckid != ID('F', 'R', 'M', '8')) {
		lsx_fail_errno(ft, SOX_EHDR, "FRM8 tag not found");
		return SOX_EHDR;
	}

	if (lsx_readqw(ft, &f8size)) {
		lsx_fail_errno(ft, SOX_EHDR, "error reading chunk size");
		return SOX_EHDR;
	}

	if (lsx_readdw(ft, &cktype) || cktype != ID('D', 'S', 'D', ' ')) {
		lsx_fail_errno(ft, SOX_EHDR, "DSD tag not found");
		return SOX_EHDR;
	}

	do {
		if (lsx_readdw(ft, &ckid) || lsx_readqw(ft, &cksize)) {
			lsx_fail_errno(ft, SOX_EHDR, "read error");
			return SOX_EHDR;
		}

		spos = lsx_tell(ft);

		switch (ckid) {
		case ID('F', 'V', 'E', 'R'):
			if (cksize != 4)
				return SOX_EHDR;
			if (lsx_readdw(ft, &fver))
				return SOX_EHDR;
			if (fver >> 24 != 1) {
				lsx_fail_errno(ft, SOX_EHDR, "unknown version");
				return SOX_EHDR;
			}
			break;

		case ID('P', 'R', 'O', 'P'):
			if (cksize < 4)
				return SOX_EHDR;
			if (lsx_readdw(ft, &cktype))
				return SOX_EHDR;
			if (cktype == ID('S', 'N', 'D', ' '))
				cksize = 4;
			break;

		case ID('F', 'S', ' ', ' '):
			if (cksize < 4)
				return SOX_EHDR;
			if (lsx_readdw(ft, &dff->sample_rate))
				return SOX_EHDR;
			break;

		case ID('C', 'H', 'N', 'L'):
			if (cksize < 4)
				return SOX_EHDR;
			if (lsx_readw(ft, &dff->num_channels))
				return SOX_EHDR;
			break;

		case ID('C', 'M', 'P', 'R'):
			if (cksize < 4)
				return SOX_EHDR;
			if (lsx_readdw(ft, &cktype))
				return SOX_EHDR;
			if (cktype != ID('D', 'S', 'D', ' ')) {
				lsx_fail_errno(ft, SOX_EHDR,
					       "unsupported compression");
				return SOX_EHDR;
			}
			break;

		case ID('D', 'S', 'D', ' '):
			if (cksize < 8)
				return SOX_EHDR;
			dff->data_size = cksize;
			cksize = 0;
			break;
		}

		cksize += cksize & 1;
		epos = lsx_tell(ft);
		if (epos != spos + cksize)
			lsx_seeki(ft, (off_t)(spos + cksize - epos), SEEK_CUR);
	} while (cksize && epos < f8size);

	if (!dff->sample_rate || !dff->num_channels || !dff->data_size) {
		lsx_fail_errno(ft, SOX_EHDR, "invalid file header");
		return SOX_EHDR;
	}

	if (ckid != ID('D', 'S', 'D', ' ')) {
		lsx_fail_errno(ft, SOX_EHDR, "unsupported data type");
		return SOX_EHDR;
	}

	dff->buf = lsx_malloc((size_t)dff->num_channels);
	if (!dff->buf)
		return SOX_ENOMEM;

	ft->signal.rate = dff->sample_rate;
	ft->signal.channels = dff->num_channels;
	ft->signal.precision = 1;
	ft->signal.length = dff->data_size * 8;

	ft->encoding.encoding = SOX_ENCODING_DSD;
	ft->encoding.bits_per_sample = 1;

	return SOX_SUCCESS;
}

static size_t dff_read(sox_format_t *ft, sox_sample_t *buf, size_t len)
{
	struct dsdiff *dff = ft->priv;
	size_t nc = dff->num_channels;
	size_t rsamp = 0;
	unsigned i, j;

	len /= nc;

	while (len >= 8) {
		if (lsx_read_b_buf(ft, dff->buf, nc) < nc)
			return rsamp * nc;

		for (i = 0; i < nc; i++) {
			unsigned d = dff->buf[i];

			for (j = 0; j < 8; j++) {
				buf[i + j * nc] = d & 128 ?
					SOX_SAMPLE_MAX : -SOX_SAMPLE_MAX;
				d <<= 1;
			}
		}

		buf += 8 * nc;
		rsamp += 8;
		len -= 8;
	}

	return rsamp * nc;
}

static int dff_stopread(sox_format_t *ft)
{
	struct dsdiff *dff = ft->priv;

	free(dff->buf);

	return SOX_SUCCESS;
}

LSX_FORMAT_HANDLER(dsdiff)
{
	static char const * const names[] = { "dff", NULL };
	static unsigned const write_encodings[] = { 0 };
	static sox_format_handler_t const handler = {
		SOX_LIB_VERSION_CODE,
		"Direct Stream Digital Interchange File Format (DSDIFF)",
		names, SOX_FILE_BIG_END,
		dff_startread, dff_read, dff_stopread,
		NULL, NULL, NULL,
		NULL, write_encodings, NULL,
		sizeof(struct dsdiff)
	};
	return &handler;
}
