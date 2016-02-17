/* Sigma-Delta modulator
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

/*
 * References:
 *
 * Derk Reefman, Erwin Janssen. 2002.
 * "Signal processing for Direct Stream Digital: A tutorial for
 * digital Sigma Delta modulation and 1-bit digital audio processing"
 * http://www.emmlabs.com/pdf/papers/DerkSigmaDelta.pdf
 *
 * P.J.A. Harpe. 2003.
 * "Trellis-type Sigma Delta Modulators for Super Audio CD applications"
 * http://www.pieterharpe.nl/docs/report_trunc.pdf
 *
 * Richard Schreier. 2000-2011.
 * "Delta Sigma Toolbox"
 * http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox
 */

#define _ISOC11_SOURCE

#include "sox_i.h"
#include "sdm.h"

#define MAX_FILTER_ORDER 8
#define PATH_HASH_SIZE 128
#define PATH_HASH_MASK (PATH_HASH_SIZE - 1)

typedef struct {
  const double  a[MAX_FILTER_ORDER];
  const double  g[MAX_FILTER_ORDER];
  int32_t       order;
  unsigned      freq;
  const char   *name;
  int           trellis_order;
  int           trellis_num;
  int           trellis_lat;
} LSX_ALIGN(32) sdm_filter_t;

typedef struct sdm_state {
  double        state[MAX_FILTER_ORDER];
  double        cost;
  uint32_t      path;
  uint8_t       next;
  uint8_t       hist;
  uint8_t       hist_used;
  struct sdm_state *parent;
  struct sdm_state *path_list;
} LSX_ALIGN(32) sdm_state_t;

typedef struct {
  sdm_state_t   sdm[2 * SDM_TRELLIS_MAX_NUM];
  sdm_state_t  *act[SDM_TRELLIS_MAX_NUM];
} sdm_trellis_t;

struct sdm {
  sdm_trellis_t trellis[2];
  sdm_state_t  *path_hash[PATH_HASH_SIZE];
  uint8_t       hist_free[2 * SDM_TRELLIS_MAX_NUM];
  unsigned      hist_fnum;
  uint32_t      trellis_mask;
  uint32_t      trellis_num;
  uint32_t      trellis_lat;
  unsigned      num_cands;
  unsigned      pos;
  unsigned      pending;
  unsigned      draining;
  unsigned      idx;
  const sdm_filter_t *filter;
  double        prev_y;
  uint64_t      conv_fail;
  uint8_t       hist[2 * SDM_TRELLIS_MAX_NUM][SDM_TRELLIS_MAX_LAT / 8];
};

static sdm_filter_t sdm_filter_safe_64 = {
  {
    7.82475202270363e-01,
    2.94699121204821e-01,
    6.53377434565800e-02,
    8.44052714292759e-03,
    4.35748599576923e-04,
  },
  {
    0, 6.98490600683106e-04,
    0, 1.97734357803445e-03,
  },
  5,
  64 * 44100,
  "safe",
  0, 0, 0,
};

static sdm_filter_t sdm_filter_safe_128 = {
  {
    7.81482424617463e-01,
    2.95962741536345e-01,
    6.66602230172735e-02,
    8.83476316538580e-03,
    5.29668589657667e-04,
  },
  {
    0, 1.74653153894942e-04,
    0, 4.94580504383930e-04,
  },
  5,
  128 * 44100,
  "safe",
  0, 0, 0,
};

static sdm_filter_t sdm_filter_safe_256 = {
  {
    7.81234241629436e-01,
    2.96279053797020e-01,
    6.69905169222529e-02,
    8.93375697799549e-03,
    5.53589139903277e-04,
  },
  {
    0, 4.36651951230006e-05,
    0, 1.23660417994961e-04,
  },
  5,
  256 * 44100,
  "safe",
  0, 0, 0,
};

static sdm_filter_t sdm_filter_fast_64 = {
  {
    9.36269395078196e-01,
    4.20603394734646e-01,
    1.10886565406908e-01,
    1.69435188274605e-02,
    1.07158824957987e-03,
  },
  {
    0, 6.98490600683106e-04,
    0, 1.97734357803445e-03,
  },
  5,
  64 * 44100,
  "fast",
  0, 0, 0,
};

static sdm_filter_t sdm_filter_fast_128 = {
  {
    9.35281560767602e-01,
    4.21733210700821e-01,
    1.12420943220744e-01,
    1.74975342420128e-02,
    1.23024771772653e-03,
  },
  {
    0, 1.74653153894942e-04,
    0, 4.94580504383930e-04,
  },
  5,
  128 * 44100,
  "fast",
  0, 0, 0,
};

static sdm_filter_t sdm_filter_fast_256 = {
  {
    9.35034612630949e-01,
    4.22016068546975e-01,
    1.12804268068047e-01,
    1.76364476199021e-02,
    1.27042574114158e-03,
  },
  {
    0, 4.36651951230006e-05,
    0, 1.23660417994961e-04,
  },
  5,
  256 * 44100,
  "fast",
  0, 0, 0,
};

static sdm_filter_t sdm_filter_hq_64 = {
  {
    7.85192429016760e-01,
    2.98899796337831e-01,
    7.05344069837308e-02,
    1.10501947024035e-02,
    1.07505576489901e-03,
    6.64106274616966e-05,
    -5.74580397150193e-09,
  },
  {
    0, 5.18282192016751e-04,
    0, 1.72954362492419e-03,
    0, 2.83233339830110e-03,
  },
  7,
  64 * 44100,
  "hq",
  0, 0, 0,
};

static sdm_filter_t sdm_filter_hq_128 = {
  {
    7.83148010097334e-01,
    3.01437231238902e-01,
    7.31891646224574e-02,
    1.20314098366875e-02,
    1.32077937193861e-03,
    9.11181979169687e-05,
    2.59895240306562e-06,
  },
  {
    0, 9.92163123766340e-05,
    0, 3.31199917300393e-04,
    0, 5.42540771343282e-04,
  },
  7,
  128 * 44100,
  "hq",
  0, 0, 0,
};

static sdm_filter_t sdm_filter_hq_256 = {
  {
    7.82785077952658e-01,
    3.01888671316811e-01,
    7.36594376027782e-02,
    1.22068270909817e-02,
    1.36572694403914e-03,
    9.57806082134936e-05,
    3.13239368043838e-06,
  },
  {
    0, 2.48046933669715e-05,
    0, 8.28068362972358e-05,
    0, 1.35653594733585e-04,
  },
  7,
  256 * 44100,
  "hq",
  0, 0, 0,
};

static const sdm_filter_t *sdm_filters[] = {
  &sdm_filter_safe_256,
  &sdm_filter_safe_128,
  &sdm_filter_safe_64,
  &sdm_filter_fast_256,
  &sdm_filter_fast_128,
  &sdm_filter_fast_64,
  &sdm_filter_hq_256,
  &sdm_filter_hq_128,
  &sdm_filter_hq_64,
  NULL,
};

static const sdm_filter_t *sdm_find_filter(const char *name, unsigned freq)
{
  int i;

  for (i = 0; sdm_filters[i]; i++)
    if (!name || !strcmp(name, sdm_filters[i]->name))
      if (sdm_filters[i]->freq <= freq)
        return sdm_filters[i];

  return NULL;
}

#include "sdm_x86.h"

#ifndef sdm_filter_calc
static double sdm_filter_calc(const double *s, double *d,
                              const sdm_filter_t *f,
                              double x, double y)
{
  const double *a = f->a;
  const double *g = f->g;
  double v;
  int i;

  d[0] = s[0] - g[0] * s[1] + x - y;
  v = x + a[0] * d[0];

  for (i = 1; i < f->order - 1; i++) {
    d[i] = s[i] + s[i - 1] - g[i] * s[i + 1];
    v += a[i] * d[i];
  }

  d[i] = s[i] + s[i - 1];
  v += a[i] * d[i];

  return v;
}
#endif

#ifndef sdm_filter_calc2
static void sdm_filter_calc2(sdm_state_t *src, sdm_state_t *dst,
                             const sdm_filter_t *f, double x)
{
  const double *a = f->a;
  double v;
  int i;

  v = sdm_filter_calc(src->state, dst[0].state, f, x, 0.0);

  for (i = 0; i < f->order; i++)
    dst[1].state[i] = dst[0].state[i];

  dst[0].state[0] += 1.0;
  dst[1].state[0] -= 1.0;

  dst[0].cost = src->cost + sqr(v + a[0]);
  dst[1].cost = src->cost + sqr(v - a[0]);
}
#endif

static inline unsigned sdm_histbuf_get(sdm_t *p)
{
  return p->hist_free[--p->hist_fnum];
}

static inline void sdm_histbuf_put(sdm_t *p, unsigned h)
{
  p->hist_free[p->hist_fnum++] = h;
}

static inline unsigned get_bit(uint8_t *p, unsigned i)
{
  return (p[i >> 3] >> (i & 7)) & 1;
}

static inline void put_bit(uint8_t *p, unsigned i, unsigned v)
{
  int b = p[i >> 3];
  int s = i & 7;
  b &= ~(1 << s);
  b |= v << s;
  p[i >> 3] = b;
}

static inline unsigned sdm_hist_get(sdm_t *p, unsigned h, unsigned i)
{
  return get_bit(p->hist[h], i);
}

static inline void sdm_hist_put(sdm_t *p, unsigned h, unsigned i, unsigned v)
{
  put_bit(p->hist[h], i, v);
}

static inline void sdm_hist_copy(sdm_t *p, unsigned d, unsigned s)
{
  memcpy(p->hist[d], p->hist[s], (size_t)(p->trellis_lat + 7) / 8);
}

static inline int64_t dbl2int64(double a)
{
  union { double d; int64_t i; } v;
  v.d = a;
  return v.i;
}

static inline int sdm_cmplt(sdm_state_t *a, sdm_state_t *b)
{
  return dbl2int64(a->cost) < dbl2int64(b->cost);
}

static inline int sdm_cmple(sdm_state_t *a, sdm_state_t *b)
{
  return dbl2int64(a->cost) <= dbl2int64(b->cost);
}

static sdm_state_t *sdm_check_path(sdm_t *p, sdm_state_t *s)
{
  unsigned index = s->path & PATH_HASH_MASK;
  sdm_state_t **hash = p->path_hash;
  sdm_state_t *t = hash[index];

  while (t) {
    if (t->path == s->path)
      return t;
    t = t->path_list;
  }

  s->path_list = hash[index];
  hash[index] = s;

  return NULL;
}

static unsigned sdm_sort_cands(sdm_t *p, sdm_trellis_t *st)
{
  sdm_state_t *r, *s, *t;
  sdm_state_t *min;
  unsigned i, j, n;

  for (i = 0; i < 2 * p->num_cands; i++) {
    s = &st->sdm[i];
    p->path_hash[s->path & PATH_HASH_MASK] = NULL;
    if (!i || sdm_cmplt(s, min))
      min = s;
  }

  for (i = 0, n = 0; i < 2 * p->num_cands; i++) {
    s = &st->sdm[i];

    if (s->next != min->next)
      continue;

    if (n == p->trellis_num && sdm_cmple(st->act[n - 1], s))
      continue;

    t = sdm_check_path(p, s);

    if (!t) {
      for (j = n; j > 0; j--) {
        t = st->act[j - 1];
        if (sdm_cmple(t, s))
          break;
        st->act[j] = t;
      }
      if (j < p->trellis_num)
        st->act[j] = s;
      if (n < p->trellis_num)
        n++;
      continue;
    }

    if (sdm_cmple(t, s))
      continue;

    for (j = 0; j < n; j++) {
      r = st->act[j];
      if (sdm_cmple(s, r))
        break;
    }

    st->act[j++] = s;

    while (r != t && j < n) {
      sdm_state_t *u = st->act[j];
      st->act[j] = r;
      r = u;
      j++;
    }
  }

  return n;
}

static inline void sdm_step(sdm_t *p, sdm_state_t *cur, sdm_state_t *next,
                            double x)
{
  const sdm_filter_t *f = p->filter;
  int i;

  sdm_filter_calc2(cur, next, f, x);

  for (i = 0; i < 2; i++) {
    next[i].path = (cur->path << 1 | i) & p->trellis_mask;
    next[i].hist = cur->hist;
    next[i].next = cur->next;
    next[i].parent = cur;
  }
}

static sox_sample_t sdm_sample_trellis(sdm_t *p, double x)
{
  sdm_trellis_t *st_cur = &p->trellis[p->idx];
  sdm_trellis_t *st_next = &p->trellis[p->idx ^ 1];
  double min_cost;
  unsigned new_cands;
  unsigned next_pos;
  unsigned output;
  unsigned i;

  next_pos = p->pos + 1;
  if (next_pos == p->trellis_lat)
    next_pos = 0;

  for (i = 0; i < p->num_cands; i++) {
    sdm_state_t *cur = st_cur->act[i];
    sdm_state_t *next = &st_next->sdm[2 * i];
    sdm_step(p, cur, next, x);
    cur->next = sdm_hist_get(p, cur->hist, next_pos);
    cur->hist_used = 0;
  }

  new_cands = sdm_sort_cands(p, st_next);
  min_cost = st_next->act[0]->cost;
  output = st_next->act[0]->next;

  for (i = 0; i < new_cands; i++) {
    sdm_state_t *s = st_next->act[i];
    if (s->parent->hist_used) {
      unsigned h = sdm_histbuf_get(p);
      sdm_hist_copy(p, h, s->hist);
      s->hist = h;
    } else {
      s->parent->hist_used = 1;
    }

    s->cost -= min_cost;
    s->next = s->parent->next;
    sdm_hist_put(p, s->hist, p->pos, s->path & 1);
  }

  for (i = 0; i < p->num_cands; i++) {
    sdm_state_t *s = st_cur->act[i];
    if (!s->hist_used)
      sdm_histbuf_put(p, s->hist);
  }

  if (new_cands < p->num_cands)
    p->conv_fail++;

  p->num_cands = new_cands;
  p->pos = next_pos;
  p->idx ^= 1;

  return output ? SOX_SAMPLE_MAX : -SOX_SAMPLE_MAX;
}

static sox_sample_t sdm_sample(sdm_t *p, double x)
{
  const sdm_filter_t *f = p->filter;
  double *s0 = p->trellis[0].sdm[p->idx].state;
  double *s1 = p->trellis[0].sdm[p->idx ^ 1].state;
  double y, v;

  v = sdm_filter_calc(s0, s1, f, x, p->prev_y);
  y = signbit(v) ? -1.0 : 1.0;

  p->idx ^= 1;
  p->prev_y = y;

  return y * SOX_SAMPLE_MAX;
}

int sdm_process(sdm_t *p, const sox_sample_t *ibuf, sox_sample_t *obuf,
                size_t *ilen, size_t *olen)
{
  sox_sample_t *out = obuf;
  size_t len = *ilen = min(*ilen, *olen);
  double x;

  if (p->trellis_mask) {
    if (p->pending < p->trellis_lat) {
      size_t pre = min(p->trellis_lat - p->pending, len);
      p->pending += pre;
      len -= pre;
      while (pre--) {
        x = *ibuf++ * (0.5 / SOX_SAMPLE_MAX);
        sdm_sample_trellis(p, x);
      }
    }
    while (len--) {
      x = *ibuf++ * (0.5 / SOX_SAMPLE_MAX);
      *out++ = sdm_sample_trellis(p, x);
    }
  } else {
    while (len--) {
      x = *ibuf++ * (0.5 / SOX_SAMPLE_MAX);
      *out++ = sdm_sample(p, x);
    }
  }

  *olen = out - obuf;

  return SOX_SUCCESS;
}

int sdm_drain(sdm_t *p, sox_sample_t *obuf, size_t *olen)
{
  if (p->trellis_mask) {
    size_t len = *olen = min(p->pending, *olen);

    if (!p->draining && p->pending < p->trellis_lat) {
      unsigned flush = p->trellis_lat - p->pending;
      while (flush--)
        sdm_sample_trellis(p, 0.0);
    }

    p->draining = 1;
    p->pending -= len;

    while (len--)
      *obuf++ = sdm_sample_trellis(p, 0.0);
  } else {
    *olen = 0;
  }

  return SOX_SUCCESS;
}

sdm_t *sdm_init(const char *filter_name,
                unsigned freq,
                unsigned trellis_order,
                unsigned trellis_num,
                unsigned trellis_latency)
{
  sdm_t *p;
  const sdm_filter_t *f;
  sdm_trellis_t *st;
  unsigned i;

  if (trellis_order > SDM_TRELLIS_MAX_ORDER) {
    lsx_fail("trellis order too high (max %d)", SDM_TRELLIS_MAX_ORDER);
    return NULL;
  }

  if (trellis_num > SDM_TRELLIS_MAX_NUM) {
    lsx_fail("trellis size too high (max %d)", SDM_TRELLIS_MAX_NUM);
    return NULL;
  }

  if (trellis_latency > SDM_TRELLIS_MAX_LAT) {
    lsx_fail("trellis latency too high (max %d)", SDM_TRELLIS_MAX_LAT);
    return NULL;
  }

  p = aligned_alloc((size_t)32, sizeof(*p));
  if (!p)
    return NULL;

  memset(p, 0, sizeof(*p));

  p->filter = sdm_find_filter(filter_name, freq);
  if (!p->filter) {
    lsx_fail("invalid filter name `%s'", filter_name);
    return NULL;
  }

  f = p->filter;
  st = &p->trellis[0];

  if (trellis_order || f->trellis_order) {
    if (trellis_order < 1)
      trellis_order = f->trellis_order ? f->trellis_order : 13;

    if (trellis_num)
      p->trellis_num = trellis_num;
    else
      p->trellis_num = f->trellis_num ? f->trellis_num : 8;

    if (trellis_latency)
      p->trellis_lat = trellis_latency;
    else
      p->trellis_lat = f->trellis_lat ? f->trellis_lat : 1024;

    p->trellis_mask = ((uint64_t)1 << trellis_order) - 1;

    for (i = 0; i < 2 * p->trellis_num; i++)
      sdm_histbuf_put(p, i);

    p->num_cands = 1;

    st->sdm[0].hist = sdm_histbuf_get(p);
    st->sdm[0].path = 0;
    st->act[0] = &st->sdm[0];
  }

  return p;
}

void sdm_close(sdm_t *p)
{
  if (p->conv_fail)
    lsx_warn("failed to converge %"PRId64" times", p->conv_fail);

  aligned_free(p);
}

typedef struct sdm_effect {
  sdm_t        *sdm;
  const char   *filter_name;
  uint32_t      trellis_order;
  uint32_t      trellis_num;
  uint32_t      trellis_lat;
} sdm_effect_t;

static int getopts(sox_effect_t *effp, int argc, char **argv)
{
  sdm_effect_t *p = effp->priv;
  lsx_getopt_t optstate;
  int c;

  lsx_getopt_init(argc, argv, "+f:t:n:l:", NULL, lsx_getopt_flag_none,
                  1, &optstate);

  while ((c = lsx_getopt(&optstate)) != -1) switch (c) {
    case 'f': p->filter_name = optstate.arg; break;
    GETOPT_NUMERIC(optstate, 't', trellis_order, 3, SDM_TRELLIS_MAX_ORDER)
    GETOPT_NUMERIC(optstate, 'n', trellis_num, 4, SDM_TRELLIS_MAX_NUM)
    GETOPT_NUMERIC(optstate, 'l', trellis_lat, 100, SDM_TRELLIS_MAX_LAT)
    default: lsx_fail("invalid option `-%c'", optstate.opt); return lsx_usage(effp);
  }

  return argc != optstate.ind ? lsx_usage(effp) : SOX_SUCCESS;
}

static int start(sox_effect_t *effp)
{
  sdm_effect_t *p = effp->priv;

  p->sdm = sdm_init(p->filter_name, effp->in_signal.rate,
                    p->trellis_order, p->trellis_num, p->trellis_lat);
  if (!p->sdm)
    return SOX_EOF;

  effp->out_signal.precision = 1;

  return SOX_SUCCESS;
}

static int flow(sox_effect_t *effp, const sox_sample_t *ibuf,
                sox_sample_t *obuf, size_t *isamp, size_t *osamp)
{
  sdm_effect_t *p = effp->priv;
  return sdm_process(p->sdm, ibuf, obuf, isamp, osamp);
}

static int drain(sox_effect_t *effp, sox_sample_t *obuf, size_t *osamp)
{
  sdm_effect_t *p = effp->priv;
  return sdm_drain(p->sdm, obuf, osamp);
}

static int stop(sox_effect_t *effp)
{
  sdm_effect_t *p = effp->priv;
  sdm_close(p->sdm);
  return SOX_SUCCESS;
}

const sox_effect_handler_t *lsx_sdm_effect_fn(void)
{
  static sox_effect_handler_t handler = {
    "sdm", "[-f filter] [-t order] [-n num] [-l latency]"
    "\n  -f       Set filter to one of: safe, fast, hq"
    "\n           Advanced options:"
    "\n  -t       Override trellis order"
    "\n  -n       Override number of trellis paths"
    "\n  -l       Override trellis latency",
    SOX_EFF_PREC, getopts, start, flow, drain, stop, 0, sizeof(sdm_effect_t),
  };
  return &handler;
}
