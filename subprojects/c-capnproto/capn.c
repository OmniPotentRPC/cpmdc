/* vim: set sw=8 ts=8 sts=8 noet: */
#include "capn.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>

#define STRUCT_PTR 0
#define LIST_PTR 1
#define FAR_PTR 2
#define DOUBLE_PTR 6

#define VOID_LIST 0
#define BIT_1_LIST 1
#define BYTE_1_LIST 2
#define BYTE_2_LIST 3
#define BYTE_4_LIST 4
#define BYTE_8_LIST 5
#define PTR_LIST 6
#define COMPOSITE_LIST 7

#define U64(val) ((uint64_t) (val))
#define I64(val) ((int64_t) (val))
#define U32(val) ((uint32_t) (val))
#define I32(val) ((int32_t) (val))
#define U16(val) ((uint16_t) (val))
#define I16(val) ((int16_t) (val))

#ifndef min
static int min(int a, int b) { return (a < b) ? a : b; }
#endif

#ifdef BYTE_ORDER
#define CAPN_LITTLE (BYTE_ORDER == LITTLE_ENDIAN)
#elif defined(__BYTE_ORDER)
#define CAPN_LITTLE (__BYTE_ORDER == __LITTLE_ENDIAN)
#else
#define CAPN_LITTLE 0
#endif

struct capn_tree *capn_tree_insert(struct capn_tree *root, struct capn_tree *n) {
	n->red = 1;
	n->link[0] = n->link[1] = NULL;

	for (;;) {
		/* parent, uncle, grandparent, great grandparent link */
		struct capn_tree *p, *u, *g, **gglink;
		int dir;

		/* Case 1: N is root */
		p = n->parent;
		if (!p) {
			n->red = 0;
			root = n;
			break;
		}

		/* Case 2: p is black */
		if (!p->red) {
			break;
		}

		g = p->parent;
		dir = (p == g->link[1]);

		/* Case 3: P and U are red, switch g to red, but must
		 * loop as G could be root or have a red parent
		 *     g    to   G
		 *    / \       / \
		 *   P   U     p   u
		 *  /         /
		 * N         N
		 */
		u = g->link[!dir];
		if (u != NULL && u->red) {
			p->red = 0;
			u->red = 0;
			g->red = 1;
			n = g;
			continue;
		}

		if (!g->parent) {
			gglink = &root;
		} else if (g->parent->link[1] == g) {
			gglink = &g->parent->link[1];
		} else {
			gglink = &g->parent->link[0];
		}

		if (dir != (n == p->link[1])) {
			/* Case 4: rotate on P, then on g
			 * here dir is /
			 *     g    to   g   to   n
			 *    / \       / \      / \
			 *   P   u     N   u    P   G
			 *  / \       / \      /|  / \
			 * 1   N     P   3    1 2 3   u
			 *    / \   / \
			 *   2   3 1   2
			 */
			struct capn_tree *two = n->link[dir];
			struct capn_tree *three = n->link[!dir];
			p->link[!dir] = two;
			g->link[dir] = three;
			n->link[dir] = p;
			n->link[!dir] = g;
			*gglink = n;
			n->parent = g->parent;
			p->parent = n;
			g->parent = n;
			if (two)
				two->parent = p;
			if (three)
				three->parent = g;
			n->red = 0;
			g->red = 1;
		} else {
			/* Case 5: rotate on g
			 * here dir is /
			 *       g   to   p
			 *      / \      / \
			 *     P   u    N   G
			 *    / \      /|  / \
			 *   N   3    1 2 3   u
			 *  / \
			 * 1   2
			 */
			struct capn_tree *three = p->link[!dir];
			g->link[dir] = three;
			p->link[!dir] = g;
			*gglink = p;
			p->parent = g->parent;
			g->parent = p;
			if (three)
				three->parent = g;
			g->red = 1;
			p->red = 0;
		}

		break;
	}

	return root;
}

void capn_append_segment(struct capn *c, struct capn_segment *s) {
	/* C11 6.5.3.2: a null message is not a segment list. */
	if (!c)
		return;
	/* C11 6.5.3.2: a null segment is not a list node. */
	if (!s)
		return;
	s->id = c->segnum++;
	s->capn = c;
	s->next = NULL;

	if (c->lastseg) {
		c->lastseg->next = s;
		c->lastseg->hdr.link[1] = &s->hdr;
		s->hdr.parent = &c->lastseg->hdr;
	} else {
		c->seglist = s;
		s->hdr.parent = NULL;
	}

	c->lastseg = s;
	c->segtree = capn_tree_insert(c->segtree, &s->hdr);
}

/* INT32-C: a segment length plus a byte count must fit in a signed int. */
static int seg_room(const struct capn_segment *s, int sz) {
	if (!s || sz < 0 || s->len < 0 || s->cap < 0)
		return 0;
	if ((uint64_t)s->len > (uint64_t)INT_MAX - (uint64_t)sz)
		return 0;
	return (uint64_t)s->len + (uint64_t)sz <= (uint64_t)s->cap;
}

static char *new_data(struct capn *c, int sz, struct capn_segment **ps) {
	struct capn_segment *s;

	/* C11 6.5.3.2: a null out-pointer is not a store. */
	if (!ps)
		return NULL;
	/* C11 6.5.3.2: a null message is not a segment list. */
	if (!c) {
		*ps = NULL;
		return NULL;
	}

	/* find a segment with sufficient data */
	for (s = c->seglist; s != NULL; s = s->next) {
		if (seg_room(s, sz)) {
			goto end;
		}
	}

	s = c->create ? c->create(c->user, c->segnum, sz) : NULL;
	if (!s || !seg_room(s, sz)) {
		*ps = NULL;
		return NULL;
	}

	capn_append_segment(c, s);
end:
	*ps = s;
	/* C11 6.5.6: do not add to a null data pointer. */
	if (!s->data)
		return NULL;
	s->len += sz;
	return s->data + s->len - sz;
}

static struct capn_segment *lookup_segment(struct capn* c, struct capn_segment *s, uint32_t id) {
	struct capn_tree **x;
	struct capn_segment *y;

	if (s && s->id == id)
		return s;
	if (!c)
		return NULL;

	if (id < c->segnum) {
		x = &c->segtree;
		y = NULL;
		while (*x) {
			y = (struct capn_segment*) *x;
			if (id == y->id) {
				return y;
			} else if (id < y->id) {
				x = &y->hdr.link[0];
			} else {
				x = &y->hdr.link[1];
			}
		}
	}

	s = c->lookup ? c->lookup(c->user, id) : NULL;
	if (!s)
		return NULL;

	if (id < c->segnum) {
		s->id = id;
		s->capn = c;
		s->next = c->seglist;
		c->seglist = s;
		s->hdr.parent = &y->hdr;
		*x = &s->hdr;
		c->segtree = capn_tree_insert(c->segtree, &s->hdr);
	} else {
		c->segnum = id;
		capn_append_segment(c, s);
	}

	return s;
}

/* INT30-C: a far word offset times eight must not wrap. */
static int far_off(uint64_t val, uint64_t need, const struct capn_segment *s, uint64_t *off) {
	uint64_t bytes = (uint64_t)(U32(val) >> 3) * 8ull;
	if (!s || s->len < 0 || bytes > (uint64_t)s->len || need > (uint64_t)s->len - bytes)
		return -1;
	*off = bytes;
	return 0;
}

static uint64_t lookup_double(struct capn_segment **s, char **d, uint64_t val) {
	uint64_t far, tag, off;
	char *p;

	/* C11 6.5.3.2: a null segment cursor is not a lookup. */
	if (!s || !*s)
		return 0;

	if ((*s = lookup_segment((*s)->capn, *s, U32(val >> 32))) == NULL) {
		return 0;
	}

	if (far_off(val, 16, *s, &off))
		return 0;
	/* C11 6.5.6: do not add to a null data pointer. */
	if (!(*s)->data)
		return 0;
	p = (*s)->data + off;

	far = capn_flip64(*(uint64_t*) p);
	tag = capn_flip64(*(uint64_t*) (p+8));

	/* the far tag should not be another double, and the tag
	 * should be struct/list and have no offset */
	if ((far&7) != FAR_PTR || U32(tag) > LIST_PTR) {
		return 0;
	}

	if ((*s = lookup_segment((*s)->capn, *s, U32(far >> 32))) == NULL) {
		return 0;
	}

	/* C11 6.5.6: a pointer before the segment is undefined.
	   The caller adds eight bytes. One less word lands on
	   the same byte from the segment start. */
	if (!(*s)->data)
		return 0;
	*d = (*s)->data;
	{
		uint32_t words = U32(far) >> 3;
		uint32_t field = (uint32_t)((int64_t)words - 1) << 2;
		return U64(field) | tag;
	}
}

static uint64_t lookup_far(struct capn_segment **s, char **d, uint64_t val) {
	uint64_t off;

	/* C11 6.5.3.2: a null segment cursor is not a lookup. */
	if (!s || !*s)
		return 0;

	if ((*s = lookup_segment((*s)->capn, *s, U32(val >> 32))) == NULL) {
		return 0;
	}

	if (far_off(val, 8, *s, &off))
		return 0;
	/* C11 6.5.6: do not add to a null data pointer. */
	if (!(*s)->data)
		return 0;

	*d = (*s)->data + off;
	return capn_flip64(*(uint64_t*)*d);
}

/* C11 6.5.6: a pointer difference is defined only inside one array.
   len is the segment bound. An address inside it keeps the byte offset. */
static int seg_data_off(const struct capn_segment *s, const char *p, uint64_t *off) {
	uintptr_t base, at;

	if (!s || !s->data || !p || s->len < 0)
		return -1;
	base = (uintptr_t) s->data;
	at = (uintptr_t) p;
	if (at < base || (uint64_t) (at - base) > (uint64_t) s->len)
		return -1;
	*off = (uint64_t) (at - base);
	return 0;
}

/* C11 6.5.3.2: a pointer word is eight bytes inside the segment. */
static int load_word(const struct capn_segment *s, const char *d, uint64_t *val) {
	uint64_t off;
	if (seg_data_off(s, d, &off) || (uint64_t)s->len - off < 8ull)
		return -1;
	*val = capn_flip64(*(uint64_t *)d);
	return 0;
}

static char *struct_ptr(struct capn_segment *s, char *d, int minsz) {
	uint64_t val;
	uint16_t datasz;

	if (load_word(s, d, &val))
		return NULL;

	switch (val&7) {
	case FAR_PTR:
		val = lookup_far(&s, &d, val);
		break;
	case DOUBLE_PTR:
		val = lookup_double(&s, &d, val);
		break;
	}

	datasz = U16(val >> 32);
	/* INT32-C: the word offset is a signed 30-bit field. Scale it in int64_t. */
	{
		uint32_t raw = U32(val) >> 2;
		int64_t words = (raw & 0x20000000u)
			? (int64_t)raw - (int64_t)0x40000000
			: (int64_t)raw;
		int64_t delta = words * 8 + 8;
		uint64_t base;
		if (seg_data_off(s, d, &base))
			return NULL;
		if (delta < 0) {
			if ((uint64_t)(-delta) > base)
				return NULL;
			d = s->data + (base - (uint64_t)(-delta));
		} else if ((uint64_t)delta > (uint64_t)s->len - base) {
			return NULL;
		} else {
			d = s->data + (base + (uint64_t)delta);
		}
	}

	if (val != 0 && (val&3) != STRUCT_PTR && datasz >= minsz && s->data <= d && d < s->data + s->len) {
		return d;
	}

	return NULL;
}

/* INT32-C: a decoded span must fit in the bytes that remain. */
static int list_end(const struct capn_segment *s, const char *d, uint64_t nbytes, char **end) {
	uint64_t off;
	if (seg_data_off(s, d, &off))
		return -1;
	if (off > (uint64_t)s->len || nbytes > (uint64_t)s->len - off)
		return -1;
	*end = (char *)d + nbytes;
	return 0;
}

static capn_ptr read_ptr(struct capn_segment *s, char *d) {
	capn_ptr ret = {CAPN_NULL};
	uint64_t val;
	char *e;

	if (load_word(s, d, &val))
		goto err;

	switch (val&7) {
	case FAR_PTR:
		val = lookup_far(&s, &d, val);
		ret.has_ptr_tag = (U32(val) >> 2) == 0;
		break;
	case DOUBLE_PTR:
		val = lookup_double(&s, &d, val);
		break;
	}

	/* INT32-C: the word offset is a signed 30-bit field. Scale it in int64_t. */
	{
		uint32_t raw = U32(val) >> 2;
		int64_t words = (raw & 0x20000000u)
			? (int64_t)raw - (int64_t)0x40000000
			: (int64_t)raw;
		int64_t delta = words * 8 + 8;
		uint64_t base;
		if (seg_data_off(s, d, &base))
			goto err;
		if (delta < 0) {
			if ((uint64_t)(-delta) > base)
				goto err;
			d = s->data + (base - (uint64_t)(-delta));
		} else if ((uint64_t)delta > (uint64_t)s->len - base) {
			goto err;
		} else {
			d = s->data + (base + (uint64_t)delta);
		}
	}

	switch (val & 3) {
	case STRUCT_PTR:
		ret.type = val ? CAPN_STRUCT : CAPN_NULL;
		goto struct_common;

	struct_common:
		ret.datasz = U32(U16(val >> 32)) * 8;
		ret.ptrs = U32(U16(val >> 48));
		if (list_end(s, d, (uint64_t)ret.datasz + 8ull * (uint64_t)ret.ptrs, &e))
			goto err;
		break;

	case LIST_PTR:
		ret.type = CAPN_LIST;
		ret.len = val >> 35;

		switch ((val >> 32) & 7) {
		case VOID_LIST:
			if (list_end(s, d, 0, &e))
				goto err;
			break;
		case BIT_1_LIST:
			ret.type = CAPN_BIT_LIST;
			if (ret.len < 0 || ret.len > INT_MAX - 7)
				goto err;
			ret.datasz = (ret.len + 7) / 8;
			if (list_end(s, d, (uint64_t)ret.datasz, &e))
				goto err;
			break;
		case BYTE_1_LIST:
			ret.datasz = 1;
			if (ret.len < 0 || list_end(s, d, (uint64_t)ret.len, &e))
				goto err;
			break;
		case BYTE_2_LIST:
			ret.datasz = 2;
			if (ret.len < 0 || list_end(s, d, (uint64_t)ret.len * 2ull, &e))
				goto err;
			break;
		case BYTE_4_LIST:
			ret.datasz = 4;
			if (ret.len < 0 || list_end(s, d, (uint64_t)ret.len * 4ull, &e))
				goto err;
			break;
		case BYTE_8_LIST:
			ret.datasz = 8;
			if (ret.len < 0 || list_end(s, d, (uint64_t)ret.len * 8ull, &e))
				goto err;
			break;
		case PTR_LIST:
			ret.type = CAPN_PTR_LIST;
			if (ret.len < 0 || list_end(s, d, (uint64_t)ret.len * 8ull, &e))
				goto err;
			break;
		case COMPOSITE_LIST: {
			uint64_t declared, elem, nbytes;
			if (ret.len < 0)
				goto err;
			declared = (uint64_t)ret.len * 8ull;
			if (list_end(s, d, 8, &e))
				goto err;

			val = capn_flip64(*(uint64_t*) d);

			d += 8;
			ret.datasz = U32(U16(val >> 32)) * 8;
			ret.ptrs = U32(U16(val >> 48));
			ret.len = U32(val) >> 2;
			ret.is_composite_list = 1;

			elem = (uint64_t)ret.datasz + 8ull * (uint64_t)ret.ptrs;
			if (ret.len < 0 || (ret.len > 0 && elem > UINT64_MAX / (uint64_t)ret.len))
				goto err;
			nbytes = elem * (uint64_t)ret.len;
			if (nbytes != declared || list_end(s, d, nbytes, &e))
				goto err;
			break;
		}
		}
		break;

	default:
		goto err;
	}

	{
		uint64_t end_off;
		if (seg_data_off(s, e, &end_off))
			goto err;
	}

	ret.data = d;
	ret.seg = s;
	return ret;
err:
	memset(&ret, 0, sizeof(ret));
	return ret;
}

void capn_resolve(capn_ptr *p) {
	/* C11 6.5.3.2: a null pointer is not a far pointer. */
	if (!p)
		return;
	if (p->type == CAPN_FAR_POINTER) {
		*p = read_ptr(p->seg, p->data);
	}
}

/* TODO: should this handle CAPN_BIT_LIST? */
/* INT32-C: a member byte offset must fit in the bytes that remain. */
static int ptr_at(const capn_ptr *p, uint64_t bytes, char **out) {
	uint64_t base;
	/* C11 6.5.8: relational compare is defined only inside one array. */
	if (!p->seg || seg_data_off(p->seg, p->data, &base))
		return -1;
	if (bytes > (uint64_t)p->seg->len - base)
		return -1;
	*out = p->seg->data + base + bytes;
	return 0;
}

capn_ptr capn_getp(capn_ptr p, int off, int resolve) {
	capn_ptr ret = {CAPN_FAR_POINTER};
	ret.seg = p.seg;

	capn_resolve(&p);

	switch (p.type) {
	case CAPN_LIST:
		/* INT32-C: the member offset is an index times the element size. */
		if (off < 0 || off >= p.len)
			goto err;
		{
			uint64_t elem = (uint64_t)p.datasz + 8ull * (uint64_t)p.ptrs;
			uint64_t bytes;
			capn_ptr inner = {CAPN_STRUCT};
			if (off > 0 && elem > UINT64_MAX / (uint64_t)off)
				goto err;
			bytes = (uint64_t)off * elem;
			if (ptr_at(&p, bytes, &inner.data))
				goto err;
			inner.is_list_member = 1;
			inner.seg = p.seg;
			inner.datasz = p.datasz;
			inner.ptrs = p.ptrs;
			return inner;
		}

	case CAPN_STRUCT:
		if (off < 0 || off >= p.ptrs)
			goto err;
		{
			uint64_t bytes = (uint64_t)p.datasz + 8ull * (uint64_t)off;
			if (ptr_at(&p, bytes, &ret.data))
				goto err;
		}
		break;

	case CAPN_PTR_LIST:
		if (off < 0 || off >= p.len)
			goto err;
		{
			uint64_t bytes = (uint64_t)off * 8ull;
			if (ptr_at(&p, bytes, &ret.data))
				goto err;
		}
		break;

	default:
		goto err;
	}

	if (resolve) {
		ret = read_ptr(ret.seg, ret.data);
	}

	return ret;

err:
	memset(&p, 0, sizeof(p));
	return p;
}

static void write_ptr_tag(char *d, capn_ptr p, int off) {
	/* INT34-C: off/8 is signed. Shift the unsigned bits. */
	uint64_t val = U64((uint32_t)I32(off/8) << 2);

	switch (p.type) {
	case CAPN_STRUCT:
		val |= STRUCT_PTR | (U64(p.datasz/8) << 32) | (U64(p.ptrs) << 48);
		break;

	case CAPN_LIST:
		if (p.is_composite_list) {
			/* INT32-C: the word count is a signed length times words per element. */
			uint64_t words = 0;
			if (p.len > 0 && p.datasz >= 0 && p.ptrs >= 0)
				words = (uint64_t)p.len * ((uint64_t)(p.datasz / 8) + (uint64_t)p.ptrs);
			val |= LIST_PTR | (U64(COMPOSITE_LIST) << 32) | (words << 35);
		} else {
			val |= LIST_PTR | (U64(p.len) << 35);

			switch (p.datasz) {
			case 8:
				val |= (U64(BYTE_8_LIST) << 32);
				break;
			case 4:
				val |= (U64(BYTE_4_LIST) << 32);
				break;
			case 2:
				val |= (U64(BYTE_2_LIST) << 32);
				break;
			case 1:
				val |= (U64(BYTE_1_LIST) << 32);
				break;
			case 0:
				val |= (U64(VOID_LIST) << 32);
				break;
			}
		}
		break;

	case CAPN_BIT_LIST:
		val |= LIST_PTR | (U64(BIT_1_LIST) << 32) | (U64(p.len) << 35);
		break;

	case CAPN_PTR_LIST:
		val |= LIST_PTR | (U64(PTR_LIST) << 32) | (U64(p.len) << 35);
		break;

	default:
		val = 0;
		break;
	}

	/* C11 6.5.3.2: a null destination is not a store. */
	if (!d)
		return;
	*(uint64_t*) d = capn_flip64(val);
}

/* C11 6.5.6 / INT31-C: a tag offset is a pointer difference narrowed to int.
   Addresses through one past cap stay comparable. Outside that, fail. */
static int tag_off_int(const struct capn_segment *s, const char *pdata, const char *at, int *off) {
	uintptr_t base, p_at, d_at;
	uint64_t po, ao;
	int64_t delta;

	if (!s || !s->data || !pdata || !at || !off || s->cap < 0)
		return -1;
	base = (uintptr_t) s->data;
	p_at = (uintptr_t) pdata;
	d_at = (uintptr_t) at;
	if (p_at < base || (uint64_t) (p_at - base) > (uint64_t) s->cap)
		return -1;
	if (d_at < base || (uint64_t) (d_at - base) > (uint64_t) s->cap)
		return -1;
	po = (uint64_t) (p_at - base);
	ao = (uint64_t) (d_at - base);
	delta = (int64_t) po - (int64_t) ao - 8;
	if (delta < (int64_t) INT_MIN || delta > (int64_t) INT_MAX)
		return -1;
	*off = (int) delta;
	return 0;
}

/* C11 6.5.6: a pointer difference is defined only inside one array.
   The segment bound is cap. An address inside it keeps the byte offset. */
static int far_byte_off(const struct capn_segment *s, const char *tgt, uint64_t *off) {
	uintptr_t base, at;

	if (!s || !s->data || !tgt || s->cap <= 0)
		return -1;
	base = (uintptr_t) s->data;
	at = (uintptr_t) tgt;
	if (at < base || (uint64_t) (at - base) >= (uint64_t) s->cap)
		return -1;
	*off = (uint64_t) (at - base);
	return 0;
}

static void write_far_ptr(char *d, struct capn_segment *s, char *tgt) {
	uint64_t off = 0;
	uint64_t val = 0;

	/* C11 6.5.3.2: a null destination is not a store. */
	if (!d)
		return;
	if (!far_byte_off(s, tgt, &off))
		val = FAR_PTR | off | (U64(s->id) << 32);
	*(uint64_t*) d = capn_flip64(val);
}

static void write_double_far(char *d, struct capn_segment *s, char *tgt) {
	uint64_t off = 0;
	uint64_t val = 0;

	/* C11 6.5.3.2: a null destination is not a store. */
	if (!d)
		return;
	if (!far_byte_off(s, tgt, &off))
		val = DOUBLE_PTR | off | (U64(s->id) << 32);
	*(uint64_t*) d = capn_flip64(val);
}

#define NEED_TO_COPY 1

static int write_ptr(struct capn_segment *s, char *d, capn_ptr p) {
	/* note p.seg can be NULL if its a ptr to static data */
	/* C11 6.5.6: do not step a null pointer, or step before the segment. */
	char *pdata = p.data;

	if (p.is_composite_list && p.data && p.seg && p.seg->data && p.seg->cap >= 8) {
		uintptr_t base = (uintptr_t) p.seg->data;
		uintptr_t at = (uintptr_t) p.data;
		if (at >= base + 8u && (uint64_t) (at - base) <= (uint64_t) p.seg->cap)
			pdata = p.data - 8;
	}

	if (p.type == CAPN_NULL || (p.type == CAPN_STRUCT && p.datasz == 0 && p.ptrs == 0)) {
		write_ptr_tag(d, p, 0);
		return 0;

	} else if (!p.seg || p.seg->capn != s->capn || p.is_list_member) {
		return NEED_TO_COPY;

	} else if (p.seg == s) {
		int off = 0;
		if (tag_off_int(s, pdata, d, &off))
			return -1;
		write_ptr_tag(d, p, off);
		return 0;

	} else if (p.has_ptr_tag) {
		/* By lucky chance, the data has a tag in front
		 * of it. This happens when new_object had to move
		 * the data to a new segment. */
		/* C11 6.5.6: form the tag address only inside the segment. */
		char *tag = NULL;
		if (pdata && p.seg && p.seg->data && p.seg->cap >= 8) {
			uintptr_t base = (uintptr_t) p.seg->data;
			uintptr_t at = (uintptr_t) pdata;
			if (at >= base + 8u && (uint64_t) (at - base) <= (uint64_t) p.seg->cap)
				tag = pdata - 8;
		}
		if (!tag)
			return -1;
		write_far_ptr(d, p.seg, tag);
		return 0;

	} else if (seg_room(p.seg, 8)) {
		/* The target segment has enough room for tag */
		char *t;
		int off = 0;
		/* C11 6.5.6: do not add to a null data pointer. */
		if (!p.seg->data)
			return -1;
		t = p.seg->data + p.seg->len;
		if (tag_off_int(p.seg, pdata, t, &off))
			return -1;
		write_ptr_tag(t, p, off);
		write_far_ptr(d, p.seg, t);
		p.seg->len += 8;
		return 0;

	} else {
		/* have to allocate room for a double far
		 * pointer */
		char *t;

		if (seg_room(s, 16)) {
			/* Try and allocate in the src segment
			 * first. This should improve lookup on
			 * read. */
			/* C11 6.5.6: do not add to a null data pointer. */
			if (!s->data)
				return -1;
			t = s->data + s->len;
			s->len += 16;
		} else {
			t = new_data(s->capn, 16, &s);
			if (!t) return -1;
		}

		write_far_ptr(t, p.seg, pdata);
		write_ptr_tag(t+8, p, 0);
		write_double_far(d, s, t);
		return 0;
	}
}

struct copy {
	struct capn_tree hdr;
	struct capn_ptr to, from;
	char *fbegin, *fend;
};

static capn_ptr new_clone(struct capn_segment *s, capn_ptr p) {
	switch (p.type) {
	case CAPN_STRUCT:
		return capn_new_struct(s, p.datasz, p.ptrs);
	case CAPN_PTR_LIST:
		return capn_new_ptr_list(s, p.len);
	case CAPN_BIT_LIST:
		return capn_new_list1(s, p.len).p;
	case CAPN_LIST:
		return capn_new_list(s, p.len, p.datasz, p.ptrs);
	default:
		return p;
	}
}

static int is_ptr_equal(const struct capn_ptr *a, const struct capn_ptr *b) {
	return a->data == b->data
		&& a->type == b->type
		&& a->len == b->len
		&& a->datasz == b->datasz
		&& a->ptrs == b->ptrs;
}

static int data_size(struct capn_ptr p) {
	uint64_t n;
	switch (p.type) {
	case CAPN_BIT_LIST:
		return (int)p.datasz;
	case CAPN_PTR_LIST:
		/* INT32-C: a pointer list is eight bytes times the length. */
		if (p.len < 0 || (uint64_t)p.len > (uint64_t)INT_MAX / 8ull)
			return -1;
		return p.len * 8;
	case CAPN_STRUCT:
		n = (uint64_t)p.datasz + 8ull * (uint64_t)p.ptrs;
		if (n > (uint64_t)INT_MAX)
			return -1;
		return (int)n;
	case CAPN_LIST:
		n = (uint64_t)p.datasz + 8ull * (uint64_t)p.ptrs;
		if (p.len < 0)
			return -1;
		if (p.len > 0 && n > (uint64_t)INT_MAX / (uint64_t)p.len)
			return -1;
		n = (uint64_t)p.len * n + (p.is_composite_list ? 8ull : 0ull);
		if (n > (uint64_t)INT_MAX)
			return -1;
		return (int)n;
	default:
		return 0;
	}
}

static int copy_ptr(struct capn_segment *seg, char *data, struct capn_ptr *t, struct capn_ptr *f, int *dep) {
	struct capn *c = seg->capn;
	struct copy *cp = NULL;
	struct capn_tree **xcp;
	int span = data_size(*f);
	char *fbegin;
	char *fend;
	int zero_sized;
	if (span < 0)
		return -1;
	/* C11 6.5.6: form a composite tag only inside the segment. */
	if (f->is_composite_list) {
		uintptr_t base, at;
		if (!f->data || !f->seg || !f->seg->data || f->seg->cap < 8)
			return -1;
		base = (uintptr_t) f->seg->data;
		at = (uintptr_t) f->data;
		if (at < base + 8u || (uint64_t) (at - base) > (uint64_t) f->seg->cap)
			return -1;
		fbegin = f->data - 8;
	} else if (!f->data) {
		if (span > 0)
			return -1;
		fbegin = NULL;
	} else {
		fbegin = f->data;
	}
	if (!fbegin) {
		fend = NULL;
		zero_sized = 1;
	} else {
		fend = fbegin + span;
		zero_sized = (fend == fbegin);
	}

	/* We always copy list members as it would otherwise be an
	 * overlapped pointer (the data is owned by the enclosing list).
	 * We do not bother with the overlapped lookup for zero sized
	 * structures/lists as they never overlap. Nor do we add them to
	 * the copy list as there is no data to be shared by multiple
	 * pointers.
	 */

	xcp = &c->copy;
	while (*xcp && !zero_sized) {
		cp = (struct copy*) *xcp;
		/* C11 6.5.8: order addresses as integers, not as pointers. */
		if ((uintptr_t) fend <= (uintptr_t) cp->fbegin) {
			xcp = &cp->hdr.link[0];
		} else if ((uintptr_t) cp->fend <= (uintptr_t) fbegin) {
			xcp = &cp->hdr.link[1];
		} else if (is_ptr_equal(f, &cp->from)) {
			/* we already have a copy so just point to that */
			return write_ptr(seg, data, cp->to);
		} else {
			/* pointer to overlapped data */
			return -1;
		}
	}

	/* no copy found - have to create a new copy */
	*t = new_clone(seg, *f);

	if (write_ptr(seg, data, *t))
		return -1;

	/* add the copy to the copy tree so we can look for overlapping
	 * source pointers and handle recursive structures */
	if (!zero_sized) {
		struct copy *n;
		struct capn_segment *cs = c->copylist;

		/* INT31-C: a negative length is not a size_t room check. */
		if (cs && (cs->len < 0 || cs->cap < 0))
			return -1;
		/* need to allocate a struct copy */
		if (!cs || (uint64_t) cs->len > (uint64_t) INT_MAX - (uint64_t) sizeof(*n)
		    || (uint64_t) cs->len + (uint64_t) sizeof(*n) > (uint64_t) cs->cap) {
			cs = c->create_local ? c->create_local(c->user, sizeof(*n)) : NULL;
			if (!cs) {
				/* can't allocate a copy structure */
				return -1;
			}
			cs->next = c->copylist;
			c->copylist = cs;
		}

		/* C11 6.5.6: do not add to a null data pointer. */
		if (!cs->data)
			return -1;
		n = (struct copy*) (cs->data + cs->len);
		cs->len += sizeof(*n);

		n->from = *f;
		n->to = *t;
		n->fbegin = fbegin;
		n->fend = fend;

		*xcp = &n->hdr;
		n->hdr.parent = &cp->hdr;

		c->copy = capn_tree_insert(c->copy, &n->hdr);
	}

	/* minimize the number of types the main copy routine has to
	 * deal with to just CAPN_LIST and CAPN_PTR_LIST. ptr list only
	 * needs t->type, t->len, t->data, t->seg, f->data, f->seg to
	 * be valid */
	switch (t->type) {
	case CAPN_STRUCT:
		if (t->datasz) {
			/* C11 7.24.1p2: a positive count needs both pointers. */
			if (!t->data || !f->data)
				return -1;
			memcpy(t->data, f->data, t->datasz);
			t->data += t->datasz;
			f->data += t->datasz;
		}
		if (t->ptrs) {
			t->type = CAPN_PTR_LIST;
			t->len = t->ptrs;
			(*dep)++;
		}
		return 0;

	case CAPN_BIT_LIST:
		/* C11 7.24.1p2: a positive count needs both pointers. */
		if (t->datasz) {
			if (!t->data || !f->data)
				return -1;
			memcpy(t->data, f->data, t->datasz);
		}
		return 0;

	case CAPN_LIST:
		if (!t->len) {
			/* empty list - nothing to copy */
		} else if (t->ptrs && t->datasz) {
			(*dep)++;
		} else if (t->datasz) {
			/* INT31-C: a negative length is not a memcpy count. */
			if (t->len < 0 || t->len > INT_MAX / t->datasz)
				return -1;
			/* C11 7.24.1p2: a positive count needs both pointers. */
			if (!t->data || !f->data)
				return -1;
			memcpy(t->data, f->data, (size_t)t->len * (size_t)t->datasz);
		} else if (t->ptrs) {
			/* INT32-C: the pointer-list length is a signed product. */
			if (t->ptrs < 0 || t->len < 0 || t->len > INT_MAX / t->ptrs)
				return -1;
			t->type = CAPN_PTR_LIST;
			t->len *= t->ptrs;
			(*dep)++;
		}
		return 0;

	case CAPN_PTR_LIST:
		if (t->len) {
			(*dep)++;
		}
		return 0;

	default:
		return -1;
	}
}

static int copy_list_member(capn_ptr* t, capn_ptr *f, int *dep) {
	/* copy struct data */
	int sz = min(t->datasz, f->datasz);
	/* C11 7.24.1p2: a positive count needs both pointers. */
	if (sz) {
		if (!t->data || !f->data)
			return -1;
		memcpy(t->data, f->data, sz);
	}
	if (t->datasz > sz) {
		if (!t->data)
			return -1;
		memset(t->data + sz, 0, t->datasz - sz);
	}
	if (t->datasz)
		t->data += t->datasz;
	if (f->datasz)
		f->data += f->datasz;

	/* reset excess pointers */
	sz = min(t->ptrs, f->ptrs);
	if (t->ptrs > sz) {
		if (!t->data)
			return -1;
		memset(t->data + sz, 0, 8*(t->ptrs - sz));
	}

	/* create a pointer list for the main loop to copy */
	if (t->ptrs) {
		t->type = CAPN_PTR_LIST;
		t->len = t->ptrs;
		(*dep)++;
	}
	return 0;
}

#define MAX_COPY_DEPTH 32

/* TODO: handle CAPN_BIT_LIST and setting from an inner bit list member */
int capn_setp(capn_ptr p, int off, capn_ptr tgt) {
	struct capn_ptr to[MAX_COPY_DEPTH], from[MAX_COPY_DEPTH];
	char *data;
	int err, dep = 0;

	capn_resolve(&p);

	/* C11 6.5.3.2: a null segment or data pointer is not a load. */
	if (tgt.type == CAPN_FAR_POINTER && (!tgt.seg || !p.seg || !tgt.data || !p.data))
		return -1;
	if (tgt.type == CAPN_FAR_POINTER && tgt.seg->capn == p.seg->capn) {
		uint64_t val = capn_flip64(*(uint64_t*) tgt.data);
		if ((val & 3) == FAR_PTR) {
			*(uint64_t*) p.data = *(uint64_t*) tgt.data;
			return 0;
		}
	}

	capn_resolve(&tgt);

	switch (p.type) {
	case CAPN_LIST:
		/* INT32-C: the member offset is an index times the element size. */
		if (off < 0 || off >= p.len || tgt.type != CAPN_STRUCT)
			return -1;
		{
			uint64_t elem = (uint64_t)p.datasz + 8ull * (uint64_t)p.ptrs;
			uint64_t bytes;
			if (off > 0 && elem > UINT64_MAX / (uint64_t)off)
				return -1;
			bytes = (uint64_t)off * elem;
			to[0] = p;
			if (ptr_at(&p, bytes, &to[0].data))
				return -1;
		}
		from[0] = tgt;
		if (copy_list_member(to, from, &dep))
			return -1;
		break;

	case CAPN_PTR_LIST:
		if (off < 0 || off >= p.len)
			return -1;
		{
			uint64_t bytes = (uint64_t)off * 8ull;
			if (ptr_at(&p, bytes, &data))
				return -1;
		}
		goto copy_ptr;

	case CAPN_STRUCT:
		if (off < 0 || off >= p.ptrs)
			return -1;
		{
			uint64_t bytes = (uint64_t)p.datasz + 8ull * (uint64_t)off;
			if (ptr_at(&p, bytes, &data))
				return -1;
		}
		goto copy_ptr;

	copy_ptr:
		err = write_ptr(p.seg, data, tgt);
		if (err != NEED_TO_COPY)
			return err;

		/* Depth first copy the source whilst using a pointer stack to
		 * maintain the ptr to set and size left to copy at each level.
		 * We also maintain a rbtree (capn->copy) of the copies indexed
		 * by the source data. This way we can detect overlapped
		 * pointers in the source (and bail) and recursive structures
		 * (and point to the previous copy).
		 */

		from[0] = tgt;
		if (copy_ptr(p.seg, data, to, from, &dep))
			return -1;
		break;

	default:
		return -1;
	}

	while (dep) {
		struct capn_ptr *tc = &to[dep-1], *tn = &to[dep];
		struct capn_ptr *fc = &from[dep-1], *fn = &from[dep];

		if (dep+1 == MAX_COPY_DEPTH) {
			return -1;
		}

		if (tc->len < 0)
			return -1;
		if (!tc->len) {
			dep--;
			continue;
		}

		if (tc->type == CAPN_LIST) {
			*fn = capn_getp(*fc, 0, 1);
			*tn = capn_getp(*tc, 0, 1);

			if (copy_list_member(tn, fn, &dep))
				return -1;

			/* C11 6.5.6: do not step a null pointer, even by zero. */
			if (fc->datasz + 8*fc->ptrs) {
				if (!fc->data)
					return -1;
				fc->data += fc->datasz + 8*fc->ptrs;
			}
			if (tc->datasz + 8*tc->ptrs) {
				if (!tc->data)
					return -1;
				tc->data += tc->datasz + 8*tc->ptrs;
			}
			tc->len--;

		} else { /* CAPN_PTR_LIST */
			*fn = read_ptr(fc->seg, fc->data);

			if (fn->type && copy_ptr(tc->seg, tc->data, tn, fn, &dep))
				return -1;

			/* C11 6.5.6: a word step needs both cursors. */
			if (!fc->data || !tc->data)
				return -1;
			fc->data += 8;
			tc->data += 8;
			tc->len--;
		}
	}

	return 0;
}

/* TODO: handle CAPN_LIST, CAPN_PTR_LIST for bit lists */

int capn_get1(capn_list1 l, int off) {
	int bit;
	/* INT34-C: a negative index is not a shift count. */
	if (l.p.type != CAPN_BIT_LIST || off < 0 || off >= l.p.len || !l.p.data)
		return 0;
	/* C11 6.5.6: the byte index must lie in the bit buffer. */
	if ((uint64_t)(off / 8) >= (uint64_t)l.p.datasz)
		return 0;
	bit = off % 8;
	return (l.p.data[off / 8] & (1 << bit)) != 0;
}

int capn_set1(capn_list1 l, int off, int val) {
	int bit;
	/* INT34-C: a negative index is not a shift count. */
	if (l.p.type != CAPN_BIT_LIST || off < 0 || off >= l.p.len || !l.p.data)
		return -1;
	/* C11 6.5.6: the byte index must lie in the bit buffer. */
	if ((uint64_t)(off / 8) >= (uint64_t)l.p.datasz)
		return -1;
	bit = off % 8;
	if (val)
		l.p.data[off / 8] |= 1 << bit;
	else
		l.p.data[off / 8] &= ~(1 << bit);
	return 0;
}

/* C11 7.24.1p2: a zero count is skipped, and a count past the data fails. */
static int bit_slice(capn_ptr p, int *off, int sz, int *bsz, int *partial) {
	if (p.type != CAPN_BIT_LIST || sz < 0 || *off < 0 || (*off & 7) != 0)
		return -1;
	if (sz > INT_MAX - 7 || p.datasz < 0)
		return -1;
	*bsz = (sz + 7) / 8;
	*off /= 8;
	if (*off < 0 || *off > p.datasz)
		return -1;
	*partial = (uint64_t)*off + (uint64_t)sz > (uint64_t)p.datasz;
	return 0;
}

int capn_getv1(capn_list1 l, int off, uint8_t *data, int sz) {
	/* Note we only support aligned reads */
	int bsz, partial, n;
	capn_ptr p = l.p;
	if (bit_slice(p, &off, sz, &bsz, &partial))
		return -1;

	if (partial) {
		n = p.datasz - off;
		if (n > 0) {
			if (!data || !p.data)
				return -1;
			memcpy(data, p.data + off, (size_t)n);
		}
		return p.len - off * 8;
	}
	if (bsz > 0) {
		if (!data || !p.data)
			return -1;
		memcpy(data, p.data + off, (size_t)bsz);
	}
	return sz;
}

int capn_setv1(capn_list1 l, int off, const uint8_t *data, int sz) {
	/* Note we only support aligned writes */
	int bsz, partial, n;
	capn_ptr p = l.p;
	if (bit_slice(p, &off, sz, &bsz, &partial))
		return -1;

	if (partial) {
		n = p.datasz - off;
		if (n > 0) {
			if (!data || !p.data)
				return -1;
			memcpy(p.data + off, data, (size_t)n);
		}
		return p.len - off * 8;
	}
	if (bsz > 0) {
		if (!data || !p.data)
			return -1;
		memcpy(p.data + off, data, (size_t)bsz);
	}
	return sz;
}

/* pull out whether we add a tag or not as a define so the unit test can
 * test double far pointers by not creating tags */
#ifndef ADD_TAG
#define ADD_TAG 1
#endif

static void new_object(capn_ptr *p, int bytes) {
	struct capn_segment *s = p->seg;
	uint64_t aligned;

	if (!s || bytes < 0) {
		memset(p, 0, sizeof(*p));
		return;
	}

	if (!bytes)
		return;

	/* INT32-C: align and the segment length must fit in a signed int. */
	if (bytes > INT_MAX - 7) {
		memset(p, 0, sizeof(*p));
		return;
	}
	aligned = ((uint64_t)bytes + 7ull) & ~7ull;
	bytes = (int)aligned;

	if (s->len >= 0 && s->cap >= 0
	    && (uint64_t)s->len <= (uint64_t)INT_MAX - (uint64_t)bytes
	    && (uint64_t)s->len + (uint64_t)bytes <= (uint64_t)s->cap) {
		/* C11 6.5.6: do not add to a null data pointer. */
		if (!s->data) {
			memset(p, 0, sizeof(*p));
			return;
		}
		p->data = s->data + s->len;
		s->len += bytes;
		return;
	}

	/* add a tag whenever we switch segments so that write_ptr can
	 * use it */
	if (ADD_TAG && bytes > INT_MAX - 8) {
		memset(p, 0, sizeof(*p));
		return;
	}
	p->data = new_data(s->capn, bytes + ADD_TAG*8, &p->seg);
	if (!p->data) {
		memset(p, 0, sizeof(*p));
		return;
	}

	if (ADD_TAG) {
		write_ptr_tag(p->data, *p, 0);
		p->data += 8;
		p->has_ptr_tag = 1;
	}
}

capn_ptr capn_root(struct capn *c) {
	capn_ptr r = {CAPN_PTR_LIST};
	/* C11 6.5.3.2: a null message is not a segment list. */
	if (!c) {
		memset(&r, 0, sizeof(r));
		return r;
	}
	r.seg = lookup_segment(c, NULL, 0);
	r.data = r.seg ? r.seg->data : new_data(c, 8, &r.seg);
	r.len = 1;

	/* C11 6.5.6: a null data pointer is not a root base. */
	if (!r.seg || !r.data || r.seg->cap < 8) {
		memset(&r, 0, sizeof(r));
	} else if (r.seg->len < 8) {
		r.seg->len = 8;
	}

	return r;
}

capn_ptr capn_new_struct(struct capn_segment *seg, int datasz, int ptrs) {
	capn_ptr p = {CAPN_STRUCT};
	uint64_t aligned, nbytes;
	p.seg = seg;
	/* INT32-C: the aligned data size plus eight times the pointer count must fit. */
	if (datasz < 0 || ptrs < 0 || datasz > INT_MAX - 7) {
		memset(&p, 0, sizeof(p));
		return p;
	}
	aligned = ((uint64_t)datasz + 7ull) & ~7ull;
	nbytes = aligned + 8ull * (uint64_t)ptrs;
	if (nbytes > (uint64_t)INT_MAX) {
		memset(&p, 0, sizeof(p));
		return p;
	}
	p.datasz = (int)aligned;
	p.ptrs = ptrs;
	new_object(&p, (int)nbytes);
	return p;
}

capn_ptr capn_new_list(struct capn_segment *seg, int sz, int datasz, int ptrs) {
	capn_ptr p = {CAPN_LIST};
	p.seg = seg;
	p.len = sz;

	if (!sz) {
		/* empty lists may as well be a len=0 void list */
	} else if (ptrs || datasz > 8) {
		p.is_composite_list = 1;
		/* INT32-C: align and the byte count must fit in a signed int. */
		if (datasz > INT_MAX - 7 || sz < 0 || ptrs < 0) {
			memset(&p, 0, sizeof(p));
			return p;
		}
		p.datasz = (datasz + 7) & ~7;
		p.ptrs = ptrs;
		{
			uint64_t nbytes = (uint64_t)p.len * ((uint64_t)p.datasz + 8ull * (uint64_t)p.ptrs) + 8ull;
			if (nbytes > (uint64_t)INT_MAX) {
				memset(&p, 0, sizeof(p));
				return p;
			}
			new_object(&p, (int)nbytes);
		}
		if (p.data) {
			uint64_t hdr = STRUCT_PTR | (U64(p.len) << 2) | (U64(p.datasz/8) << 32) | (U64(ptrs) << 48);
			*(uint64_t*) p.data = capn_flip64(hdr);
			p.data += 8;
		}
	} else if (sz < 0 || datasz < 0) {
		memset(&p, 0, sizeof(p));
		return p;
	} else if (datasz > 4) {
		int bytes;
		p.datasz = 8;
		/* INT32-C: length times eight must fit in a signed int. */
		if ((uint64_t)p.len * 8ull > (uint64_t)INT_MAX) {
			memset(&p, 0, sizeof(p));
			return p;
		}
		bytes = p.len * 8;
		new_object(&p, bytes);
	} else if (datasz > 2) {
		int bytes;
		p.datasz = 4;
		if ((uint64_t)p.len * 4ull > (uint64_t)INT_MAX) {
			memset(&p, 0, sizeof(p));
			return p;
		}
		bytes = p.len * 4;
		new_object(&p, bytes);
	} else {
		int bytes;
		p.datasz = datasz;
		if (p.datasz > 0 && (uint64_t)p.len > (uint64_t)INT_MAX / (uint64_t)p.datasz) {
			memset(&p, 0, sizeof(p));
			return p;
		}
		bytes = p.len * p.datasz;
		new_object(&p, bytes);
	}

	return p;
}

capn_list1 capn_new_list1(struct capn_segment *seg, int sz) {
	capn_list1 l = {{CAPN_BIT_LIST}};
	l.p.seg = seg;
	l.p.len = sz;
	/* INT32-C: the bit-list byte count adds seven before dividing by eight. */
	if (sz < 0 || sz > INT_MAX - 7) {
		memset(&l.p, 0, sizeof(l.p));
		return l;
	}
	l.p.datasz = (sz + 7) / 8;
	new_object(&l.p, l.p.datasz);
	return l;
}

capn_ptr capn_new_ptr_list(struct capn_segment *seg, int sz) {
	capn_ptr p = {CAPN_PTR_LIST};
	p.seg = seg;
	p.len = sz;
	p.ptrs = 0;
	p.datasz = 0;
	/* INT32-C: length times eight must fit in a signed int. */
	if (sz < 0 || (uint64_t)sz * 8ull > (uint64_t)INT_MAX) {
		memset(&p, 0, sizeof(p));
		return p;
	}
	new_object(&p, sz * 8);
	return p;
}

capn_ptr capn_new_string(struct capn_segment *seg, const char *str, int sz) {
	capn_ptr p = {CAPN_LIST};
	uint64_t n;
	p.seg = seg;
	p.datasz = 1;
	/* C11 7.24.1p2: a null source is not a memcpy argument. */
	if (!str) {
		memset(&p, 0, sizeof(p));
		return p;
	}
	/* INT32-C: the length plus the trailing NUL must fit in a signed int. */
	if (sz >= 0) {
		n = (uint64_t)sz;
	} else {
		n = strlen(str);
	}
	if (n >= (uint64_t)INT_MAX) {
		memset(&p, 0, sizeof(p));
		return p;
	}
	p.len = (int)n + 1;
	new_object(&p, p.len);
	if (p.data) {
		memcpy(p.data, str, (size_t)p.len - 1);
	}
	return p;
}

capn_text capn_get_text(capn_ptr p, int off, capn_text def) {
	capn_ptr m = capn_getp(p, off, 1);
	capn_text ret = def;
	/* C11 6.5.6: a null list is not an index base. */
	if (m.type == CAPN_LIST && m.datasz == 1 && m.len > 0 && m.data
	    && m.data[m.len - 1] == 0) {
		ret.seg = m.seg;
		ret.str = m.data;
		ret.len = m.len - 1;
	}
	return ret;
}

int capn_set_text(capn_ptr p, int off, capn_text tgt) {
	capn_ptr m = {CAPN_NULL};
	if (tgt.seg) {
		uintptr_t base, at;
		uint64_t span;
		/* INT32-C: the length plus the trailing NUL must fit in a signed int. */
		if (tgt.len < 0 || tgt.len == INT_MAX)
			return -1;
		/* C11 6.5.3.2: a null string is not list data. */
		if (!tgt.str)
			return -1;
		/* C11 6.5.6: the string bytes must lie in the named segment. */
		if (!tgt.seg->data || tgt.seg->cap < 0)
			return -1;
		base = (uintptr_t) tgt.seg->data;
		at = (uintptr_t) tgt.str;
		if (at < base)
			return -1;
		span = (uint64_t) (at - base);
		if (span > (uint64_t) tgt.seg->cap
		    || (uint64_t) tgt.len + 1ull > (uint64_t) tgt.seg->cap - span)
			return -1;
		m.type = CAPN_LIST;
		m.seg = tgt.seg;
		m.data = (char*)tgt.str;
		m.len = tgt.len + 1;
		m.datasz = 1;
	} else if (tgt.str) {
		m = capn_new_string(p.seg, tgt.str, tgt.len);
	}
	return capn_setp(p, off, m);
}

capn_data capn_get_data(capn_ptr p, int off) {
	capn_data ret;
	ret.p = capn_getp(p, off, 1);
	if (ret.p.type != CAPN_LIST || ret.p.datasz != 1) {
		memset(&ret, 0, sizeof(ret));
	}
	return ret;
}

#define SZ 8
#include "capn-list.inc"
#undef SZ

#define SZ 16
#include "capn-list.inc"
#undef SZ

#define SZ 32
#include "capn-list.inc"
#undef SZ

#define SZ 64
#include "capn-list.inc"
#undef SZ
