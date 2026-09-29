
#define CAT2(A,B) A ## B
#define CAT(A,B) CAT2(A, B)
#define UINT_T CAT(CAT(uint, SZ), _t)
#define LIST_T CAT(capn_list, SZ)
#define FLIP CAT(capn_flip, SZ)

UINT_T CAT(capn_get,SZ) (LIST_T l, int off) {
	char *d;
	capn_ptr p = l.p;
	/* C11 6.5.6: a null list is not an element base. */
	if (off < 0 || off >= p.len || !p.data)
		return 0;

	switch (p.type) {
	case CAPN_LIST:
		if (p.datasz < SZ/8)
			return 0;
		{
			uint64_t stride = (uint64_t)p.datasz + 8ull * (uint64_t)p.ptrs;
			uint64_t byte_off = (uint64_t)off * stride;
			/* INT32-C: the index times the stride must fit in a signed int. */
			if (byte_off > (uint64_t)INT_MAX)
				return 0;
			d = p.data + (size_t)byte_off;
		}
		return FLIP(*(UINT_T*)d);

	case CAPN_PTR_LIST:
		d = struct_ptr(p.seg, p.data + 8*off, SZ/8);
		if (d) {
			return FLIP(*(UINT_T*)d);
		} else {
			return 0;
		}

	default:
		return 0;
	}
}

int CAT(capn_getv,SZ) (LIST_T l, int off, UINT_T *to, int sz) {
	int i;
	capn_ptr p = l.p;
	/* C11 6.5.6: a null list is not an element base. */
	if (off < 0 || !p.data || sz < 0 || p.len < 0)
		return -1;
	/* INT31-C: len minus off must not go negative. */
	if ((uint64_t)off > (uint64_t)p.len)
		sz = 0;
	else if ((uint64_t)sz > (uint64_t)p.len - (uint64_t)off)
		sz = (int)((uint64_t)p.len - (uint64_t)off);

	switch (p.type) {
	case CAPN_LIST:
		if (p.datasz == SZ/8 && !p.ptrs && (SZ == 8 || CAPN_LITTLE)) {
			uint64_t nbytes = (uint64_t)sz * (uint64_t)(SZ/8);
			/* INT32-C: the element count times the width must fit in a signed int. */
			if (nbytes > (uint64_t)INT_MAX)
				return -1;
			memcpy(to, p.data + off, (size_t)nbytes);
			return sz;
		} else if (p.datasz < SZ/8) {
			return -1;
		}

		for (i = 0; i < sz; i++) {
			uint64_t stride = (uint64_t)p.datasz + 8ull * (uint64_t)p.ptrs;
			uint64_t byte_off = ((uint64_t)i + (uint64_t)off) * stride;
			/* INT32-C: the index times the stride must fit in a signed int. */
			if (byte_off > (uint64_t)INT_MAX)
				return -1;
			char *d = p.data + (size_t)byte_off;
			to[i] = FLIP(*(UINT_T*)d);
		}
		return sz;

	case CAPN_PTR_LIST:
		for (i = 0; i < sz; i++) {
			char *d = struct_ptr(p.seg, p.data + 8*(i+off), SZ/8);
			if (d) {
				to[i] = FLIP(*(UINT_T*)d);
			} else {
				return -1;
			}
		}
		return sz;

	default:
		return -1;
	}
}

int CAT(capn_set,SZ) (LIST_T l, int off, UINT_T v) {
	char *d;
	capn_ptr p = l.p;
	/* C11 6.5.6: a null list is not an element base. */
	if (off < 0 || off >= p.len || !p.data)
		return -1;

	switch (p.type) {
	case CAPN_LIST:
		if (p.datasz < SZ/8)
			return -1;
		{
			uint64_t stride = (uint64_t)p.datasz + 8ull * (uint64_t)p.ptrs;
			uint64_t byte_off = (uint64_t)off * stride;
			/* INT32-C: the index times the stride must fit in a signed int. */
			if (byte_off > (uint64_t)INT_MAX)
				return -1;
			d = p.data + (size_t)byte_off;
		}
		*(UINT_T*) d = FLIP(v);
		return 0;

	case CAPN_PTR_LIST:
		d = struct_ptr(p.seg, p.data + 8*off, SZ/8);
		if (!d) {
			return -1;
		}
		*(UINT_T*) d = FLIP(v);
		return 0;

	default:
		return -1;
	}
}

int CAT(capn_setv,SZ) (LIST_T l, int off, const UINT_T *from, int sz) {
	int i;
	capn_ptr p = l.p;
	/* C11 6.5.6: a null list is not an element base. */
	if (off < 0 || !p.data || sz < 0 || p.len < 0)
		return -1;
	/* INT31-C: len minus off must not go negative. */
	if ((uint64_t)off > (uint64_t)p.len)
		sz = 0;
	else if ((uint64_t)sz > (uint64_t)p.len - (uint64_t)off)
		sz = (int)((uint64_t)p.len - (uint64_t)off);

	switch (p.type) {
	case CAPN_LIST:
		if (p.datasz == SZ/8 && !p.ptrs && (SZ == 8 || CAPN_LITTLE)) {
			uint64_t nbytes = (uint64_t)sz * (uint64_t)(SZ/8);
			/* INT32-C: the element count times the width must fit in a signed int. */
			if (nbytes > (uint64_t)INT_MAX)
				return -1;
			memcpy(p.data + off, from, (size_t)nbytes);
			return sz;
		} else if (p.datasz < SZ/8) {
			return -1;
		}

		for (i = 0; i < sz; i++) {
			uint64_t stride = (uint64_t)p.datasz + 8ull * (uint64_t)p.ptrs;
			uint64_t byte_off = ((uint64_t)i + (uint64_t)off) * stride;
			/* INT32-C: the index times the stride must fit in a signed int. */
			if (byte_off > (uint64_t)INT_MAX)
				return -1;
			char *d = p.data + (size_t)byte_off;
			*(UINT_T*) d = FLIP(from[i]);
		}
		return sz;

	case CAPN_PTR_LIST:
		for (i = 0; i < sz; i++) {
			char *d = struct_ptr(p.seg, p.data + 8*(i+off), SZ/8);
			if (d) {
				*(UINT_T*) d = FLIP(from[i]);
			} else {
				return -1;
			}
		}
		return sz;

	default:
		return -1;
	}
}

LIST_T CAT(capn_new_list,SZ) (struct capn_segment *seg, int sz) {
	LIST_T l = {{CAPN_LIST}};
	l.p.seg = seg;
	l.p.len = sz;
	l.p.datasz = SZ/8;
	/* INT32-C: the element count times the width must fit in a signed int. */
	if (sz < 0 || (uint64_t)sz > (uint64_t)INT_MAX / (uint64_t)(SZ/8)) {
		memset(&l.p, 0, sizeof(l.p));
		return l;
	}
	new_object(&l.p, sz*(SZ/8));
	return l;
}

#undef CAT2
#undef CAT
#undef UINT_T
#undef LIST_T
#undef FLIP

