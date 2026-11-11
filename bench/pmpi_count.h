/* SPDX-License-Identifier: MIT */
/* Copyright (c) cpmdc developers */
#ifndef CPMDC_PMPI_COUNT_H
#define CPMDC_PMPI_COUNT_H

/* PMPI shim that counts collective calls on this rank.
 * Link it ahead of the MPI library, or preload the shared build.
 * PMPI_COUNT_REPORT=1 prints this rank's totals from MPI_Finalize.
 */

enum {
  PMPI_COUNT_ALLREDUCE,
  PMPI_COUNT_REDUCE,
  PMPI_COUNT_BCAST,
  PMPI_COUNT_BARRIER,
  PMPI_COUNT_GATHER,
  PMPI_COUNT_GATHERV,
  PMPI_COUNT_ALLGATHER,
  PMPI_COUNT_ALLGATHERV,
  PMPI_COUNT_SCATTER,
  PMPI_COUNT_ALLTOALL,
  PMPI_COUNT_REDUCE_SCATTER_BLOCK,
  PMPI_COUNT_SCAN,
  PMPI_COUNT_EXSCAN,
  PMPI_COUNT_IALLREDUCE,
  PMPI_COUNT_IBCAST,
  PMPI_COUNT_IBARRIER,
  PMPI_COUNT_COMM_DUP,
  PMPI_COUNT_COMM_SPLIT,
  PMPI_COUNT_KINDS
};

const char *pmpi_count_name(int kind);

/* Copies this rank's running totals into out[PMPI_COUNT_KINDS]. */
void pmpi_count_snapshot(long long *out);

#endif
