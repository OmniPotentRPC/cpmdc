/* SPDX-License-Identifier: MIT */
/* Copyright (c) cpmdc developers */
#include "pmpi_count.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static long long g_counts[PMPI_COUNT_KINDS];

static const char *const g_names[PMPI_COUNT_KINDS] = {
    "Allreduce", "Reduce",     "Bcast",     "Barrier",
    "Gather",    "Gatherv",    "Allgather", "Allgatherv",
    "Scatter",   "Alltoall",   "Reduce_scatter_block",
    "Scan",      "Exscan",     "Iallreduce", "Ibcast",
    "Ibarrier",  "Comm_dup",   "Comm_split"};

const char *pmpi_count_name(int kind) {
  if (kind < 0 || kind >= PMPI_COUNT_KINDS)
    return "unknown";
  return g_names[kind];
}

void pmpi_count_snapshot(long long *out) {
  if (out)
    memcpy(out, g_counts, sizeof(g_counts));
}

int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  g_counts[PMPI_COUNT_ALLREDUCE] += 1;
  return PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

int MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
  g_counts[PMPI_COUNT_REDUCE] += 1;
  return PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
}

int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
              MPI_Comm comm) {
  g_counts[PMPI_COUNT_BCAST] += 1;
  return PMPI_Bcast(buffer, count, datatype, root, comm);
}

int MPI_Barrier(MPI_Comm comm) {
  g_counts[PMPI_COUNT_BARRIER] += 1;
  return PMPI_Barrier(comm);
}

int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
               MPI_Comm comm) {
  g_counts[PMPI_COUNT_GATHER] += 1;
  return PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                     recvtype, root, comm);
}

int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, const int recvcounts[], const int displs[],
                MPI_Datatype recvtype, int root, MPI_Comm comm) {
  g_counts[PMPI_COUNT_GATHERV] += 1;
  return PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                      displs, recvtype, root, comm);
}

int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm) {
  g_counts[PMPI_COUNT_ALLGATHER] += 1;
  return PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                        recvtype, comm);
}

int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, const int recvcounts[], const int displs[],
                   MPI_Datatype recvtype, MPI_Comm comm) {
  g_counts[PMPI_COUNT_ALLGATHERV] += 1;
  return PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                         displs, recvtype, comm);
}

int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
                MPI_Comm comm) {
  g_counts[PMPI_COUNT_SCATTER] += 1;
  return PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                      recvtype, root, comm);
}

int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 MPI_Comm comm) {
  g_counts[PMPI_COUNT_ALLTOALL] += 1;
  return PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                       recvtype, comm);
}

int MPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf,
                             int recvcount, MPI_Datatype datatype, MPI_Op op,
                             MPI_Comm comm) {
  g_counts[PMPI_COUNT_REDUCE_SCATTER_BLOCK] += 1;
  return PMPI_Reduce_scatter_block(sendbuf, recvbuf, recvcount, datatype, op,
                                   comm);
}

int MPI_Scan(const void *sendbuf, void *recvbuf, int count,
             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  g_counts[PMPI_COUNT_SCAN] += 1;
  return PMPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
}

int MPI_Exscan(const void *sendbuf, void *recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  g_counts[PMPI_COUNT_EXSCAN] += 1;
  return PMPI_Exscan(sendbuf, recvbuf, count, datatype, op, comm);
}

int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                   MPI_Request *request) {
  g_counts[PMPI_COUNT_IALLREDUCE] += 1;
  return PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
}

int MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root,
               MPI_Comm comm, MPI_Request *request) {
  g_counts[PMPI_COUNT_IBCAST] += 1;
  return PMPI_Ibcast(buffer, count, datatype, root, comm, request);
}

int MPI_Ibarrier(MPI_Comm comm, MPI_Request *request) {
  g_counts[PMPI_COUNT_IBARRIER] += 1;
  return PMPI_Ibarrier(comm, request);
}

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm) {
  g_counts[PMPI_COUNT_COMM_DUP] += 1;
  return PMPI_Comm_dup(comm, newcomm);
}

int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) {
  g_counts[PMPI_COUNT_COMM_SPLIT] += 1;
  return PMPI_Comm_split(comm, color, key, newcomm);
}

int MPI_Finalize(void) {
  const char *env = getenv("PMPI_COUNT_REPORT");
  if (env && strcmp(env, "1") == 0) {
    int rank = 0, k;
    long long total = 0;
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (k = 0; k < PMPI_COUNT_KINDS; ++k) {
      if (g_counts[k] == 0)
        continue;
      total += g_counts[k];
      fprintf(stderr, "pmpi_count rank %d %s %lld\n", rank, g_names[k],
              g_counts[k]);
    }
    fprintf(stderr, "pmpi_count rank %d total %lld\n", rank, total);
  }
  return PMPI_Finalize();
}
