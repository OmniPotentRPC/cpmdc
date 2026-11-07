/* SPDX-License-Identifier: MIT */
#define _POSIX_C_SOURCE 200809L

#include "cpmdc_pwscale.h"
#include "pwscale_fft.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/utsname.h>
#include <unistd.h>

#if defined(CPMDC_PWSCALE_HAVE_MPI)
#include <mpi.h>
#endif

static const char *HEADER =
    "kind,mode,dist,scaling,backend,parallel,workers,run,natoms,npw,cell_bohr,alpha,"
    "sigma,seconds,plan_s,energy_h,max_abs_dE,max_abs_dF,g_visited,atom_visits,"
    "alloc_bytes,h2d_s,d2h_s,counts_match,omp_threads,mpi_ranks,team_observed,"
    "cpu_model,nproc,uname_m,gpu\n";

typedef struct hw_info {
  char model[192];
  char machine[128];
  char gpu[192];
  long nproc;
} hw_info;

typedef struct row {
  const char *kind;
  const char *mode;
  const char *dist;
  const char *scaling;
  const char *backend;
  const char *parallel;
  int workers;
  int run;
  int natoms;
  int npw;
  double cell;
  double alpha;
  double sigma;
  int has_seconds;
  double seconds;
  int has_plan;
  double plan_s;
  int has_energy;
  double energy_h;
  int has_de;
  double max_abs_dE;
  int has_df;
  double max_abs_dF;
  long long g_visited;
  long long atom_visits;
  long long alloc_bytes;
  int has_h2d;
  double h2d_s;
  int has_d2h;
  double d2h_s;
  int counts_match;
  int omp_threads;
  int mpi_ranks;
  int team_observed;
  const hw_info *hw;
} row;

static void csv_text(FILE *f, const char *s) {
  const char *p;
  fputc('"', f);
  for (p = s ? s : ""; *p; ++p) {
    if (*p == '"')
      fputc('"', f);
    fputc(*p, f);
  }
  fputc('"', f);
}

static void csv_num_or_blank(FILE *f, int present, double value) {
  if (!present)
    return;
  fprintf(f, "%.9g", value);
}

static void write_row(FILE *f, const row *r) {
  csv_text(f, r->kind);
  fputc(',', f);
  csv_text(f, r->mode);
  fputc(',', f);
  csv_text(f, r->dist);
  fputc(',', f);
  csv_text(f, r->scaling);
  fputc(',', f);
  csv_text(f, r->backend);
  fputc(',', f);
  csv_text(f, r->parallel);
  fputc(',', f);
  fprintf(f, "%d,%d,%d,%d,%.9g,%.9g,%.9g,", r->workers, r->run, r->natoms, r->npw,
          r->cell, r->alpha, r->sigma);
  csv_num_or_blank(f, r->has_seconds, r->seconds);
  fputc(',', f);
  csv_num_or_blank(f, r->has_plan, r->plan_s);
  fputc(',', f);
  csv_num_or_blank(f, r->has_energy, r->energy_h);
  fputc(',', f);
  csv_num_or_blank(f, r->has_de, r->max_abs_dE);
  fputc(',', f);
  csv_num_or_blank(f, r->has_df, r->max_abs_dF);
  fputc(',', f);
  fprintf(f, "%lld,%lld,%lld,", r->g_visited, r->atom_visits, r->alloc_bytes);
  csv_num_or_blank(f, r->has_h2d, r->h2d_s);
  fputc(',', f);
  csv_num_or_blank(f, r->has_d2h, r->d2h_s);
  fputc(',', f);
  fprintf(f, "%d,%d,%d,%d,", r->counts_match, r->omp_threads, r->mpi_ranks,
          r->team_observed);
  csv_text(f, r->hw->model);
  fputc(',', f);
  fprintf(f, "%ld,", r->hw->nproc);
  csv_text(f, r->hw->machine);
  fputc(',', f);
  csv_text(f, r->hw->gpu);
  fputc('\n', f);
  fflush(f);
}

static void load_hw(hw_info *hw) {
  FILE *cpu;
  char line[256];
  memset(hw, 0, sizeof(*hw));
  snprintf(hw->model, sizeof(hw->model), "unknown");
  snprintf(hw->machine, sizeof(hw->machine), "unknown");
  snprintf(hw->gpu, sizeof(hw->gpu), "none");
  hw->nproc = sysconf(_SC_NPROCESSORS_ONLN);
  {
    struct utsname uts;
    if (uname(&uts) == 0)
      snprintf(hw->machine, sizeof(hw->machine), "%s", uts.machine);
  }
  cpu = fopen("/proc/cpuinfo", "r");
  if (cpu) {
    while (fgets(line, sizeof(line), cpu)) {
      if (strncmp(line, "model name", 10) == 0) {
        char *p = strchr(line, ':');
        if (p) {
          p += 1;
          while (*p == ' ' || *p == '\t')
            p += 1;
          snprintf(hw->model, sizeof(hw->model), "%s", p);
          hw->model[strcspn(hw->model, "\r\n")] = '\0';
        }
        break;
      }
    }
    fclose(cpu);
  }
  if (cpmdc_pwscale_cuda_device(hw->gpu, (int)sizeof(hw->gpu)) != 1)
    snprintf(hw->gpu, sizeof(hw->gpu), "none");
}

static int parse_threads(const char *text, int *out, int cap) {
  int n = 0;
  const char *p = text;
  while (*p && n < cap) {
    char *end = NULL;
    long v = strtol(p, &end, 10);
    if (end == p || v < 1 || v > 4096)
      return -1;
    out[n++] = (int)v;
    if (*end == ',')
      p = end + 1;
    else if (*end == '\0')
      break;
    else
      return -1;
  }
  return n;
}

static double peak_gap(const double *got, const double *ref, int n) {
  int i;
  double peak = 0.0;
  for (i = 0; i < n; ++i) {
    double gap = fabs(got[i] - ref[i]);
    if (gap > peak)
      peak = gap;
  }
  return peak;
}

static int forces_ok(const double *got, const double *ref, int n) {
  int i;
  for (i = 0; i < n; ++i) {
    if (!cpmdc_pwscale_within(got[i], ref[i]))
      return 0;
  }
  return 1;
}

static void maybe_poison(cpmdc_pwscale_stats *got, int run) {
  const char *env = getenv("CPMDC_PWSCALE_POISON_RUN");
  if (env && atoi(env) == run)
    got->energy_h += 1.0;
}

/* Fill the row from this sample. Returns 1 when the sample agrees. */
static int record_force_sample(row *line, const cpmdc_pwscale_stats *got,
                               const cpmdc_pwscale_stats *serial,
                               const double *forces, const double *ref, int n3,
                               int requested, int use_mpi, double *de, double *df) {
  int counts = 0;
  int ok = cpmdc_pwscale_sample_check(got->energy_h, forces, serial->energy_h, ref,
                                      n3, got->g_visited, serial->g_visited,
                                      got->atom_visits, serial->atom_visits, de, df,
                                      &counts);
  line->has_energy = 1;
  line->energy_h = got->energy_h;
  line->has_de = 1;
  line->max_abs_dE = *de;
  line->has_df = 1;
  line->max_abs_dF = *df;
  line->g_visited = got->g_visited;
  line->atom_visits = got->atom_visits;
  line->alloc_bytes = got->alloc_bytes;
  line->counts_match = counts;
  line->omp_threads = got->workers;
  line->mpi_ranks = got->mpi_ranks;
  if (use_mpi) {
    line->workers = got->mpi_ranks;
    line->team_observed = got->mpi_ranks == requested && got->workers == 1;
  } else {
    line->workers = got->workers;
    line->team_observed = got->workers == requested;
  }
  if (!line->team_observed)
    return 0;
  return ok;
}

static int run_force(const char *mode, const char *dist_name, const char *scaling,
                     const char *axis, int dist, int natoms, int npw,
                     double cell, double alpha, double sigma, int workers, int runs,
                     int use_mpi, FILE *out, const hw_info *hw, int rank) {
  cpmdc_pwscale_problem prob;
  double *pos = NULL, *f_serial = NULL, *f_par = NULL;
  int *z = NULL;
  cpmdc_pwscale_stats serial, got;
  row line;
  int rc, run, n3, bad = 0;
  double de = 0.0, df = 0.0;

  memset(&line, 0, sizeof(line));
  prob.natoms = natoms;
  prob.npw = npw;
  prob.cell_bohr = cell;
  prob.alpha = alpha;
  prob.sigma = sigma;
  n3 = natoms * 3;
  pos = (double *)calloc((size_t)n3, sizeof(double));
  f_serial = (double *)calloc((size_t)n3, sizeof(double));
  f_par = (double *)calloc((size_t)n3, sizeof(double));
  z = (int *)calloc((size_t)natoms, sizeof(int));
  if (!pos || !f_serial || !f_par || !z)
    bad = 1;
#if defined(CPMDC_PWSCALE_HAVE_MPI)
  if (use_mpi) {
    int ready = bad ? 0 : 1;
    int all_ready = 0;
    MPI_Allreduce(&ready, &all_ready, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!all_ready)
      bad = 1;
  }
#endif
  if (bad)
    goto done;
  cpmdc_pwscale_fill_geometry(&prob, pos, z);
  rc = cpmdc_pwscale_eval(&prob, pos, z, f_serial, 1, dist, &serial);
  if (rc != CPMDC_PWSCALE_OK) {
    bad = 1;
    goto done;
  }
  if (use_mpi)
    rc = cpmdc_pwscale_eval_mpi(&prob, pos, z, f_par, 1, dist, &got, NULL, 0);
  else
    rc = cpmdc_pwscale_eval(&prob, pos, z, f_par, workers, dist, &got);
  if (rc != CPMDC_PWSCALE_OK) {
    bad = 1;
    goto done;
  }
  memset(&line, 0, sizeof(line));
  line.mode = mode;
  line.dist = dist_name;
  line.scaling = scaling;
  line.backend = "cpu";
  line.parallel = axis;
  line.natoms = natoms;
  line.npw = npw;
  line.cell = cell;
  line.alpha = alpha;
  line.sigma = sigma;
  line.hw = hw;
  if (!record_force_sample(&line, &got, &serial, f_par, f_serial, n3, workers,
                           use_mpi, &de, &df))
    bad = 1;
#if defined(CPMDC_PWSCALE_HAVE_MPI)
  if (use_mpi) {
    int any = 0;
    MPI_Allreduce(&bad, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    bad = any;
  }
#endif
  if (rank == 0 && out) {
    line.kind = "check";
    line.run = 0;
    write_row(out, &line);
  }
  if (bad) {
    if (rank == 0)
      fprintf(stderr, "agreement failed for %s workers %d dE %.3e dF %.3e\n", mode,
              line.workers, de, df);
    goto done;
  }
  for (run = 1; run <= runs; ++run) {
    int sample_ok;
    if (use_mpi)
      rc = cpmdc_pwscale_eval_mpi(&prob, pos, z, f_par, 1, dist, &got, NULL, 0);
    else
      rc = cpmdc_pwscale_eval(&prob, pos, z, f_par, workers, dist, &got);
    if (rc != CPMDC_PWSCALE_OK) {
      bad = 1;
      break;
    }
    maybe_poison(&got, run);
    sample_ok = record_force_sample(&line, &got, &serial, f_par, f_serial, n3,
                                    workers, use_mpi, &de, &df);
    if (!sample_ok)
      bad = 1;
    if (rank == 0 && out) {
      line.kind = "time";
      line.run = run;
      line.has_seconds = 1;
      line.seconds = got.seconds;
      write_row(out, &line);
    }
    if (!sample_ok)
      break;
  }
#if defined(CPMDC_PWSCALE_HAVE_MPI)
  if (use_mpi) {
    int any = 0;
    MPI_Allreduce(&bad, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    bad = any;
  }
#endif
done:
  free(pos);
  free(f_serial);
  free(f_par);
  free(z);
  return bad ? 1 : 0;
}

static int record_fft_sample(row *line, int weak, int workers, int points,
                             int actual_threads, long long points_executed,
                             double roundtrip, double versus) {
  const long long batch = weak ? (workers < 1 ? 1 : workers) : 1;
  const long long expected = batch * (long long)points;
  line->g_visited = points_executed;
  line->counts_match = points_executed == expected;
  line->has_de = 1;
  line->max_abs_dE = roundtrip;
  line->has_df = 1;
  line->max_abs_dF = versus;
  line->mpi_ranks = 1;
  if (weak) {
    line->omp_threads = actual_threads;
    line->workers = actual_threads;
    line->team_observed = actual_threads == workers;
  } else {
    /* FFTW's internal team is not visible. workers stays the plan request. */
    line->omp_threads = 0;
    line->workers = workers;
    line->team_observed = 0;
  }
  if (weak && !line->team_observed)
    return 0;
  if (!line->counts_match || roundtrip > 1e-8 || versus > 1e-8)
    return 0;
  return 1;
}

static int run_fft_mode(const char *mode, const char *scaling, int weak, int n,
                        int workers, int runs, FILE *out, const hw_info *hw) {
  row line;
  int run, rc, bad = 0, actual_threads = 0;
  double seconds = 0, plan_s = 0, roundtrip = 0, versus = 0;
  long long alloc_bytes = 0, points_executed = 0;
  const int points = n * n * n;

  rc = cpmdc_pwscale_fft(n, workers, weak, &seconds, &plan_s, &roundtrip, &versus,
                         &alloc_bytes, &actual_threads, &points_executed);
  if (rc != CPMDC_PWSCALE_OK)
    return rc;
  memset(&line, 0, sizeof(line));
  line.mode = mode;
  line.dist = "fft";
  line.scaling = scaling;
  line.backend = "cpu";
  line.parallel = "threads";
  line.npw = points;
  line.alloc_bytes = alloc_bytes;
  line.hw = hw;
  if (!record_fft_sample(&line, weak, workers, points, actual_threads,
                         points_executed, roundtrip, versus))
    bad = 1;
  line.kind = "check";
  line.run = 0;
  line.has_plan = 1;
  line.plan_s = plan_s;
  if (out)
    write_row(out, &line);
  if (bad) {
    fprintf(stderr,
            "fft agreement failed n %d workers %d roundtrip %.3e versus %.3e "
            "points %lld team %d\n",
            n, workers, roundtrip, versus, points_executed, actual_threads);
    return 1;
  }
  for (run = 1; run <= runs; ++run) {
    rc = cpmdc_pwscale_fft(n, workers, weak, &seconds, &plan_s, &roundtrip, &versus,
                           &alloc_bytes, &actual_threads, &points_executed);
    if (rc != CPMDC_PWSCALE_OK)
      return rc;
    if (!record_fft_sample(&line, weak, workers, points, actual_threads,
                           points_executed, roundtrip, versus))
      bad = 1;
    line.kind = "time";
    line.run = run;
    line.has_seconds = 1;
    line.seconds = seconds;
    line.has_plan = 1;
    line.plan_s = plan_s;
    line.alloc_bytes = alloc_bytes;
    if (out)
      write_row(out, &line);
    if (bad) {
      fprintf(stderr,
              "fft timed sample failed n %d run %d roundtrip %.3e versus %.3e "
              "points %lld\n",
              n, run, roundtrip, versus, points_executed);
      return 1;
    }
  }
  return 0;
}

static int write_todo(FILE *out, const hw_info *hw, const char *reason) {
  row line;
  if (!out)
    return 0;
  memset(&line, 0, sizeof(line));
  line.kind = "todo";
  line.mode = "cuda";
  line.dist = "pw";
  line.scaling = "strong";
  line.backend = "cuda";
  line.parallel = "device";
  line.workers = 0;
  line.hw = hw;
  write_row(out, &line);
  fprintf(stderr, "%s\n", reason);
  return 2;
}

static void usage(void) {
  fprintf(stderr,
          "cpmdc-pwscale --mode <force-strong|force-weak-atoms|force-weak-pw|"
          "fft-strong|fft-weak|cuda> --out FILE [--threads 1,2,4,8] "
          "[--ranks] [--runs 5] [--natoms N] [--npw N] [--atoms-per N] "
          "[--pw-per N] [--fft-n N] [--cell B] [--alpha A] [--sigma S]\n");
}

int main(int argc, char **argv) {
  const char *mode = NULL;
  const char *out_path = NULL;
  int threads[16];
  int nthreads = 0;
  int ranks = 0;
  int runs = 5;
  int natoms = 64;
  int npw = 8192;
  int atoms_per = 64;
  int pw_per = 2048;
  int fft_n = 64;
  double cell = 12.0, alpha = 0.4, sigma = 0.9;
  int i, rc = 2, rank = 0, size = 1;
  hw_info hw;
  FILE *out = NULL;

  for (i = 1; i < argc; ++i) {
    const char *a = argv[i];
    const char *v = (i + 1 < argc) ? argv[i + 1] : NULL;
    if (strcmp(a, "--mode") == 0 && v) {
      mode = v;
      i += 1;
    } else if (strcmp(a, "--out") == 0 && v) {
      out_path = v;
      i += 1;
    } else if (strcmp(a, "--threads") == 0 && v) {
      nthreads = parse_threads(v, threads, 16);
      if (nthreads < 1)
        return 2;
      i += 1;
    } else if (strcmp(a, "--ranks") == 0) {
      ranks = 1;
    } else if (strcmp(a, "--runs") == 0 && v) {
      runs = atoi(v);
      i += 1;
    } else if (strcmp(a, "--natoms") == 0 && v) {
      natoms = atoi(v);
      i += 1;
    } else if (strcmp(a, "--npw") == 0 && v) {
      npw = atoi(v);
      i += 1;
    } else if (strcmp(a, "--atoms-per") == 0 && v) {
      atoms_per = atoi(v);
      i += 1;
    } else if (strcmp(a, "--pw-per") == 0 && v) {
      pw_per = atoi(v);
      i += 1;
    } else if (strcmp(a, "--fft-n") == 0 && v) {
      fft_n = atoi(v);
      i += 1;
    } else if (strcmp(a, "--cell") == 0 && v) {
      cell = atof(v);
      i += 1;
    } else if (strcmp(a, "--alpha") == 0 && v) {
      alpha = atof(v);
      i += 1;
    } else if (strcmp(a, "--sigma") == 0 && v) {
      sigma = atof(v);
      i += 1;
    } else if (strcmp(a, "--help") == 0) {
      usage();
      return 0;
    } else {
      usage();
      return 2;
    }
  }
  if (!mode || !out_path || runs < 1) {
    usage();
    return 2;
  }
  if (nthreads == 0) {
    threads[0] = 1;
    nthreads = 1;
  }

#if defined(CPMDC_PWSCALE_HAVE_MPI)
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  load_hw(&hw);
  {
    int open_ok = 1;
    if (rank == 0) {
      out = fopen(out_path, "a+");
      if (!out) {
        fprintf(stderr, "cannot open output\n");
        open_ok = 0;
      } else if (fseek(out, 0, SEEK_END) != 0 || ftell(out) == 0) {
        fputs(HEADER, out);
      }
    }
#if defined(CPMDC_PWSCALE_HAVE_MPI)
    MPI_Bcast(&open_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (!open_ok) {
      rc = 2;
      goto finish;
    }
  }

  if (strcmp(mode, "cuda") == 0) {
    if (rank == 0) {
      if (!cpmdc_pwscale_cuda_compiled())
        rc = write_todo(out, &hw, "cuda backend is not part of this build");
      else if (!cpmdc_pwscale_cuda_device(NULL, 0))
        rc = write_todo(out, &hw, "cuda backend found no device");
      else {
        cpmdc_pwscale_problem prob;
        double *pos, *f_cpu, *f_gpu;
        int *z;
        cpmdc_pwscale_stats cpu, gpu;
        int n3 = natoms * 3;
        prob.natoms = natoms;
        prob.npw = npw;
        prob.cell_bohr = cell;
        prob.alpha = alpha;
        prob.sigma = sigma;
        pos = (double *)calloc((size_t)n3, sizeof(double));
        f_cpu = (double *)calloc((size_t)n3, sizeof(double));
        f_gpu = (double *)calloc((size_t)n3, sizeof(double));
        z = (int *)calloc((size_t)natoms, sizeof(int));
        cpmdc_pwscale_fill_geometry(&prob, pos, z);
        cpmdc_pwscale_eval(&prob, pos, z, f_cpu, 1, CPMDC_PWSCALE_DIST_PW, &cpu);
        rc = 0;
        for (i = 0; i < runs; ++i) {
          row line;
          int grc = cpmdc_pwscale_eval_cuda(&prob, pos, z, f_gpu, &gpu);
          if (grc != CPMDC_PWSCALE_OK) {
            rc = write_todo(out, &hw, "cuda backend failed during evaluation");
            break;
          }
          if (!cpmdc_pwscale_within(gpu.energy_h, cpu.energy_h) ||
              !forces_ok(f_gpu, f_cpu, n3)) {
            fprintf(stderr, "cuda result is outside tolerance\n");
            rc = 1;
            break;
          }
          memset(&line, 0, sizeof(line));
          line.kind = "time";
          line.mode = "cuda";
          line.dist = "pw";
          line.scaling = "strong";
          line.backend = "cuda";
          line.parallel = "device";
          line.workers = 1;
          line.run = i + 1;
          line.natoms = natoms;
          line.npw = npw;
          line.cell = cell;
          line.alpha = alpha;
          line.sigma = sigma;
          line.has_seconds = 1;
          line.seconds = gpu.seconds;
          line.has_energy = 1;
          line.energy_h = gpu.energy_h;
          line.has_de = 1;
          line.max_abs_dE = fabs(gpu.energy_h - cpu.energy_h);
          line.has_df = 1;
          line.max_abs_dF = peak_gap(f_gpu, f_cpu, n3);
          line.g_visited = gpu.g_visited;
          line.atom_visits = gpu.atom_visits;
          line.alloc_bytes = gpu.alloc_bytes;
          line.has_h2d = 1;
          line.h2d_s = gpu.h2d_seconds;
          line.has_d2h = 1;
          line.d2h_s = gpu.d2h_seconds;
          line.counts_match = gpu.g_visited == cpu.g_visited &&
                              gpu.atom_visits == cpu.atom_visits;
          line.hw = &hw;
          write_row(out, &line);
        }
        free(pos);
        free(f_cpu);
        free(f_gpu);
        free(z);
      }
    } else {
      rc = 0;
    }
    goto finish;
  }

  if (strcmp(mode, "force-strong") == 0 || strcmp(mode, "force-weak-atoms") == 0 ||
      strcmp(mode, "force-weak-pw") == 0) {
    int dist = CPMDC_PWSCALE_DIST_PW;
    const char *dist_name = "pw";
    const char *scaling = "strong";
    if (strcmp(mode, "force-weak-atoms") == 0) {
      dist = CPMDC_PWSCALE_DIST_ATOMS;
      dist_name = "atoms";
      scaling = "weak";
    } else if (strcmp(mode, "force-weak-pw") == 0) {
      scaling = "weak";
    }
    if (ranks) {
      int use_natoms = natoms;
      int use_npw = npw;
      if (strcmp(mode, "force-weak-atoms") == 0)
        use_natoms = atoms_per * size;
      else if (strcmp(mode, "force-weak-pw") == 0)
        use_npw = pw_per * size;
      rc = run_force(mode, dist_name, scaling, "ranks", dist, use_natoms, use_npw,
                     cell, alpha, sigma, size, runs, 1, out, &hw, rank);
    } else {
      if (size != 1) {
        if (rank == 0)
          fprintf(stderr, "thread sweep expects one rank\n");
        rc = 2;
        goto finish;
      }
      rc = 0;
      for (i = 0; i < nthreads; ++i) {
        int use_natoms = natoms;
        int use_npw = npw;
        int workers = threads[i];
        if (strcmp(mode, "force-weak-atoms") == 0)
          use_natoms = atoms_per * workers;
        else if (strcmp(mode, "force-weak-pw") == 0)
          use_npw = pw_per * workers;
        rc = run_force(mode, dist_name, scaling, "threads", dist, use_natoms, use_npw,
                       cell, alpha, sigma, workers, runs, 0, out, &hw, rank);
        if (rc != 0)
          break;
      }
    }
    goto finish;
  }

  if (strcmp(mode, "fft-strong") == 0 || strcmp(mode, "fft-weak") == 0) {
    int weak = strcmp(mode, "fft-weak") == 0;
    const char *scaling = weak ? "weak" : "strong";
    if (ranks || size != 1) {
      if (rank == 0)
        fprintf(stderr, "fft sweep is threaded\n");
      rc = 2;
      goto finish;
    }
    rc = 0;
    for (i = 0; i < nthreads; ++i) {
      rc = run_fft_mode(mode, scaling, weak, fft_n, threads[i], runs, out, &hw);
      if (rc != 0)
        break;
    }
    goto finish;
  }

  usage();
  rc = 2;

finish:
  if (out)
    fclose(out);
#if defined(CPMDC_PWSCALE_HAVE_MPI)
  MPI_Finalize();
#endif
  return rc;
}
