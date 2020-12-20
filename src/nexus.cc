/*
 * Copyright (c) 2017-2019, Carnegie Mellon University and
 *     Los Alamos National Laboratory.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
 * HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <assert.h>

#include "nexus_internal.h"

#define DEFAULT_NX_LIMIT 4   /* default max# pending hg addr lookup req's */

/*
 * nexus_bootstrap: bootstrap nexus library (collective call)
 */
nexus_ctx_t nexus_bootstrap(progressor_handle_t *nethand,
                            progressor_handle_t *localhand) {
  nexus_ctx_t nctx;
  char *env;
  progressor_handle_t *nasmhand = NULL;  /* used if localhand == NULL */
  hg_class_t *nasmcls = NULL;
  hg_context_t *nasmctx = NULL;

  /*
   * sanity check args: both nethand and localhand must be in listen mode
   */
  if (!HG_Class_is_listening(mercury_progressor_hgclass(nethand)) ||
      (localhand &&
       !HG_Class_is_listening(mercury_progressor_hgclass(localhand))) ) {
    fprintf(stderr, "nexus_bootstrap: error: handles not listening\n");
    return(NULL);
  }

  nctx = new nexus_ctx;   /* c++ malloc+init here */
  nctx->local2global = NULL;
  nctx->rank2node = NULL;
  nctx->node2rep = NULL;
  nctx->localcomm = MPI_COMM_NULL;
  nctx->repcomm = MPI_COMM_NULL;
  nctx->hg_remote = NULL;
  nctx->hg_local = NULL;
  nctx->internal_local = 0;
  /* now safe to call nx_destroy() on nctx if there is an error */

  env = getenv("NEXUS_LOOKUP_LIMIT");    /* allow env to override default */
  if (env == NULL || env[0] == 0) {
    nctx->nx_limit = DEFAULT_NX_LIMIT;
  } else {
    nctx->nx_limit = atoi(env);
    if (nctx->nx_limit <= 0) {
      nctx->nx_limit = 1;
    }
  }

  /*
   * if we are not given a local handle, we default to generating
   * a new na+sm one.  the app can access this using nexus_localprogressor().
   */
  if (localhand == NULL) {
    char naurl[32];
    snprintf(naurl, sizeof(naurl), "na+sm://%d/0", getpid());
    nasmcls = HG_Init(naurl, HG_TRUE);
    if (nasmcls)
      nasmctx = HG_Context_create(nasmcls);
    if (nasmcls && nasmctx)
      nasmhand = mercury_progressor_init(nasmcls, nasmctx);

    if (!nasmhand) {
      if (nasmctx) HG_Context_destroy(nasmctx);
      if (nasmcls) HG_Finalize(nasmcls);
      goto error;
    }
    nctx->internal_local = 1; /* if we create it, we have to finalize it */
  }

  /*
   * install our comm and do our MPI setup
   */
  nctx->mycomm = MPI_COMM_WORLD;   /* don't hardwire it, may change later */
  if (nx_mpisetup(nctx) < 0)
    goto error;

  /*
   * install progress handles
   */
  if ((nctx->hg_remote = mercury_progressor_duphandle(nethand)) == NULL) {
    fprintf(stderr, "nexus_bootstrap: error: duphandle nethand failed\n");
    goto error;
  }
  if (nasmhand) {
    nctx->hg_local = nasmhand;  /* transfer ownership to nctx */
    nasmhand = NULL;
  } else {
    nctx->hg_local = mercury_progressor_duphandle(localhand);
    if (nctx->hg_local == NULL) {
      fprintf(stderr, "nexus_bootstrap: error: duphandle localhand failed\n");
      goto error;
    }
  }

  /*
   * build lmap: maps local per global rank to peer's local address
   */
  if (nx_build_lmap(nctx) < 0)
    goto error;

  if (!nctx->grank)
    fprintf(stdout, "NX: LOCAL %s (NX-LIMIT=%d)\n",
            (mercury_progressor_hgcontext(nctx->hg_local) !=
             mercury_progressor_hgcontext(nctx->hg_remote)) ? "DONE"
                                                            : "VIA REMOTE",
             nctx->nx_limit);

  /*
   * build rmap: maps remote node number to its rep's remote address.
   */
  if (nx_build_rmap(nctx) < 0)
    goto error;
  if (!nctx->grank) fprintf(stdout, "NX: REMOTE DONE\n");

  /*
   * done!
   */
#ifdef NEXUS_DEBUG
  fprintf(stderr, "NX-%d: ALL DONE\n", nctx->grank);
#endif
  return(nctx);

error:
  if (nasmhand) mercury_progressor_freehandle(nasmhand);
  if (localhand == NULL && nctx->internal_local) {
    HG_Context_destroy(nasmctx);
    HG_Finalize(nasmcls);
  }
  nx_destroy(nctx, 0);
  return(NULL);
}

/*
 * nexus_destroy: shutdown and free a nexus
 */
void nexus_destroy(nexus_ctx_t nctx) {
  nx_destroy(nctx, 1);   /* heavy lifting done in internal function */
}

/*
 * nx_addr_get: private utility function for nexus_next_hop.  get address
 * or return HG_ADDR_NULL if no address found.
 */
namespace {
hg_addr_t nx_addr_get(nexus_map_t* map, int key) {
  nexus_map_t::iterator it = map->find(key);
  return (it == map->end()) ? HG_ADDR_NULL : it->second;
}
}  // namespace

/*
 * nexus_next_hop: lookup next hop info in nexus
 */
nexus_ret_t nexus_next_hop(nexus_ctx_t nctx, int dest, int* rank,
                           hg_addr_t* addr) {
  int srcrep, destrep;
  int destn;
  assert(nctx != NULL);

  /* stop here if we are the final hop */
  if (nctx->grank == dest) return NX_DONE;

  /* if dest is local we return its local address */
  if ((*addr = nx_addr_get(&nctx->lmap, dest)) != HG_ADDR_NULL) {
    *rank = dest; /* the next stop is the final dest */
    /* we are either the original src, or a dest rep */
    return NX_ISLOCAL;
  }

  /* dest >> dest node */
  destn = nctx->rank2node[dest];
  /* dest node >> src rep's global rank */
  srcrep = nctx->local2global[(destn % nctx->lsize)];
  /* dest node >> dest rep's global rank */
  destrep = nctx->node2rep[destn];
#ifdef NEXUS_DEBUG
  fprintf(stderr, "NX-%d: dest=%d, destnode=%d, srcrep=%d, destrep=%d\n",
          nctx->grank, dest, destn, srcrep, destrep);
#endif

  if (nctx->grank != srcrep) {
    *addr = nx_addr_get(&nctx->lmap, srcrep);
    if (*addr == HG_ADDR_NULL) return NX_NOTFOUND;
    *rank = srcrep; /* the next stop is src rep */
    /* we are the original src */
    return NX_SRCREP;
  } else {
    *addr = nx_addr_get(&nctx->rmap, destn);
    if (*addr == HG_ADDR_NULL) return NX_NOTFOUND;
    *rank = destrep; /* the next stop is dest rep */
    /* we are the src rep */
    return NX_DESTREP;
  }
}

nexus_ret_t nexus_global_barrier(nexus_ctx_t nctx) {
  assert(nctx != NULL);
  int rv = MPI_Barrier(nctx->mycomm);
  if (rv != MPI_SUCCESS) return NX_ERROR;
  return NX_SUCCESS;
}

int nexus_global_rank(nexus_ctx_t nctx) {
  assert(nctx != NULL);
  return nctx->grank;
}

int nexus_global_size(nexus_ctx_t nctx) {
  assert(nctx != NULL);
  return nctx->gsize;
}

int nexus_node_id(nexus_ctx_t nctx) {
  assert(nctx != NULL);
  return nctx->nodeid;
}

int nexus_total_nodes(nexus_ctx_t nctx) {
  assert(nctx != NULL);
  return nctx->nnodes;
}

nexus_ret_t nexus_local_barrier(nexus_ctx_t nctx) {
  assert(nctx != NULL);
  int rv = MPI_Barrier(nctx->localcomm);
  if (rv != MPI_SUCCESS) return NX_ERROR;
  return NX_SUCCESS;
}

int nexus_local_rank(nexus_ctx_t nctx) {
  assert(nctx != NULL);
  return nctx->lrank;
}

int nexus_local_size(nexus_ctx_t nctx) {
  assert(nctx != NULL);
  return nctx->lsize;
}

int nexus_local2global(nexus_ctx_t nctx, int lrank) {
  assert(nctx != NULL);
  return nctx-> local2global[lrank];
}

int nexus_rank2node(nexus_ctx_t nctx, int grank) {
  assert(nctx != NULL);
  return nctx->rank2node[grank];
}

int nexus_node2rep(nexus_ctx_t nctx, int nodeid) {
  assert(nctx != NULL);
  return nctx-> node2rep[nodeid];
}

nexus_ret_t nexus_set_grank(nexus_ctx_t nctx, int rank) {
  nctx->grank = rank;
  return NX_SUCCESS;
}

progressor_handle_t *nexus_localprogressor(nexus_ctx_t nctx) {
  assert(nctx != NULL);
  return nctx->hg_local;
}

progressor_handle_t *nexus_remoteprogressor(nexus_ctx_t nctx) {
  assert(nctx != NULL);
  return nctx->hg_remote;
}

/*
 * nexus_dump: dump nexus tables to stderr or files
 */
void nexus_dump(nexus_ctx_t nctx, char *outfile) {
  char *fname, *addr;
  int fnamelen, lcv;
  hg_size_t addr_alloc_sz, sz;
  FILE *fp;
  nexus_map_t::iterator it;
  hg_class_t *cls;
  hg_return_t hret;

  if (outfile) {
    fnamelen = strlen(outfile) + 32;  /* add room for suffix */
    fname = (char *)malloc(fnamelen);
    if (!fname) {
      fprintf(stderr, "nexus_dump: fname malloc failed\n");
      return;
    }
  } else {
    fname = NULL;
  }

  addr_alloc_sz = (nctx->gaddrsz > nctx->laddrsz) ? nctx->gaddrsz
                                                  : nctx->laddrsz;
  addr = (char *)malloc(addr_alloc_sz);
  if (!addr) {
    fprintf(stderr, "nexus_dump: addr malloc failed\n");
    if (fname) free(fname);
    return;
  }

  if (outfile) {
    snprintf(fname, fnamelen, "%s.%d.id", outfile, nctx->grank);
    fp = fopen(fname, "w");
    if (!fp) {
      perror("nexus_dump");
      goto done;
    }
  } else {
    fp = stderr;
  }

  fprintf(fp, "NX-%d: %d %d %d %d %d %d\n", nctx->grank, nctx->grank,
          nctx->gsize, nctx->lrank, nctx->lsize, nctx->nodeid,
          nctx->nnodes);
  fprintf(fp, "NX-%d: local %s\n", nctx->grank,
          mercury_progressor_addrstring(nctx->hg_local));
  fprintf(fp, "NX-%d: remote %s\n", nctx->grank,
          mercury_progressor_addrstring(nctx->hg_remote));
  fprintf(fp, "NX-%d: grank2node", nctx->grank);
  for (lcv = 0 ; lcv < nctx->gsize ; lcv++) {
    fprintf(fp, " %d", nctx->rank2node[lcv]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "NX-%d: local2global", nctx->grank);
  for (lcv = 0 ; lcv < nctx->lsize ; lcv++) {
    fprintf(fp, " %d", nctx->local2global[lcv]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "NX-%d: node2rep", nctx->grank);
  if (nctx->node2rep) {
    for (lcv = 0 ; lcv < nctx->nnodes ; lcv++) {
      fprintf(fp, " %d", nctx->node2rep[lcv]);
    }
  }
  fprintf(fp, "\n");

  if (fp != stderr)
    fclose(fp);

  if (outfile) {
    snprintf(fname, fnamelen, "%s.%d.lmap", outfile, nctx->grank);
    fp = fopen(fname, "w");
    if (!fp) {
      perror("nexus_dump");
      goto done;
    }
  } else {
    fp = stderr;
  }

  cls = mercury_progressor_hgclass(nctx->hg_local);
  it = nctx->lmap.begin();
  for (; it != nctx->lmap.end(); ++it) {
    sz = addr_alloc_sz;
    hret = HG_Addr_to_string(cls, addr, &sz, it->second);
    if (hret != HG_SUCCESS) strcpy(addr, "n/a");
    fprintf(fp, "NX-%d: lmap %d %s\n", nctx->grank, it->first, addr);
  }
  if (fp != stderr)
    fclose(fp);

  if (outfile) {
    snprintf(fname, fnamelen, "%s.%d.rmap", outfile, nctx->grank);
    fp = fopen(fname, "w");
    if (!fp) {
      perror("nexus_dump");
      goto done;
    }
  } else {
    fp = stderr;
  }

  cls = mercury_progressor_hgclass(nctx->hg_remote);
  it = nctx->rmap.begin();
  for (; it != nctx->rmap.end(); ++it) {
    sz = addr_alloc_sz;
    hret = HG_Addr_to_string(cls, addr, &sz, it->second);
    if (hret != HG_SUCCESS) strcpy(addr, "n/a");
    fprintf(fp, "NX-%d: rmap %d %s\n", nctx->grank, it->first, addr);
  }
  if (fp != stderr)
    fclose(fp);


done:
  if (fname) free(fname);
  if (addr) free(addr);
  return;
}
