/*
 * Copyright (c) 2016-2017 Carnegie Mellon University.
 *
 * All rights reserved.
 *
 * Use of this source code is governed by a BSD-style license that can be
 * found in the LICENSE file. See the AUTHORS file for names of contributors.
 */

#pragma once

#include <map>
#include <mercury.h>

typedef struct nexus_ctx *nexus_ctx_t;
typedef struct nexus_iter *nexus_iter_t;

/* error codes */
typedef enum {
    NX_SUCCESS = 0, /* operation succeeded */
    NX_ERROR,       /* operation resulted in error */
    NX_NOTFOUND,    /* address not found */
    NX_ISLOCAL,     /* dest is local */
    NX_SRCREP,      /* dest is srcrep */
    NX_DESTREP,     /* dest is dstrep */
    NX_INVAL,       /* invalid parameter */
    NX_DONE,        /* already at destination */
} nexus_ret_t;

/**
 * nexus_bootstrap: bootstraps the Nexus library
 * @param string of the network subnet to be preferred for Mercury endpoints
 * @param string of the Mercury protocol plugin to be preferred
 * @return nexus context or NULL on error
 */
nexus_ctx_t nexus_bootstrap(char *subnet, char *proto);

/**
 * Destroys the Nexus library freeing all allocated resources
 * @param nexus context
 */
void nexus_destroy(nexus_ctx_t nctx);

/**
 * Returns next Mercury address in route to dest or error
 * @param nexus context
 * @param MPI rank of destination
 * @param MPI rank of next hop (returned iff address is not NULL)
 * @param Mercury address of next hop (returned)
 * @return NX_SUCCESS or an error code
 */
nexus_ret_t nexus_next_hop(nexus_ctx_t nctx, int dest,
                           int *rank, hg_addr_t *addr);

/**
 * Return global rank of this process (assumes nexus is up)
 * @param nctx context
 * @return global rank
 */
int nexus_global_rank(nexus_ctx_t nctx);

/**
 * Sets the global rank of the process (for debug purposes)
 * @param nexus context
 * @param new MPI rank
 * @return NX_SUCCESS or an error code
 */
nexus_ret_t nexus_set_grank(nexus_ctx_t nctx, int rank);

/**
 * Return mercury class for local (na+sm)
 *
 * @param nctx context
 */
hg_class_t *nexus_hgclass_local(nexus_ctx_t nctx);

/**
 * Return mercury class for remote (bmi+tcp or something like that)
 *
 * @param nctx context
 */
hg_class_t *nexus_hgclass_remote(nexus_ctx_t nctx);

/**
 * Return mercury context for local (bmi+tcp or something like that)
 *
 * @param nctx context
 */
hg_context_t *nexus_hgcontext_local(nexus_ctx_t nctx);

/**
 * Return mercury context for remote (bmi+tcp or something like that)
 *
 * @param nctx context
 */
hg_context_t *nexus_hgcontext_remote(nexus_ctx_t nctx);

/**
 * Allocate a new iterator.  nctx must remain active while iter is
 * allocated.  must free iterator when done.
 *
 * @param nctx context
 * @param local set non-zero if you want a local map iterator
 */
nexus_iter_t nexus_iter(nexus_ctx_t nctx, int local);

/**
 * Free a previously allocated iterator
 *
 * @param nitp pointer to iterator handle (we set to null)
 */
void nexus_iter_free(nexus_iter_t *nitp);

/**
 * Return non-zero if we are at the end of the map.
 *
 * @param nit iterator handle
 */
int nexus_iter_atend(nexus_iter_t nit);

/**
 * Advance the iterator
 *
 * @param nit iterator handle
 */
void nexus_iter_advance(nexus_iter_t nit);

/**
 * Return current hgaddr of iterator
 *
 * @param nit iterator handle
 */
hg_addr_t nexus_iter_addr(nexus_iter_t nit);

/**
 * Return current global rank of iterator
 *
 * @param nit iterator handle
 */
int nexus_iter_globalrank(nexus_iter_t nit);

/**
 * Return current subrank of iterator.  subrank is local rank for local
 * maps and node number for remote maps.
 *
 * @param nit iterator handle
 */
int nexus_iter_subrank(nexus_iter_t nit);
