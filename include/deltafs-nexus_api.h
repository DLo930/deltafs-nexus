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
#pragma once

#include <mercury.h>
#include <mercury-progressor/mercury-progressor.h>
#include <map>

typedef struct nexus_ctx* nexus_ctx_t;
typedef struct nexus_iter* nexus_iter_t;

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

#ifdef __cplusplus
extern "C" {
#endif

/**
 * nexus_bootstrap: bootstrap nexus library.  on success nexus will
 * dup the handles to keep a reference to mercury.  a collective call.
 *
 * @param nethand progressor handle for network (non-local) traffic
 * @param localhand progressor handle for local traffic (e.g. na+sm)
 * @return nexus context or NULL on error
 */
nexus_ctx_t nexus_bootstrap(progressor_handle_t *nethand,
                            progressor_handle_t *localhand);

/**
 * Destroys the Nexus library freeing all allocated resources, including
 * freeing its dup'd progressor_handle_t structures.
 *
 * @param nexus context
 */
void nexus_destroy(nexus_ctx_t nctx);

/**
 * Returns next Mercury address in route to dest or error
 * @param nexus context
 * @param MPI rank of destination
 * @param MPI rank of next hop (returned)
 * @param Mercury address of next hop (returned)
 * @return NX_SUCCESS or an error code
 */
nexus_ret_t nexus_next_hop(nexus_ctx_t nctx, int dest, int* rank,
                           hg_addr_t* addr);

/**
 * Blocks until all processes in the global communicator have reached this
 * routine.
 * @param nctx context
 * @return NX_SUCCESS or an error code
 */
nexus_ret_t nexus_global_barrier(nexus_ctx_t nctx);

/**
 * Return global rank of this process (assumes nexus is up)
 * @param nctx context
 * @return global rank
 */
int nexus_global_rank(nexus_ctx_t nctx);

/**
 * Return size of the global communicator (assumes nexus is up)
 * @param nctx context
 * @return size of the global communicator
 */
int nexus_global_size(nexus_ctx_t nctx);

/**
 * Return node id of this process (assumes nexus is up)
 * @param nctx context
 * @return node id
 */
int nexus_node_id(nexus_ctx_t nctx);

/**
 * Return number of nodes (assumes nexus is up)
 * @param nctx context
 * @return number of nodes
 */
int nexus_total_nodes(nexus_ctx_t nctx);

/**
 * Blocks until all processes in the local communicator have reached this
 * routine.
 * @param nctx context
 * @return NX_SUCCESS or an error code
 */
nexus_ret_t nexus_local_barrier(nexus_ctx_t nctx);

/**
 * Return local rank of this process (assumes nexus is up)
 * @param nctx context
 * @return local rank
 */
int nexus_local_rank(nexus_ctx_t nctx);

/**
 * Return size of the local communicator (assumes nexus is up)
 * @param nctx context
 * @return size of the global communicator
 */
int nexus_local_size(nexus_ctx_t nctx);

/**
 * Return global rank given a local rank (assumes nexus is up)
 * @param nctx context
 * @param lrank local rank
 * @return global rank of the given local rank
 */
int nexus_local2global(nexus_ctx_t nctx, int lrank);

/**
 * Return node id given a global rank (assumes nexus is up)
 * @param nctx context
 * @param grank global rank
 * @return node id of the given global rank
 */
int nexus_rank2node(nexus_ctx_t nctx, int grank);

/**
 * Return global rank of a given node's rep (assumes nexus is up)
 * @param nctx context
 * @param nodeid node id
 * @return global rank of the given node
 */
int nexus_node2rep(nexus_ctx_t nctx, int nodeid);

/**
 * Sets the global rank of the process (for debug purposes)
 * @param nexus context
 * @param new MPI rank
 * @return NX_SUCCESS or an error code
 */
nexus_ret_t nexus_set_grank(nexus_ctx_t nctx, int rank);

/**
 * Return nctx's progressor handle for local communication.
 * result valid until nexus_destroy() is called.  note that
 * calling bootstrap with a NULL "localhand" causes nexus to
 * generate an internal instance of na+sm for local RPCs.
 * this instance will be finalized when nexus_destroy() is
 * called (even if you dup the handle returned by nexus_localprogressor()).
 *
 * @param nctx context
 */
progressor_handle_t *nexus_localprogressor(nexus_ctx_t nctx);

/**
 * Return nctx's progressor handle for remote communication.
 * result valid until nexus_destroy() is called.
 *
 * @param nctx context
 */
progressor_handle_t *nexus_remoteprogressor(nexus_ctx_t nctx);

/**
 * Dump nexus tables (for debugging)
 *
 * @param nctx context
 * @param outfile output file (NULL means dump to stderr)
 */
void nexus_dump(nexus_ctx_t nctx, char *outfile);

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
void nexus_iter_free(nexus_iter_t* nitp);

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
 * Return current subrank of iterator.  subrank is 0 for local
 * maps and node number for remote maps.
 *
 * @param nit iterator handle
 */
int nexus_iter_subrank(nexus_iter_t nit);

#ifdef __cplusplus
}
#endif
