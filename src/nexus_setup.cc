/*
 * Copyright (c) 2017, Carnegie Mellon University.
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

#include <ifaddrs.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <netdb.h>
#include <pthread.h>

#include "nexus_internal.h"

typedef struct {
    hg_context_t *hgctx;
    int bgdone;
} bgthread_dat_t;

/*
 * Network support pthread. Need to call progress to push the network and then
 * trigger to run the callback.
 */
static void *nexus_bgthread(void *arg)
{
    bgthread_dat_t *bgdat = (bgthread_dat_t *)arg;
    hg_return_t hret;

#ifdef NEXUS_DEBUG
    fprintf(stdout, "Network thread running\n");
#endif

    /* while (not done sending or not done recving */
    while (!bgdat->bgdone) {
        unsigned int count = 0;

        do {
            hret = HG_Trigger(bgdat->hgctx, 0, 1, &count);
        } while (hret == HG_SUCCESS && count);

        if (hret != HG_SUCCESS && hret != HG_TIMEOUT)
            msg_abort("nexus_bgthread: HG_Trigger failed");

        hret = HG_Progress(bgdat->hgctx, 100);
        if (hret != HG_SUCCESS && hret != HG_TIMEOUT)
            msg_abort("nexus_bgthread: HG_Progress failed");
    }

    return NULL;
}

/*
 * Put together the remote Mercury endpoint address from bootstrap parameters.
 * Writes the server URI into *uri on success. Aborts on error.
 */
static void prepare_addr(nexus_ctx_t *nctx, int minport, int maxport,
                         char *subnet, char *proto, char *uri)
{
    struct ifaddrs *ifaddr, *cur;
    int family, ret, rank, size, port;
    char ip[16];

    /* Query local socket layer to get our IP addr */
    if (getifaddrs(&ifaddr) == -1)
        msg_abort("getifaddrs failed");

    for (cur = ifaddr; cur != NULL; cur = cur->ifa_next) {
        if (cur->ifa_addr != NULL) {
            family = cur->ifa_addr->sa_family;

            if (family == AF_INET) {
                if (getnameinfo(cur->ifa_addr, sizeof(struct sockaddr_in), ip,
                                sizeof(ip), NULL, 0, NI_NUMERICHOST) == -1)
                    msg_abort("getnameinfo failed");

                if (strncmp(subnet, ip, strlen(subnet)) == 0)
                    break;
            }
        }
    }

    if (cur == NULL)
        msg_abort("no ip addr");

    freeifaddrs(ifaddr);

    /* sanity check on port range */
    if (maxport - minport < 0)
        msg_abort("bad min-max port");
    if (minport < 1)
        msg_abort("bad min port");
    if (maxport > 65535)
        msg_abort("bad max port");

    MPI_Comm_rank(nctx->localcomm, &rank);
    MPI_Comm_size(nctx->localcomm, &size);
    port = minport + (rank % (1 + maxport - minport));
    for (; port <= maxport; port += size) {
        int so, n = 1;
        struct sockaddr_in addr;

        /* test port availability */
        so = socket(PF_INET, SOCK_STREAM, 0);
        setsockopt(so, SOL_SOCKET, SO_REUSEADDR, &n, sizeof(n));
        if (so != -1) {
            addr.sin_family = AF_INET;
            addr.sin_addr.s_addr = INADDR_ANY;
            addr.sin_port = htons(port);
            n = bind(so, (struct sockaddr*)&addr, sizeof(addr));
            close(so);
            if (n == 0)
                break; /* done */
        } else {
            msg_abort("socket");
        }
    }

    if (port > maxport) {
        int so, n = 1;
        struct sockaddr_in addr;
        socklen_t addr_len;

        port = 0;
        fprintf(stderr, "Warning: no free ports available within the specified "
                "range\n>>> auto detecting ports ...\n");
        so = socket(PF_INET, SOCK_STREAM, 0);
        setsockopt(so, SOL_SOCKET, SO_REUSEADDR, &n, sizeof(n));
        if (so != -1) {
            addr.sin_family = AF_INET;
            addr.sin_addr.s_addr = INADDR_ANY;
            addr.sin_port = htons(0);
            n = bind(so, (struct sockaddr*)&addr, sizeof(addr));
            if (n == 0) {
                n = getsockname(so, (struct sockaddr*)&addr, &addr_len);
                if (n == 0)
                    port = ntohs(addr.sin_port); /* okay */
            }
            close(so);
        } else {
            msg_abort("socket");
        }
    }

    if (port == 0)
        msg_abort("no free ports");

    /* add proto */
    sprintf(uri, "%s://%s:%d", proto, ip, port);
#ifdef NEXUS_DEBUG
    fprintf(stdout, "Info: Using address %s\n", uri);
#endif
}

typedef struct hg_lookup_out {
    hg_return_t hret;
    hg_addr_t addr;
    pthread_mutex_t cb_mutex;
    pthread_cond_t cb_cv;
} hg_lookup_out_t;

static hg_return_t hg_lookup_cb(const struct hg_cb_info *info)
{
    hg_lookup_out_t *out = (hg_lookup_out_t *)info->arg;
    out->hret = info->ret;
    if (out->hret != HG_SUCCESS)
        out->addr = HG_ADDR_NULL;
    else
        out->addr = info->info.lookup.addr;

    pthread_mutex_lock(&out->cb_mutex);
    pthread_cond_signal(&out->cb_cv);
    pthread_mutex_unlock(&out->cb_mutex);

    return HG_SUCCESS;
}

static hg_return_t hg_lookup(nexus_ctx_t *nctx, hg_context_t *hgctx,
                             char *hgaddr, hg_addr_t *addr)
{
    hg_lookup_out_t *out = NULL;
    hg_return_t hret;

    /* Init addr metadata */
    out = (hg_lookup_out_t *)malloc(sizeof(*out));
    if (out == NULL)
        return HG_NOMEM_ERROR;

    /* rank is set, perform lookup */
    pthread_mutex_init(&out->cb_mutex, NULL);
    pthread_cond_init(&out->cb_cv, NULL);
    pthread_mutex_lock(&out->cb_mutex);
    hret = HG_Addr_lookup(hgctx, &hg_lookup_cb, out, hgaddr, HG_OP_ID_IGNORE);
    if (hret != HG_SUCCESS)
        goto err;

    /* Lookup posted, wait until finished */
    pthread_cond_wait(&out->cb_cv, &out->cb_mutex);
    pthread_mutex_unlock(&out->cb_mutex);

    if (out->hret != HG_SUCCESS) {
        hret = out->hret;
    } else {
        hret = HG_SUCCESS;
        *addr = out->addr;
    }

err:
    pthread_cond_destroy(&out->cb_cv);
    pthread_mutex_destroy(&out->cb_mutex);
    free(out);
    return hret;
}

typedef struct {
    int pid;
    int hgid;
    int grank;
    int lrank;
} ldata_t;

static void discover_local_info(nexus_ctx_t *nctx)
{
    int ret;
    char hgaddr[128];
    ldata_t ldat;
    ldata_t *hginfo;
    hg_return_t hret;
    pthread_t bgthread; /* network background thread */
    bgthread_dat_t *bgarg;

    MPI_Comm_rank(nctx->localcomm, &(nctx->lrank));
    MPI_Comm_size(nctx->localcomm, &(nctx->lsize));

    /* Initialize local Mercury listening endpoints */
    snprintf(hgaddr, sizeof(hgaddr), "na+sm://%d/0", getpid());
#ifdef NEXUS_DEBUG
    fprintf(stderr, "Initializing for %s\n", hgaddr);
#endif

    nctx->local_hgcl = HG_Init(hgaddr, HG_TRUE);
    if (!nctx->local_hgcl)
        msg_abort("HG_init failed for local endpoint");

    nctx->local_hgctx = HG_Context_create(nctx->local_hgcl);
    if (!nctx->local_hgctx)
        msg_abort("HG_Context_create failed for local endpoint");

    /* Start the network thread */
    bgarg = (bgthread_dat_t *)malloc(sizeof(*bgarg));
    if (!bgarg)
        msg_abort("malloc failed");

    bgarg->hgctx = nctx->local_hgctx;
    bgarg->bgdone = 0;

    ret = pthread_create(&bgthread, NULL, nexus_bgthread, (void*)bgarg);
    if (ret != 0)
        msg_abort("pthread_create failed");

    /* Exchange PID, ID, global rank, local rank among local ranks */
    ldat.pid = getpid();
    ldat.hgid = 0;
    ldat.grank = nctx->grank;
    ldat.lrank = nctx->lrank;

    hginfo = (ldata_t *)malloc(sizeof(ldata_t) * (nctx->lsize));
    if (!hginfo)
        msg_abort("malloc failed");

    MPI_Allgather(&ldat, sizeof(ldata_t), MPI_BYTE, hginfo,
                  sizeof(ldata_t), MPI_BYTE, nctx->localcomm);

    /* Build local => global rank map */
    nctx->localranks = (int *)malloc(sizeof(int) * (nctx->lsize));
    if (!nctx->localranks)
        msg_abort("malloc failed");

    for (int i = 0; i < nctx->lsize; i++) {
        int eff_i = (nctx->lrank + i) % nctx->lsize;
        hg_addr_t localaddr;

        /* Find the local root */
        if (hginfo[eff_i].lrank == 0)
            nctx->lroot = hginfo[eff_i].grank;

        /* Update mapping */
        nctx->localranks[hginfo[eff_i].lrank] = hginfo[eff_i].grank;

#ifdef NEXUS_DEBUG
        fprintf(stdout, "[%d] Idx %d: pid %d, id %d, grank %d, lrank %d\n",
                nctx->grank, eff_i, hginfo[eff_i].pid, hginfo[eff_i].hgid,
                hginfo[eff_i].grank, hginfo[eff_i].lrank);
#endif

        snprintf(hgaddr, sizeof(hgaddr), "na+sm://%d/%d",
                 hginfo[eff_i].pid, hginfo[eff_i].hgid);

        if (hginfo[eff_i].grank == nctx->grank) {
            hret = HG_Addr_self(nctx->local_hgcl, &localaddr);
        } else {
            hret = hg_lookup(nctx, nctx->local_hgctx, hgaddr, &localaddr);
        }

        if (hret != HG_SUCCESS) {
            fprintf(stderr, "Tried to lookup %s\n", hgaddr);
            msg_abort("hg_lookup failed");
        }

        /* Add to local map */
        nctx->laddrs[hginfo[eff_i].grank] = localaddr;
#ifdef NEXUS_DEBUG
        print_hg_addr(nctx->local_hgcl, hgaddr, localaddr);
#endif
    }

    free(hginfo);

    /* Sync before terminating background threads */
    MPI_Barrier(nctx->localcomm);

    /* Terminate network thread */
    bgarg->bgdone = 1;
    pthread_join(bgthread, NULL);

    free(bgarg);
}

typedef struct {
    char addr[60];
    int grank;
} rdata_t;

static void discover_remote_info(nexus_ctx_t *nctx, char *hgaddr)
{
    int ret, rep_size, rep_rank;
    hg_return_t hret;
    pthread_t bgthread; /* network background thread */
    bgthread_dat_t *bgarg;
    rdata_t rdat;
    rdata_t *hginfo;

    MPI_Comm_rank(nctx->repcomm, &rep_rank);
    MPI_Comm_size(nctx->repcomm, &rep_size);

    /* Build rank => rep mapping */
    nctx->rankreps = (int *)malloc(sizeof(int) * nctx->gsize);
    if (!nctx->rankreps)
        msg_abort("malloc failed");

    MPI_Allgather(&nctx->lroot, 1, MPI_INT, nctx->rankreps,
                  1, MPI_INT, MPI_COMM_WORLD);

    /* Initialize remote Mercury listening endpoints */
    nctx->remote_hgcl = HG_Init(hgaddr, HG_TRUE);
    if (!nctx->remote_hgcl)
        msg_abort("HG_Init failed for remote endpoint");

    nctx->remote_hgctx = HG_Context_create(nctx->remote_hgcl);
    if (!nctx->remote_hgctx)
        msg_abort("HG_Context_create failed for remote endpoint");

    /* Start the network thread */
    bgarg = (bgthread_dat_t *)malloc(sizeof(*bgarg));
    if (!bgarg)
        msg_abort("malloc failed");

    bgarg->hgctx = nctx->remote_hgctx;
    bgarg->bgdone = 0;

    ret = pthread_create(&bgthread, NULL, nexus_bgthread, (void*)bgarg);
    if (ret != 0)
        msg_abort("pthread_create failed");

    /* Exchange PID, ID, global rank, local rank among local ranks */
    strncpy(rdat.addr, hgaddr, sizeof(rdat.addr));
    rdat.addr[60] = '\0';
    rdat.grank = nctx->grank;

    hginfo = (rdata_t *)malloc(sizeof(rdata_t) * rep_size);
    if (!hginfo)
        msg_abort("malloc failed");

    MPI_Allgather(&rdat, sizeof(rdata_t), MPI_BYTE, hginfo,
                  sizeof(rdata_t), MPI_BYTE, nctx->repcomm);

    for (int i = 0; i < rep_size; i++) {
        int eff_i = (rep_rank + i) % rep_size;
        hg_addr_t remoteaddr;

#ifdef NEXUS_DEBUG
        fprintf(stdout, "[%d] eff_i = %d, addr = %s, grank = %d\n",
                i, eff_i, hginfo[eff_i].addr, hginfo[eff_i].grank);
#endif

        if (hginfo[eff_i].grank == nctx->grank) {
            hret = HG_Addr_self(nctx->remote_hgcl, &remoteaddr);
        } else {
            hret = hg_lookup(nctx, nctx->remote_hgctx, hginfo[eff_i].addr,
                             &remoteaddr);
        }

        if (hret != HG_SUCCESS) {
            fprintf(stderr, "Tried to lookup %s\n", hginfo[eff_i].addr);
            msg_abort("hg_lookup failed");
        }

        /* Add to remote map */
        nctx->gaddrs[hginfo[eff_i].grank] = remoteaddr;
    }

    /* Sync before terminating background threads */
    MPI_Barrier(nctx->repcomm);

    /* Terminate network thread */
    bgarg->bgdone = 1;
    pthread_join(bgthread, NULL);

    free(hginfo);
    free(bgarg);
}

nexus_ret_t nexus_bootstrap(nexus_ctx_t *nctx, int minport, int maxport,
                            char *subnet, char *proto)
{
    char hgaddr[128];

    /* Grab MPI rank info */
    MPI_Comm_rank(MPI_COMM_WORLD, &(nctx->grank));
    MPI_Comm_size(MPI_COMM_WORLD, &(nctx->gsize));

    if (!nctx->grank)
        fprintf(stdout, "Nexus: started bootstrap\n");

    init_local_comm(nctx);
    discover_local_info(nctx);

    if (!nctx->grank)
        fprintf(stdout, "Nexus: done local info discovery\n");

    prepare_addr(nctx, minport, maxport, subnet, proto, hgaddr);
    init_rep_comm(nctx);
    discover_remote_info(nctx, hgaddr);

    if (!nctx->grank)
        fprintf(stdout, "Nexus: done remote info discovery\n");

#ifdef NEXUS_DEBUG
    fprintf(stdout, "[%d] grank = %d, lrank = %d, gsize = %d, lsize = %d\n",
            nctx->grank, nctx->grank, nctx->lrank, nctx->gsize, nctx->lsize);
#endif /* NEXUS_DEBUG */

    return NX_SUCCESS;
}

nexus_ret_t nexus_destroy(nexus_ctx_t *nctx)
{
    map<int, hg_addr_t>::iterator it;

    /* Free local Mercury addresses */
    for (it = nctx->laddrs.begin(); it != nctx->laddrs.end(); it++)
        if (it->second != HG_ADDR_NULL)
            HG_Addr_free(nctx->local_hgcl, it->second);

    /* Sync before tearing down local endpoints */
    MPI_Barrier(nctx->localcomm);
    MPI_Comm_free(&nctx->localcomm);

    /* Destroy Mercury local endpoints */
    HG_Context_destroy(nctx->local_hgctx);
    HG_Finalize(nctx->local_hgcl);

    if (!nctx->grank)
        fprintf(stdout, "Nexus: done local info cleanup\n");

    /* Free remote Mercury addresses */
    for (it = nctx->gaddrs.begin(); it != nctx->gaddrs.end(); it++)
        if (it->second != HG_ADDR_NULL)
            HG_Addr_free(nctx->remote_hgcl, it->second);

    /* Sync before tearing down remote endpoints */
    MPI_Barrier(nctx->repcomm);

    /* Destroy Mercury remote endpoints */
    HG_Context_destroy(nctx->remote_hgctx);
    HG_Finalize(nctx->remote_hgcl);

    if (!nctx->grank)
        fprintf(stdout, "Nexus: done remote info cleanup\n");

    free(nctx->rankreps);
    MPI_Comm_free(&nctx->repcomm);
    return NX_SUCCESS;
}
