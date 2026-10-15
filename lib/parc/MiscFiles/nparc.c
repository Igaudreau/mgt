/* Parc.c:  parallel array support library (for Fortran)
*/
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"
#include "dam.h"
#include "parc.h"

/* Implementation layer
*/
/* Provide a default implementation (if undefined by compiler) */
#ifdef PARC_SMP
#undef PARC_MPI
#else
#ifdef PARC_MPI
#undef PARC_SMP
#else
#define PARC_SMP
#endif
#endif
/* Specify all sources? */
#define RECVSEQ
/* DUAL message structure (head,tail) */
#define PARC_DUAL
#undef BIG_HEAD_SIZE
#ifdef BIG_HEAD_SIZE
#define PARC_HEAD_SIZE 512
#else
#define PARC_HEAD_SIZE 128
#endif
#define PARC_HEAD_TAG 1
#define PARC_TAIL_TAG 2
/* MPI implementation */
#ifdef PARC_MPI
#include "mpi.h"
#define PARC_MPI_FATAL
#ifdef PARC_DUAL
#define PARC_MPI_LITE
#else
#undef PARC_MPI_LITE
#endif
#endif
/* SMP implementation
*/
#ifdef PARC_SMP
#include "smp.h"
#endif

/* Constants :
*/
#define PARC_EMPTY 0
#define MAXDIST 8
#define MAXCK 128
#define MAXLEV 31
#define MAXCACHE 256
#define MAXDIM 7
#define TREE_NULL -1

/* Data :
*/
static int Argc=0;
static char **Argv;
static bool alive=FALSE,dead=FALSE;
static int limpe=1,mype=0;
static DamStruct texdam,indam[MAXPE],outdam[MAXPE];
#ifdef PARC_DUAL
static uword tmphead[PARC_HEAD_SIZE+1000];
#endif
static int ntex=0;
static int pestatus=0;
static int nsent=0,nrecd=0;
static uword tv_sec0=0;
static bool notv_sec0=TRUE;
static int ck_count=0,ck_lev=0;
static char ck_names[MAXCK][SZNAM];
static float ck_cpu0[MAXCK],ck_cpu[MAXCK],ck_clock0[MAXCK],
  ck_clock[MAXCK];
static int ck_calls[MAXCK],ck_sent[MAXCK],ck_recvd[MAXCK],
  ck_sent0[MAXCK],ck_recvd0[MAXCK];
static int ncache=0;
static CacheStruct caches[MAXCACHE];
/*
  Generic timers
*/
/* parc_clock:  wall clock time (from system) */
void parc_clock(float *t) {
   struct timeval f;

   gettimeofday(&f, (struct timezone *)NULL);
/* Capture suitable offset to avoid loss of floating precision */
   if (notv_sec0) {
      tv_sec0 = (uword)(f.tv_sec);
      notv_sec0 = FALSE;
   }
   *t = ((float)(f.tv_sec-tv_sec0)) + ((float)f.tv_usec)*0.000001;
}
#ifdef PARC_SMP
/*
  Private SMP specific functions.
*/
/* parc_mp_exit:  tidy up and leave MP (SMP version) */
static void parc_mp_exit() {
   smp_exit();
}
/* parc_mp_register:  register with MP (SMP version) */
static bool parc_mp_register(int *mype, int *limpe, bool *error) {
   int nstep,n,i;
/* Return in error state */
   if (*error) {
      return(!(*error));
   }
/* Initialise SMP */
   if (!smp_register(Argc, Argv, limpe, mype)) {
      fprintf(stderr, "error:  smp_register failed (parc_mp_register)\n");
   } else {
   }
   return(!(*error));
}
/* parc_mp_sends:  implementation of sends (SMP version)  */
static void parc_mp_sends(int nsend, int sends[], bool *error) {
   int i,ipe;
   DamPtr s;
   int nword,nhead,ntail;
/* For each desitination */
   for (i=0; i < nsend; i++) {
      ipe = sends[i];
      s = &outdam[ipe];
/*    Load buffer length into first word of buffer */
      s->mb[0] = s->mp;
/*    Head and tail */
      nword = s->mb[0];
      if (nword > PARC_HEAD_SIZE) {
         nhead = PARC_HEAD_SIZE;
         ntail = nword - nhead;
      } else {
         nhead = PARC_HEAD_SIZE;
         ntail = 0;
      }
/*    Accounting statistics */
      nsent = nsent + nhead + ntail;
/*    Send buffer as a fixed-length header and trailer */
      if (!smp_send(ipe, PARC_HEAD_SIZE, s->mb)) {
         fprintf(stderr, "error[%d]:  smp_send failed\n", mype);
/*    No second message needed */
      } else if (ntail <= 0) {
/*    Second message */
      } else if (!smp_send(ipe, ntail, &(s->mb[nhead]))) {
         fprintf(stderr, "error[%d]:  smp_send failed\n", mype);
      } else {
      }
   }
}
/* parc_mp_recvs:  implementation of receives (SMP version) */
static void parc_mp_recvs(int nrecv, int recvs[], bool *error) {
   int ipe,i,j,nword,mpi_ios;
   bool found;
   DamPtr s;
   int nhead,ntail,nitem;
   for (i=0; i < nrecv; i++) {
/*    Receive any (header) message into temporary buffer */
      ipe = recvs[i];
      if (!smp_recv(ipe, PARC_HEAD_SIZE, tmphead)) {
         fprintf(stderr, "error[%d]:  smp_recv failed\n", mype);
/*    Find source and ultimate message size */
      } else {
/*       Ultimate message size */
         nword = (int)tmphead[0];
/*       Number of words in header */
         if (nword > PARC_HEAD_SIZE) {
            nhead = PARC_HEAD_SIZE;
            ntail = nword - nhead;
         } else {
            nhead = PARC_HEAD_SIZE;
            ntail = 0;
         }
/*       Transfer data from temporary buffer to input buffer */
         s = &indam[ipe];
         s->mp = 0;
         for (j=0; j < nhead; j++) {
            s->mb[s->mp++] = tmphead[j];
         }
         nrecd = nrecd + nhead + ntail;
/*       Optional second message */
         if (ntail <= 0) {
         } else if (!dam_space(nhead+ntail, s, error)) {
            fprintf(stderr,"parc_mp_recvs:  dam_space failed\n");
         } else if (!smp_recv(ipe, ntail, &(s->mb[nhead]))) {
            fprintf(stderr, "error[%d]:  smp_recv failed\n", mype);
         } else {
         }
      }
/*    Prepare buffer for read */
      s->mp = 0;
      s->mn = nword;
   }
}
/* parc_mp_barrier (SMP version):  barrier */
void parc_mp_barrier(int *ios) {
   smp_barrier();
}
#endif

#ifdef PARC_MPI
/*
  Private MPI specific functions.
*/
/* parc_mp_perror:  print MP error message (MPI version) */
static void parc_mp_perror(int mpi_ios) {
#ifndef PARC_MPI_LITE
   char buf[MPI_MAX_ERROR_STRING];
   int nchar;
#endif
#ifdef PARC_MPI_LITE
   fprintf(stderr, "error[%d]:  cannot print MPI error string\n",
     mype);
#else
   if (MPI_Error_string(mpi_ios, buf, &nchar) == MPI_SUCCESS) {
      fprintf(stderr,"mpi[%d]:  %s\n", mype, buf);
   } else {
      fprintf(stderr, "error[%d]:  failed to print MPI error\n",
        mype);
   }
#endif
}
/* parc_mp_exit:  tidy up and leave MP (MPI version) */
static void parc_mp_exit() {
   MPI_Finalize();
}

/* parc_mp_register:  register with MP (MPI version) */
static bool parc_mp_register(int *mype, int *limpe, bool *error) {
   int mpi_ios;
   int nstep,n;
/* Return in error state */
   if (*error) {
      return(!(*error));
   }
/* Initialise MPI */
   if ((mpi_ios=MPI_Init(&Argc, &Argv)) != MPI_SUCCESS) {
      fprintf(stderr, "error:  MPI_Init failed (parc_mp_register)\n");
      parc_mp_perror(mpi_ios);
      *error = TRUE;
#ifndef PARC_MPI_LITE
/* Set MPI error handler characteristics */
   } else if ((mpi_ios=MPI_Errhandler_set(MPI_COMM_WORLD,
     MPI_ERRORS_ARE_FATAL)) != MPI_SUCCESS) {
      fprintf(stderr, "error:  setting MPI error handler\n");
      parc_mp_perror(mpi_ios);
      *error = TRUE;
#endif
   } else if ((mpi_ios=MPI_Comm_size(MPI_COMM_WORLD, limpe)) !=
      MPI_SUCCESS) {
      fprintf(stderr, "error:  MPI_Comm_size failed (parc_register)\n");
      parc_mp_perror(mpi_ios);
      *error = TRUE;
   } else if ((mpi_ios=MPI_Comm_rank(MPI_COMM_WORLD, mype)) !=
       MPI_SUCCESS) {
      fprintf(stderr, "error:  MPI_Comm_rank failed (parc_register)\n");
      parc_mp_perror(mpi_ios);
      *error = TRUE;
   } else {
   }
   return(!(*error));
}
/* parc_mp_sends:  implementation of sends (MPI version)  */
static void parc_mp_sends(int nsend, int sends[], bool *error) {
   int i,ipe,mpi_ios;
   DamPtr s;
#ifdef PARC_DUAL
   int nword,nhead,ntail;
#endif

/* For each desitination */
   for (i=0; i < nsend; i++) {
      ipe = sends[i];
      s = &outdam[ipe];
#ifdef PARC_DUAL
/*    Load buffer length into first word of buffer */
      s->mb[0] = s->mp;
/*    Head and tail */
      nword = s->mb[0];
      if (nword > PARC_HEAD_SIZE) {
         nhead = PARC_HEAD_SIZE;
         ntail = nword - nhead;
      } else {
         nhead = PARC_HEAD_SIZE;
         ntail = 0;
      }
/*    Accounting statistics */
      nsent = nsent + nhead + ntail;
/*    Send buffer as a fixed-length header and trailer */
      if ((mpi_ios=MPI_Send(s->mb, PARC_HEAD_SIZE, mpi_uword, ipe, 
        PARC_HEAD_TAG, MPI_COMM_WORLD)) != MPI_SUCCESS) {
         fprintf(stderr, "error[%d]:  MPI_Sends failed\n", mype);
         parc_mp_perror(mpi_ios);
         *error = TRUE;
/*    No second message needed */
      } else if (ntail <= 0) {
/*    Second message */
      } else if ((mpi_ios=MPI_Send(&(s->mb[nhead]), ntail,
          mpi_uword, ipe, PARC_TAIL_TAG, 
          MPI_COMM_WORLD)) != MPI_SUCCESS) {
         fprintf(stderr, "error[%d]:  MPI_Sends failed\n", mype);
         parc_mp_perror(mpi_ios);
         *error = TRUE;
      } else {
      }
#else
/*    Send buffer in one message */
/*    Accounting statistics */
      nsent = nsent + s->mp;
      if ((mpi_ios=MPI_Send(s->mb, s->mp, mpi_uword, ipe, mype, 
        MPI_COMM_WORLD)) != MPI_SUCCESS) {
         fprintf(stderr, "error[%d]:  MPI_Sends failed\n", mype);
         parc_mp_perror(mpi_ios);
         *error = TRUE;
      }
#endif
   }
}
/* parc_mp_recvs:  implementation of receives (MPI version) */
static void parc_mp_recvs(int nrecv, int recvs[], bool *error) {
   int ipe,i,j,nword,mpi_ios;
   bool found;
   DamPtr s;
#ifdef PARC_DUAL
   int nhead,ntail,nitem;
#endif
   MPI_Status status;
   for (i=0; i < nrecv; i++) {
#ifdef PARC_DUAL
/*    Receive any (header) message into temporary buffer */
#ifdef RECVSEQ
      ipe = recvs[i];
#else
      ipe = MPI_ANY_SOURCE;
#endif
      if ((mpi_ios=MPI_Recv(tmphead, PARC_HEAD_SIZE, mpi_uword, 
        ipe,  PARC_HEAD_TAG,
        MPI_COMM_WORLD, &status)) != MPI_SUCCESS) {
         fprintf(stderr, "error[%d]:  MPI_Recv failed\n", mype);
         parc_mp_perror(mpi_ios);
         *error = TRUE;
/*    Find source and ultimate message size */
      } else {
/*       Source PE */
#ifndef RECVSEQ
         ipe = status.MPI_SOURCE;
#endif
         MPI_Get_count(&status, mpi_uword, &nitem);
/*       Ultimate message size */
         nword = (int)tmphead[0];
/*       Number of words in header */
         if (nword > PARC_HEAD_SIZE) {
            nhead = PARC_HEAD_SIZE;
            ntail = nword - nhead;
         } else {
            nhead = PARC_HEAD_SIZE;
            ntail = 0;
         }
/*       Transfer data from temporary buffer to input buffer */
         s = &indam[ipe];
         s->mp = 0;
         for (j=0; j < nhead; j++) {
            s->mb[s->mp++] = tmphead[j];
         }
         nrecd = nrecd + nhead + ntail;
         recvs[i] = ipe;
/*       Optional second message */
         if (ntail <= 0) {
         } else if (!dam_space(nhead+ntail, s, error)) {
            fprintf(stderr,"parc_mp_recvs:  dam_space failed\n");
         } else if ((mpi_ios=MPI_Recv(&(s->mb[nhead]), ntail, 
           mpi_uword, ipe,  PARC_TAIL_TAG,
           MPI_COMM_WORLD, &status)) != MPI_SUCCESS) {
            fprintf(stderr, "error[%d]:  MPI_Recv failed\n", mype);
            parc_mp_perror(mpi_ios);
            *error = TRUE;
         } else {
         }
      }
#else
/*    Probe then receive message in one go */
/*    Probe for any available source */
      if ((mpi_ios=MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, 
        MPI_COMM_WORLD, &status)) != MPI_SUCCESS) {
         fprintf(stderr, "error[%d]:  MPI_Probe failed\n", mype);
         parc_mp_perror(mpi_ios);
         *error = TRUE;
      }
/*    Get size of source message */
      if ((mpi_ios=MPI_Get_count(&status, mpi_uword, &nword))
         != MPI_SUCCESS) {
         fprintf(stderr, "error[%d]:  MPI_Get_count failed\n", mype);
         parc_mp_perror(mpi_ios);
         *error = TRUE;
      }
/*    Get source PE */
      ipe = status.MPI_SOURCE;
      nrecd = nrecd + nword;
      recvs[i] = ipe;
      s = &indam[ipe];
      s->mp = 0;
      s->mn = 0;
      if (!dam_space(nword, s, error)) {
         fprintf(stderr, "error[%d]:  failed to allocate memory\n",
           mype);
         *error = TRUE;
      } else if ((mpi_ios=MPI_Recv(s->mb, nword, mpi_uword, ipe, 
        MPI_ANY_TAG, MPI_COMM_WORLD, &status)) != MPI_SUCCESS) {
         fprintf(stderr, "error[%d]:  MPI_Recv failed\n", mype);
         parc_mp_perror(mpi_ios);
         *error = TRUE;
      }
#endif
/*    Prepare buffer for read */
      s->mp = 0;
      s->mn = nword;
   }
}
/* parc_mp_barrier:  barrier (MPI version) */
void parc_mp_barrier(int *ios) {
   int mpi_ios;
   if ((mpi_ios=MPI_Barrier(MPI_COMM_WORLD)) != MPI_SUCCESS) {
      *ios = 0;
   } else {
      *ios = 1;
   }
}
#endif
/*
   Private generic functions
*/
/* parc_sr:  recursive utility for send-recv */
static void parc_sr(int offset, int count,
  int nsend, int sends[], bool sendok[],
  int nrecv, int recvs[], bool recvok[],
  bool *error) {
   int nleft,tide,nright,nsendx,sendxs[MAXPE],nrecvx,recvxs[MAXPE],
     isend,ipe,irecv;
/* Split the PE range into two */
/* count for left-hand part */
   nleft = count/2;
/* Offset for right-hand part */
   tide = offset + nleft;
/* count for right hand part */
   nright = count - nleft;
/* LEFT numbered PEs send, RIGHT numbered PEs receive */
/* We are LEFT number PE */
   if ((mype >= offset) && (mype < (offset+nleft))) {
      nsendx = 0;
      for (isend=0; isend < nsend; isend++) {
         ipe = sends[isend];
         if ((ipe >= tide) && (ipe < (tide+nright))) {
            sendxs[nsendx++] = sends[isend];
            sendok[isend] = TRUE;
         }
      }
      parc_mp_sends(nsendx, sendxs, error);
/* We are RIGHT numbered PE */
   } else {
      nrecvx = 0;
      for (irecv=0; irecv < nrecv; irecv++) {
         ipe = recvs[irecv];
         if ((ipe >= offset) && (ipe < (offset+nleft))) {
#ifdef RECVSEQ
            recvxs[nrecvx++] = ipe;
#else
            nrecvx++;
#endif
            recvok[irecv] = TRUE;
         }
      }
      parc_mp_recvs(nrecvx, recvxs, error);
   }
/* RIGHT numbered PEs send, LEFT numbered PEs receive */
/* We are RIGHT number PE */
   if ((mype >= tide) && (mype < (tide+nright))) {
      nsendx = 0;
      for (isend=0; isend < nsend; isend++) {
         ipe = sends[isend];
         if ((ipe >= offset) && (ipe < (offset+nleft))) {
            sendxs[nsendx++] = ipe;
            sendok[isend] = TRUE;
         }
      }
      parc_mp_sends(nsendx, sendxs, error);
/* We are LEFT numbered PE */
   } else {
      nrecvx = 0;
      for (irecv=0; irecv < nrecv; irecv++) {
         ipe = recvs[irecv];
         if ((ipe >= tide) && (ipe < (tide+nright))) {
#ifdef RECVSEQ
            recvxs[nrecvx++] = recvs[irecv];
#else
            nrecvx++;
#endif
            recvok[irecv] = TRUE;
         }
      }
      parc_mp_recvs(nrecvx, recvxs, error);
   }
/* Recurse only if MYPE is contained in either part */
   if (nleft < 2) {
   } else if (mype < offset) {
   } else if (mype-offset >= nleft) {
   } else {
      parc_sr(offset, nleft, nsend, sends, sendok, nrecv, recvs, recvok,
        error);
   }
   if (nright < 2) {
   } else if (mype < tide) {
   } else if (mype-tide >= nright) {
   } else {
      parc_sr(tide, nright, nsend, sends, sendok, nrecv, recvs, recvok,
        error);
   }
}
/* parc_sendrecv:  implementation of senrecv (generic,improved) */
static void parc_sendrecv(int nsend, int sends[], 
  int nrecv, int recvs[], bool *error) {
   int count,offset,i;
   bool sendok[MAXPE],recvok[MAXPE];
/* Blank out check arrays */
   for (i=0; i < nsend; i++) {
      sendok[i] = FALSE;
   }
   for (i=0; i < nrecv; i++) {
      recvok[i] = FALSE;
   }
/* if more than one PE, recurse */
   offset = 0;
   count = limpe;
   if (count > 1) {
      parc_sr(offset, count, nsend, sends, sendok, nrecv, recvs,
        recvok, error);
   }
/* Check that we have send and recieved what we expected */
   for (i=0; i < nsend; i++) {
      if (!sendok[i]) {
         fprintf(stderr, "error[%d]:  did not send to %d (sendrecv)\n",
           mype, i);
         *error = TRUE;
      }
   }
   for (i=0; i < nrecv; i++) {
      if (!recvok[i]) {
         fprintf(stderr, "error[%d]:  did not recv from %d (sendrecv)\n",
           mype, i);
         *error = TRUE;
      }
   }
}
/* parc_tex2buf:  transfer local text to given buffer */
static bool parc_tex2buf(DamPtr t, int status, int nt, DamPtr d, bool *error) {
   int nchar,nword,it,i;
/* Rewind buffer to receive text */
   d->mp = 0;
   d->mn = 0;
/* Rewind text to read */
   t->mp = 0;
/* Number of words in new buffer */
#ifdef PARC_DUAL
   nword = t->mn + 3;
#else
   nword = t->mn + 2;
#endif
   if (dam_space(nword, d, error)) {
#ifdef PARC_DUAL
/*    Blank out location for buffer length */
      d->ib[d->mp++] = PARC_EMPTY;
#endif
/*    Store status */
      d->ib[d->mp++] = status;
/*    Store number of text messages */
      d->ib[d->mp++] = nt;
/*    For each text message */
      for (it = 0; it < nt; it++) {
         nchar = t->ib[t->mp++];
         d->ib[d->mp++] = nchar;
         for (i=0; i < nchar; i++) {
            d->ib[d->mp++] = t->ib[t->mp++];
         }
      }
   }
   return(!(*error));
}
/* parc_prtex:  print merged texts from siblings */
static bool parc_prtex(bool *error) {
   int ipe,i,nchar,it,nt,status;
   DamPtr s;
/* For each buffer */
   for (ipe=0; ipe < limpe; ipe++) {
      s = &indam[ipe];
/*    Rewind source buffer for read */
      s->mp = 0;
#ifdef PARC_DUAL
/*    Number of words in buffer (ignore) */
      i = s->ib[s->mp++];
#endif
/*    Status */
      status = s->ib[s->mp++];
/*    Number of text messages */
      nt = s->ib[s->mp++];
/*    For each text message */
      for (it=0; it < nt; it++) {
         nchar = s->ib[s->mp++];
         fprintf(stdout, "[%d] ", ipe);
         for (i=0; i < nchar; i++) {
            fprintf(stdout, "%c", (char)s->ib[s->mp++]);
         }
         fprintf(stdout, "\n");
      }
   }
   return(!(*error));
}
/* parc_sync_masmer:  master merges sibling text-status */
static void parc_sync_masmer(int *worst, bool *error) {
   int ipe,sibstatus;
#ifdef PARC_DUAL
   int i;
#endif
   DamPtr s;
/* Transfer local text-status to input buffer (ploy) */
   parc_tex2buf(&texdam, pestatus, ntex, &indam[0], error);
/* Zap current text-status */
   texdam.mp = 0;
   texdam.mn = 0;
   ntex = 0;
   pestatus = 0;
/* Prepare local master buffer for read */
   indam[0].mn = indam[0].mp;
   indam[0].mp = 0;
/* Find worst (largest positive) status value */
   *worst = 0;
   for (ipe=0; ipe < limpe; ipe++) {
      s = &indam[ipe];
#ifdef PARC_DUAL
/*    Number of words in buffer (ignore) */
      i = s->ib[s->mp++];
#endif
/*    Sibling status */
      sibstatus = s->ib[s->mp++];
      if (sibstatus > (*worst)) {
         *worst = sibstatus;
      }
   }
/* Print texts */
   parc_prtex(error);
}
/* parc_sync_mascas:  master adds broadcast worst status for siblings */
static void parc_sync_mascas(int worst, bool *error) {
   int ipe;
   DamPtr s;
/* Zero outdam for master */
      s = &outdam[0];
      s->mp = 0;
      s->mn = 0;
/* Broadcast worst status to siblings */
   for (ipe=1; ipe < limpe; ipe++) {
      s = &outdam[ipe];
/*    Rewind buffer for write */
      s->mp = 0;
      s->mn = 0;
/*    Add worst status for broadcast */
#ifdef PARC_DUAL
      if (dam_space(2, s, error)) {
         s->ib[s->mp++] = PARC_EMPTY;
         s->ib[s->mp++] = worst;
      }
#else
      if (dam_space(1, s, error)) {
         s->ib[s->mp++] = worst;
      }
#endif
   }
}
/* parc_sync_sibsend:  sibling sends in text-status */
static void parc_sync_sibsend(bool *error) {
   DamPtr s;
/* Prepare output buffer to receive text-status data */
   s = &outdam[0];
/* Transfer text to output buffer */
   parc_tex2buf(&texdam, pestatus, ntex, s, error);
/* Zap current text-status */
   texdam.mp = 0;
   texdam.mn = 0;
   ntex = 0;
   pestatus = 0;
/* Prepare buffer for write */
   s->mn = s->mp;
}
/* parc_sync_sibrec:  sibling grabs global status */
static void parc_sync_sibrec(int *ios, bool *error) {
   DamPtr s;
#ifdef PARC_DUAL
   int i;
#endif
   s = &indam[0];
#ifdef PARC_DUAL
/* Skip over buffer length */
   i = s->ib[s->mp++];
#endif
/* Grab and return status */
   *ios = s->ib[s->mp++];
}

/*
  Public performance functions
*/
/*
  Public Registration etc.
*/
/* parc_args:  acquire a copy of command line args */
void parc_args(int argc, char *argv[]) {
   int i;
   Argc = argc;
   Argv = argv;
}
/* parc_error:  report error condition */
void parc_error(char txt[], int *ios) {
   int i,nchar;
   bool error=FALSE;
   DamPtr s;
/* Adjust incoming status if not set by caller */
   if (*ios <= 0) {
      *ios = 1;
   }
/* If alive broadcast status */
   if (alive) {
      nchar = strlen(txt);
      while (nchar > 0) {
         if (txt[nchar-1] != ' ') {
            break;
         }
         nchar--;
      }
      s = &texdam;
      if (dam_space(nchar+1, s, &error)) {
         s->ib[s->mp++] = nchar;
         for (i=0; i < nchar; i++) {
            s->ib[s->mp++] = (int)txt[i];
         }
         ntex++;
      }
      if (*ios > pestatus) {
         pestatus = (*ios);
      }
   }
}
/* parc_warn:  warning message */
void parc_warn(char txt[]) {
   int i,nchar;
   bool error=FALSE;
   DamPtr s;
/* If alive broadcast status */
   if (alive) {
      nchar = strlen(txt);
      while (nchar > 0) {
         if (txt[nchar-1] != ' ') {
            break;
         }
         nchar--;
      }
      s = &texdam;
      if (dam_space(nchar+1, s, &error)) {
         s->ib[s->mp++] = nchar;
         for (i=0; i < nchar; i++) {
            s->ib[s->mp++] = (int)txt[i];
         }
         ntex++;
      } else {
         fprintf(stderr, "error:  enlarging text buffer\n");
      }
   }
}
/* parc_barrier:  simple barrier */
void parc_barrier(int *ios) {
   parc_mp_barrier(ios);
}
/* parc_sync:  simple syncrhonisation of text-status */
void parc_sync(int *ios) {
   int ipe,worst,nsend,sends[MAXPE],
     nrecv,recvs[MAXPE],i;
   bool error=FALSE;
   float t0,t1;
   if (*ios > pestatus) {
      pestatus = (*ios);
   }
   if (alive) {
/*    Master PE */
      if (mype==0) {
         parc_clock(&t0);
/*       Receive text-status from siblings */
#ifdef RECVSEQ
         nrecv = 0;
         for (i=1; i < limpe; i++) {
            recvs[nrecv++] = i;
         }
#else
         nrecv = limpe-1;
#endif
         parc_mp_recvs(nrecv, recvs, &error);
/*       Merge text-status and establish worst status */
         parc_sync_masmer(&worst, &error);
         parc_clock(&t1);
/*       Prefix output buffers for siblings with worst status */
         parc_clock(&t0);
         parc_sync_mascas(worst, &error);
/*       Send data to siblings */
         nsend = 0;
         for (ipe=1; ipe < limpe; ipe++) {
            sends[nsend++] = ipe;
         }
         parc_mp_sends(nsend, sends, &error);
         parc_clock(&t1);
         *ios = worst;
/*    Sibling PEs */
      } else {
         parc_sync_sibsend(&error);
/*       Send output to master */
         parc_clock(&t0);
         nsend = 0;
         sends[nsend++] = 0;
         parc_mp_sends(nsend, sends, &error);
         parc_clock(&t1);
/*       Receive worst status from master */
         parc_clock(&t0);
         nrecv = 1;
#ifdef RECVSEQ
         recvs[0] = 0;
#endif
         parc_mp_recvs(nrecv, recvs, &error);
         parc_clock(&t1);
         parc_sync_sibrec(ios, &error);
      }
   }
}
/* parc_rank:  return rank */
void parc_rank(int *mype_r, int *limpe_r, int *ios) {
/* Check incoming status */
   if (*ios != 0) {
      return;
   }
/* Check if dead (meaning attempt to rank after crash) */
   if (dead) {
      fprintf(stderr, "error:  attempt to rank when dead\n");
      *ios = 1;
/* Check if already alive, just return parameters */
   } else if (alive) {
      *mype_r = mype;
      *limpe_r = limpe;
   }
}
/* parc_register:  register for parallel execution */
void parc_register(int *mype_r, int *limpe_r, int *ios) {
   bool error=FALSE;
   int mpi_ios,i,ipe;
/* Check incoming status */
   if (*ios != 0) {
      return;
   }
/* Check if dead (meaning attempt to register after crash) */
   if (dead) {
      fprintf(stderr, "error:  attempt to register when dead\n");
      *ios = 1;
/* Check if already alive, just return parameters */
   } else if (alive) {
      *mype_r = mype;
      *limpe_r = limpe;
/* Check on type sizes */
   } else if (sizeof(float)!=sizeof(uword)) {
      fprintf(stderr, "error:  sizes of uword and float differ\n");
      *ios = 1;
/* Create text message data structure */
   } else if (!dam_init(&texdam, &error)) {
      fprintf(stderr, "error:  initialising text data\n");
      *ios = 1;
/* Register with implementation layer */
   } else if (!parc_mp_register(&mype, &limpe, &error)) {
      fprintf(stderr, "error:  failed to register with MPI\n");
      *ios = 1;
   } else {
/*    Initialise the various buffers */
      for (ipe=0; ipe < limpe; ipe++) {
         dam_init(&indam[ipe], &error);
         if (error) {
            fprintf(stderr, "error:  initializing indam[%d]\n",ipe);
            break;
         }
         dam_init(&outdam[ipe], &error);
         if (error) {
            fprintf(stderr, "error:  initializing outdam[%d]\n",ipe);
            break;
         }
      }
/*    if still no error then we are alive and kicking */
      if (error) {
         fprintf(stderr, "error:  initialising io buffers\n");
         *ios = 1;
      } else {
         alive = TRUE;
         *limpe_r = limpe;
         *mype_r = mype;
      }
   }
}
/* parc_exit:  tidy up and leave PARC */
void parc_exit() {
   if (alive) {
      fprintf(stdout, "parc_exit[%d]:  nsent=%d,nrecd=%d words\n",
        mype, nsent, nrecd);
      parc_mp_exit();
      alive = FALSE;
      dead = TRUE;
   } else {
      fprintf(stderr, "error:  PARC already exited\n");
   }
}
/*
   Broadcast functions
*/
/* parc_broadcast:  broadcast (non CHARACTER) data */
static void parc_broadcast(uword buf[], int *type, int *ndim, int dims[],
  int *ios) {
   int ipe,worst,nsend,sends[MAXPE],nrecv,recvs[MAXPE],i,ni;
   bool error=FALSE;
   DamPtr s;
   float t0,t1;
   int nw,idim;
#ifdef DEBUG_BROADCAST
   fprintf(stdout, "parc_broadcast[%d]:  starts\n", mype);
#endif
/* Check number dimensions etc. */
   if ((*ndim) < 0) {
      fprintf(stderr, "error[%d]:  ndim < 1 (parc_broadcast)\n", mype);
      *ios = 1;
/* Calculate number of words to be broadcast */
   } else {
      nw = 1;
      for (idim=0; idim < (*ndim); idim++) {
         nw = nw*dims[idim];
         if (dims[idim] < 1) {
            fprintf(stderr, "error[%d]:  dims[%d]=%d (parc_broadcast)\n",
              mype, idim, dims[idim]);
            *ios = 1;
            break;
         }
      }
   }
/* Blank out data if we have an error */
   if (*ios != 0) {
      nw = 0;
   }
/* Grab a copy of the incoming status */
   if (*ios > pestatus) {
      pestatus = (*ios);
   }
   if (alive) {
/*    Master PE */
      if (mype==0) {
         parc_clock(&t0);
/*       Receive text-status from siblings */
#ifdef RECVSEQ
         nrecv = 0;
         for (i=1; i < limpe; i++) {
            recvs[nrecv++] = i;
         }
#else
         nrecv = limpe-1;
#endif
         parc_mp_recvs(nrecv, recvs, &error);
/*       Merge text-status and establish worst status */
         parc_sync_masmer(&worst, &error);
         parc_clock(&t1);
/*       Prefix output buffers for siblings with worst status */
         parc_clock(&t0);
         parc_sync_mascas(worst, &error);
/*       Append broadcast data */
         for (ipe=1; ipe < limpe; ipe++) {
            s = &outdam[ipe];
            if (dam_space(nw+1, s, &error)) {
               s->ib[s->mp++] = nw;
               for (i=0; i < nw; i++) {
                  s->mb[s->mp++] = buf[i];
               }
            }
         }
/*       Send data to siblings */
         nsend = 0;
         for (ipe=1; ipe < limpe; ipe++) {
            sends[nsend++] = ipe;
         }
         parc_mp_sends(nsend, sends, &error);
         parc_clock(&t1);
         *ios = worst;
/*    Sibling PEs */
      } else {
         parc_sync_sibsend(&error);
/*       Send output to master */
         parc_clock(&t0);
         nsend = 0;
         sends[nsend++] = 0;
         parc_mp_sends(nsend, sends, &error);
         parc_clock(&t1);
/*       Receive input from master */
         parc_clock(&t0);
         nrecv = 1;
#ifdef RECVSEQ
         recvs[0] = 0;
#endif
         parc_mp_recvs(nrecv, recvs, &error);
         parc_clock(&t1);
/*      Acquire merged status from master input */
         parc_sync_sibrec(ios, &error);
/*      Acquire broadcast data from master input */
         s = &indam[0];
         ni = s->ib[s->mp++];
         for (i=0; i < ni; i++) {
            buf[i] = s->mb[s->mp++];
         }
      }
   }
}
/* parc_broadcastr:  broadcast (REAL) data */
void parc_broadcastr(float buf[], int *ndim, int dims[],
  int *ios) {
   int type=PARC_FLOAT_TYPE;
   parc_broadcast((uword*)buf, &type, ndim, dims, ios);
}
/* parc_broadcast0c:  broadcast (scalar CHARACTER) data */
void parc_broadcast0c(char buf[],  int *ios) {
   int ipe,worst,nsend,sends[MAXPE],nrecv,recvs[MAXPE],i,ni;
   bool error=FALSE;
   DamPtr s;
   float t0,t1;
   int nw;
/* Calculate number of words to be broadcast */
   nw = strlen(buf) + 1;
/* Blank data if we have an error */
   if (*ios != 0) {
      nw = 0;
   }
/* Grab a copy of the incoming status */
   if (*ios > pestatus) {
      pestatus = (*ios);
   }
   if (alive) {
/*    Master PE */
      if (mype==0) {
/*       Receive text-status from siblings */
#ifdef RECVSEQ
         nrecv = 0;
         for (i=1; i < limpe; i++) {
            recvs[nrecv++] = i;
         }
#else
         nrecv = limpe-1;
#endif
         parc_mp_recvs(nrecv, recvs, &error);
/*       Merge text-status and establish worst status */
         parc_sync_masmer(&worst, &error);
/*       Prefix output buffers for siblings with worst status */
         parc_clock(&t0);
         parc_sync_mascas(worst, &error);
/*       Append broadcast data */
         for (ipe=0; ipe < limpe; ipe++) {
            s = &outdam[ipe];
            if (dam_space(nw+1, s, &error)) {
               s->ib[s->mp++] = nw;
               for (i=0; i < nw; i++) {
                  s->mb[s->mp++] = (uword)buf[i];
               }
            }
         }
/*       Send data to siblings */
         nsend = 0;
         for (ipe=1; ipe < limpe; ipe++) {
            sends[nsend++] = ipe;
         }
         parc_mp_sends(nsend, sends, &error);
         *ios = worst;
/*    Sibling PEs */
      } else {
         parc_sync_sibsend(&error);
/*       Send output to master */
         nsend = 0;
         sends[nsend++] = 0;
         parc_mp_sends(nsend, sends, &error);
/*       Receive input from master */
         nrecv = 1;
#ifdef RECVSEQ
         recvs[0] = 0;
#endif
         parc_mp_recvs(nrecv, recvs, &error);
         parc_clock(&t1);
/*      Acquire merged status from master input */
         parc_sync_sibrec(ios, &error);
/*      Acquire broadcast data from master input */
         s = &indam[0];
         ni = s->ib[s->mp++];
         for (i=0; i < ni; i++) {
            buf[i] = (char)s->mb[s->mp++];
         }
      }
   }
}
/* parc_xdata:  data exchange (non CHARACTER) */
static void parc_xdata(int *type, int *npe, 
  int outs[], uword outbuf[], 
  int ins[], uword inbuf[],
  int *ios) {
#ifdef PARC_SMP
  /* SMP use old Knuth algorithm */
   int ipe,worst,nsend,sends[MAXPE],nrecv,recvs[MAXPE],i,ni,nwpe,k;
   int gouts[MAXPE][MAXPE],nsendx,sendxs[MAXPE],nrecvx,recvxs[MAXPE],
     mo,mi,jpe,irecv;
   bool error=FALSE;
   float *fbuf;
   DamPtr s;
   float t0,t1;
   int nw,idim;
/* Check for bad args and update ios if necessary */
   if ((*npe) > limpe) {
      fprintf(stderr, "error[%d]:  npe > limpe (parc_xdata)\n", mype);
      *ios = 1;
   }
/* Grab a copy of the incoming status */
   if (*ios > pestatus) {
      pestatus = (*ios);
   }
   if (alive) {
/*    Master PE */
      if (mype==0) {
         parc_clock(&t0);
/*       Receive text-status from siblings */
#ifdef RECVSEQ
         nrecv = 0;
         for (i=1; i < limpe; i++) {
            recvs[nrecv++] = i;
         }
#else
         nrecv = limpe-1;
#endif
         parc_mp_recvs(nrecv, recvs, &error);
/*       Merge text-status and establish worst status */
         parc_sync_masmer(&worst, &error);
/*       Zero global counts array up to limpe */
         for (ipe=0; ipe < limpe; ipe++) {
            for (jpe=0; jpe < limpe; jpe++) {
               gouts[ipe][jpe] = 0;
            }
         }
/*       Copy outs data from local PE */
         for (jpe=0; jpe < (*npe); jpe++) {
            gouts[0][jpe] = outs[jpe];
         }
/*       Acquire outs data from sibling input */
         for (ipe=1; ipe < (*npe); ipe++) {
            s = &indam[ipe];
            ni = s->ib[s->mp++];
            for (jpe=0; jpe < (*npe); jpe++) {
               gouts[ipe][jpe] = s->ib[s->mp++];
            }
         }
/*       Count senders to master */
         nrecvx = 0;
         for (ipe=1; ipe < (*npe); ipe++) {
            if (gouts[ipe][0] > 0) {
               recvxs[nrecvx++] = ipe;
            }
         }
/*       Prefix output buffers for siblings with worst status */
         parc_sync_mascas(worst, &error);
/*       Append broadcast data */
         for (ipe=0; ipe < limpe; ipe++) {
            s = &outdam[ipe];
            if (dam_space((*npe)+1, s, &error)) {
               if (ipe < (*npe)) {
                  s->ib[s->mp++] = (*npe);
                  for (jpe=0; jpe < (*npe); jpe++) {
                     s->mb[s->mp++] = gouts[jpe][ipe];
                  }
               } else {
                  s->ib[s->mp++] = 0;
               }
            }
         }
/*       Send data to siblings */
         nsend = 0;
         for (ipe=1; ipe < limpe; ipe++) {
            sends[nsend++] = ipe;
         }
         parc_mp_sends(nsend, sends, &error);
         *ios = worst;
/*    Sibling PEs */
      } else {
         parc_sync_sibsend(&error);
/*       Append out data */
         s = &outdam[0];
         if (dam_space((*npe)+1, s, &error)) {
            if (mype < (*npe)) {
               s->ib[s->mp++] = (*npe);
               for (ipe=0; ipe < (*npe); ipe++) {
                  s->ib[s->mp++] = outs[ipe];
               }
            } else {
               s->ib[s->mp++] = 0;
            }
         }
/*       Send output to master */
         nsend = 0;
         sends[nsend++] = 0;
         parc_mp_sends(nsend, sends, &error);
         parc_clock(&t1);
/*       Receive input from master */
         parc_clock(&t0);
         nrecv = 1;
#ifdef RECVSEQ
         recvs[0] = 0;
#endif
         parc_mp_recvs(nrecv, recvs, &error);
/*      Acquire merged status from master input */
         parc_sync_sibrec(ios, &error);
/*      Acquire broadcast data from master input */
         s = &indam[0];
/*       Number of PEs in activity */
         ni = s->ib[s->mp++];
         nrecvx = 0;
         for (ipe=0; ipe < ni; ipe++) {
            ins[ipe] = s->ib[s->mp++];
/*          Increment number of expected senders */
            if (ins[ipe] > 0) {
               recvxs[nrecvx++] = ipe;
            }
         }
      }
/*    Fill output buffers */
      nsendx = 0;
      mo = 0;
      for (ipe=0; ipe < (*npe); ipe++) {
         if (ipe != mype) {
            if (outs[ipe] > 0) {
               s = &outdam[ipe];
               s->mp = 0;
               if (dam_space(outs[ipe]+2, s, &error)) {
                  s->ib[s->mp++] = PARC_EMPTY;
                  s->ib[s->mp++] = outs[ipe];
                  for (i=0; i < outs[ipe]; i++) {
                     s->mb[s->mp++] = outbuf[mo++];
                  }
               }
               sendxs[nsendx++] = ipe;
            }
         }
      }
/*    Now do a global send-receive using the primitive */
      parc_sendrecv(nsendx, sendxs, nrecvx, recvxs, &error);
/*    Zero argument array */
      for (ipe=0; ipe < limpe; ipe++) {
         ins[ipe] = 0;
      }
/*    Acquire incoming data */
      mi = 0;
      for (irecv=0; irecv < nrecvx; irecv++) {
         ipe = recvxs[irecv];
         s = &indam[ipe];
         s->mp = 0;
         ni = s->ib[s->mp++];
         ins[ipe] = s->ib[s->mp++];
         for (i=0; i < ins[ipe]; i++) {
            inbuf[mi++] = s->mb[s->mp++];
         }
      }
   }
}
#else
   /* MPI us MPI_Alltoallv */
  int scounts[MAXPE], rcounts[MAXPE], rdispls[MAXPE], sdispls[MAXPE], 
    nproc, soff, roff;
  /* Prepare to send no of outs to each processor */
  for (nproc = 0; nproc < *npe; nproc++){
      scounts[nproc] = 1;
      rcounts[nproc] = *npe;
      rdispls[nproc] = nproc;
      sdispls[nproc] = nproc;
  }
  /* Send outs */
  MPI_Alltoallv(outs, scounts, sdispls, MPI_INT,
    ins, rcounts, rdispls, MPI_INT, MPI_COMM_WORLD);
  /* Prepare for main message */
  roff = 0;
  soff = 0;
  for (nproc = 0; nproc < *npe; nproc++){
    scounts[nproc] = outs[nproc];
    rcounts[nproc] = ins[nproc];
    rdispls[nproc] = roff;
    sdispls[nproc] = soff;
    roff += ins[nproc];
    soff += outs[nproc];
  }
  /* Send main message */
  MPI_Alltoallv(outbuf, scounts, sdispls, mpi_uword,
    inbuf, rcounts, rdispls, mpi_uword, MPI_COMM_WORLD);
}
#endif
/* parc_xdatar:  data exchange (REAL) */
void parc_xdatar(int *npe, 
  int outs[], float outbuf[], 
  int ins[], float inbuf[],
  int *ios) {
   int type=PARC_FLOAT_TYPE;
   int i,k,ipe;
   parc_xdata(&type, npe, outs, (uword *)outbuf,
     ins, (uword *)inbuf, ios);
}
/* parc_size: print sizes */
void parc_size()
{
  int npe;
  printf("mype=%d\n", mype);
  for (npe = 0; npe < limpe; npe++)
    printf("npe=%d, indam->mm=%d, outdam->mm=%d\n",npe, 
	   indam[npe].mm, outdam[npe].mm);  
}
/* parc_reset: reset buffers */
void parc_reset()
{
  int npe;
  DamPtr s;
  uword *newb;
  /* Reset buffers for each processor */
  for (npe = 0; npe < limpe; npe++){
    /* Input buffers */
    s = &indam[npe];
    free(s->mb);
    if ((newb=(uword*)malloc((size_t)
      (DAM_MAX_BLOCK*sizeof(uword))))==(uword*)NULL){
      printf("error:  allocating memory block (parc_reset)\n");
    }
    else{
      s->mb = newb;
      s->ib = (lint*)s->mb;
      s->fb = (float*)s->mb;
      s->mm = DAM_MAX_BLOCK;
      s->mp = 0;
      s->mn = 0;
    }
    /* Output buffers */
    s = &outdam[npe];
    free(s->mb);
    if ((newb=(uword*)malloc((size_t)
      (DAM_MAX_BLOCK*sizeof(uword))))==(uword*)NULL){
      printf("error:  allocating memory block (parc_reset)\n");
    }
    else{
      s->mb = newb;
      s->ib = (lint*)s->mb;
      s->fb = (float*)s->mb;
      s->mm = DAM_MAX_BLOCK;
      s->mp = 0;
      s->mn = 0;
    }
  }
}
