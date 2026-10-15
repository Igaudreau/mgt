/* smp.c:  SMP message passing library (for debug)
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include "const.h"
#include "smp.h"

/* Implementation layer
*/
/* Constants :
*/
#define MAXPE 256
#define ICONIC
#define MAXARG 512
#define SZARG 513
#define SZPATH 513

/* Data :
*/
static int mype=0,npe=1,smpid=0;
static bool alive=FALSE,debug=FALSE,valgrind=FALSE;
static int fdi[MAXPE],fdo[MAXPE];
static pid_t pids[MAXPE];
static uword tv_sec0=0;
static bool notv_sec0=TRUE;

/* parc_clock:  wall clock time (from system) */
void smp_clock(float *t) {
   struct timeval f;
   gettimeofday(&f, (struct timezone *)NULL);
/* Capture suitable offset to avoid loss of floating precision */
   if (notv_sec0) {
      tv_sec0 = (uword)(f.tv_sec);
      notv_sec0 = FALSE;
   }
   *t = ((float)(f.tv_sec-tv_sec0)) + ((float)f.tv_usec)*0.000001;
}
/* smp_spawn:  spawn sibling instances */
static bool smp_spawn(int argc, char *argv[]) {
   bool error=FALSE;
   char args[MAXARG][SZARG];
   char *newargv[MAXARG];
   char path[SZPATH];
   int ipe,i,newargc,nchar,offset,dum;
   pid_t pid;
/* Check that we are called by master */
   if (mype != 0) {
      fprintf(stderr, "error[%d]:  master must spawn (smp_spawn)\n",
        mype);
      error = TRUE;
      return(!error);
   }
/* For each sibling */
   for (ipe=1; ipe < npe; ipe++) {
/*    Initialise outgoing argument counter */
      newargc = 0;
/*    Path is always an xterm */
      nchar = strlen("xterm");
      if (nchar+1 > SZPATH) {
         fprintf(stderr, "error[%d]:  path overflow (smp_spawn)\n",
           mype);
         error = TRUE;
         break;
      }
      sprintf(path, "xterm");
      sprintf(args[newargc++], "xterm");
/*    Some other flags for xterm */
#ifdef ICONIC
      sprintf(args[newargc++], "-iconic");
#endif
      sprintf(args[newargc++], "-geometry");
      sprintf(args[newargc++], "72x24");
      sprintf(args[newargc++], "-fg");
      sprintf(args[newargc++], "black");
      sprintf(args[newargc++], "-bg");
      sprintf(args[newargc++], "yellow");
      sprintf(args[newargc++], "-fn");
      sprintf(args[newargc++], "7x13bold");
      sprintf(args[newargc++], "-title");
      offset = 0;
      for (i=strlen(argv[0])-1; i >= 0; i--) {
         if (argv[0][i]=='/') {
            offset = i+1;
            break;
         }
      }
      sprintf(args[newargc++], "%s[%d]", &(argv[0][offset]), ipe);
/*    -e flag */
      sprintf(args[newargc++], "-e");
/*    Run GDB as the xterm application */
      if (debug) {
         sprintf(args[newargc++], "gdb");
         sprintf(args[newargc++], "-quiet");
         sprintf(args[newargc++], "--command=.smp/%d/gdb.run.%d", smpid, ipe);
/*    Run with valgrind */
      } else if (valgrind) {
         sprintf(args[newargc++], "valgrind");
      }
/*    Pass on program name argument (only) to execvp  */
      nchar = strlen(argv[0]);
      if (nchar+1 > SZARG) {
         fprintf(stderr, "error[%d]  arg too long (smp_spawn)\n",
           mype);
         error = TRUE;
         break;
      }
      sprintf(args[newargc++], "%s", argv[0]);
/*    Add remaining program arguments (avoid -mype=xx) */
      if (!debug) {
         for (i=1; i < argc; i++) {
            nchar = strlen(argv[i]);
            if (sscanf(argv[i], "-mype=%d\n", &dum)==1) {
               sprintf(args[newargc++], "-mype=%d", ipe);
            } else if (nchar+1 > SZARG) {
               fprintf(stderr, "error[%d]  arg too long (smp_spawn)\n",
                 mype);
               error = TRUE;
               break;
            } else {
               sprintf(args[newargc++], "%s", argv[i]);
            }
         }
      }
/*    Construct newarg pointers */
      for (i=0; i < newargc; i++) {
         newargv[i] = &(args[i][0]);
      }
/*    Add a null pointer */
      newargv[newargc] = (char*)NULL;
/*    Fork process */
      if ((pid=fork()) < 0) {
         fprintf(stderr, "error[%d]:  forking sibling\n", mype);
         error = TRUE;
         break;
/*    We are in the sibling process */
      } else if (pid == 0) {
         execvp(path, newargv);
         exit(EXIT_FAILURE);
/*    We are still in the master processes */
      } else {
         pids[ipe] = pid;
      }
   }
   return(!error);
}
/* smp_mkfifos:  create matrix of FIFOs */
static bool smp_mkfifos() {
   mode_t mode_dir,mode_fifo;
   char name[SZNAM],file[SZPATH];
   int i,j;
   bool error=FALSE;
/* Return without action if not master */
   if (mype != 0) {
      return(!error);
   }
/* Permissions mode */
   mode_fifo = S_IRUSR | S_IWUSR;
   mode_dir = S_IRUSR | S_IWUSR | S_IXUSR;
/* Attempt to creathe .smp subdirectory (OK if exists) */
   mkdir(".smp", mode_dir);
   sprintf(file, ".smp/%d", smpid);
   mkdir(file, mode_dir);
/* Full connectivity matrix (send from "i" to "j") */
   for (i=0; i < npe; i++) {
      for (j=0; j < npe; j++) {
/*       Construct simple name */
         sprintf(name, "%s/.fifo_%d_%d", file, i, j);
/*       Attempt to make the FIFO */
         if (mkfifo(name,mode_fifo) < 0) {
/*       Attempt to change permissions (anyway) */
         } else if (chmod(name, mode_fifo) != 0) {
            fprintf(stderr, "error:  failed to chmod fifo=%s\n",
              name);
            break;
            error = TRUE;
         }
      }
      if (error) {
         break;
      }
   }
   return(!error);
}
/* smp_openfifos:  open FIFOs for siblings */
static bool smp_openfifos() {
   char nameo[SZPATH],namei[SZPATH];
   int i;
   bool error=FALSE;
/* Let us try using O_RDWR to avoid blocking */
   for (i=0; i < npe; i++) {
/*    Construct FIFO names */
      sprintf(nameo,".smp/%d/.fifo_%d_%d", smpid, mype, i);
      sprintf(namei,".smp/%d/.fifo_%d_%d", smpid, i, mype);
/*    Avoid talking to myself */
      if (mype==i) {
         fdi[i] = -1;
         fdo[i] = -1;
/*    Open FIFO for writing to sibling */
      } else if ((fdo[i]=open(nameo,O_RDWR)) < 0) {
         fprintf(stderr, "error[%d]:  opening sibling fifo=%s (out)\n",
           mype, nameo);
         perror("reason: ");
         error = TRUE;
         break;
/*    Open FIFO for reading from sibling */
      } else if ((fdi[i]=open(namei,O_RDWR)) < 0) {
         fprintf(stderr, "error[%d]:  opening sibling fifo=%s (in)\n",
           mype, namei);
         perror("reason: ");
         error = TRUE;
         break;
      }
   }
   return(!error);
}
/* smp_send:  send message to sibling */
bool smp_send(int ipe, int nw, uword buf[]) {
   bool error=FALSE;
   int offset,nsent;
   uchar *cbuf;
   unsigned int nbyte;
/* Check that we are alive */
   if (!alive) {
      fprintf(stderr, "error:  called when unregistered (smp_send)\n");
      error = TRUE;
      return(!error);
   }
/* Grab character pointer for buf */
   cbuf = (uchar*)buf;
/* Offset into buffer (words) */
   offset = 0;
/* Number of bytes remaining to be sent */
   nbyte = nw*sizeof(uword);
/* Repeatedly try until we succeed */
   while (offset < nw*sizeof(uword)) {
/*    Attempt to write what remains of the buffer */
      nsent=write(fdo[ipe],(void*)(&(cbuf[offset])),nbyte); 
/*    Error detected */
      if (nsent == (-1)) {
/*       Only react if not interrupted */
         if (errno != EINTR) {
            fprintf(stderr, "error[%d]:  write error on ipe=%d\n",
              mype, ipe);
            perror("reason: ");
            error = TRUE;
            break;
         } else {
            fprintf(stderr, "warning[%d]:  write interrupt on ipe=%d\n",
              mype, ipe);
         }
/*    Otherwise update offset and number of bytes remaining */
      } else {
         nbyte = nbyte - nsent;
         offset = offset + nsent;
      }
   }
   return(!error);
}
/* smp_recv:  receive message from sibling */
bool smp_recv(int ipe, int nw, uword buf[]) {
   bool error=FALSE;
   int offset,nrecv;
   uchar *cbuf;
   unsigned int nbyte;
/* Check that we are alive */
   if (!alive) {
      fprintf(stderr, "error:  called when unregistered (smp_recv)\n");
      error = TRUE;
      return(!error);
   }
/* uchar pointer to buffer */
   cbuf = (uchar*)buf;
/* Offset into buffer (words) */
   offset = 0;
/* Number of bytes remaining to be received */
   nbyte = nw*sizeof(uword);
/* Repeatedly try until we succeed */
   while (offset < nw*sizeof(uword)) {
/*    Attempt to read what remains of the buffer */
      nrecv=read(fdi[ipe],(void*)(&(cbuf[offset])),nbyte); 
/*    Error detected */
      if (nrecv == (-1)) {
/*       Only react if not interrupted */
         if (errno != EINTR) {
            fprintf(stderr, "error[%d]:  read error on ipe=%d\n",
              mype, ipe);
            perror("reason: ");
            error = TRUE;
            break;
         } else {
            fprintf(stderr, "warning[%d]:  read interrupt on ipe=%d\n",
              mype, ipe);
         }
/*    Otherwise update offset and number of bytes remaining */
      } else {
         nbyte = nbyte - nrecv;
         offset = offset + nrecv;
      }
   }
   return(!error);
}
/* smp_exit:  exit smp, clean up mess */
void smp_exit() {
   int ipe;
/* If alive, we try to kill off the siblings */
   if (alive) {
/*    master tries to kill off siblings (the xterms) */
      if (mype==0) {
         for (ipe=1; ipe < npe; ipe++) {
            if (kill(pids[ipe], SIGTERM) != 0) {
               fprintf(stderr, "error[%d]:  killing ipe=%d\n",
                 mype, ipe);
            }
         }
      }
   }
}
/* smp_abort:  exit smp, clean up mess */
void smp_abort() {
   int ipe;
/* If alive, we try to kill off the siblings */
   if (alive) {
/*    master tries to kill off siblings (the xterms) */
      if (mype==0) {
         for (ipe=1; ipe < npe; ipe++) {
            if (kill(pids[ipe], SIGTERM) != 0) {
               fprintf(stderr, "error[%d]:  killing ipe=%d\n",
                 mype, ipe);
            }
         }
      }
   }
}
/* smp_barrier:  create a pseudo-synchronisation barrier */
void smp_barrier() {
   uword buf[1];
   int ipe;
   if (mype==0) {
      for (ipe=1; ipe < npe; ipe++) {
         smp_recv(ipe, 1, buf);
      }
      for (ipe=1; ipe < npe; ipe++) {
         smp_send(ipe, 1, buf);
      }
   } else {
      smp_send(0, 1, buf);
      smp_recv(0, 1, buf);
   }
}
/* smp_register:  register and get rank */
bool smp_register(int argc, char *argv[], int* npe0, int* mype0) {
   bool error=FALSE;
   int i;
/* Check that we are dead */
   if (alive) {
      fprintf(stderr, "error:  called when registered (smp_register)\n");
      error = TRUE;
      return(!error);
   }
/* Missing -mype argument implies master */
   mype = 0;
/* Scan through supplied command-line arguments for -mype, -npe etc. */
   for (i=0; i < argc; i++) {
      if (sscanf(argv[i], "-mype=%d", &mype) == 1) {
         *mype0 = mype;
      } else if (sscanf(argv[i], "-npe=%d", &npe) == 1) {
         *npe0 = npe;
      } else if (sscanf(argv[i], "-smpid=%d", &smpid) == 1) {
      } else if (strcmp(argv[i], "-debug") == 0) {
         debug = TRUE;
         valgrind = FALSE;
      } else if (strcmp(argv[i], "-valgrind") == 0) {
         debug = FALSE;
         valgrind = TRUE;
      }
   }
/* Check on number of PEs */
   if ((npe < 1)||(npe>MAXPE)) {
      fprintf(stderr, "error[%d]:  npe=%d bad (smp_register)\n", 
        mype, npe);
      error = TRUE;
/* Check mype for sanity */
   } else if ((mype < 0)||(mype > (npe-1))) {
      fprintf(stderr, "error[%d]:  mype=%d bad for npe=%d (smp_register)\n",
        mype, mype, npe);
      error = TRUE;
/* Create FIFO matrix (returns without error if not master) */
   } else if (!smp_mkfifos()) {
      fprintf(stderr, "error[%d]:  creating FIFOs (smp_register)\n",
       mype);
      error = TRUE;
/* Open FIFO channels */
   } else if (!smp_openfifos()) {
      fprintf(stderr, "error[%d]:  creating FIFOs (smp_register)\n",
       mype);
      error = TRUE;
/* If we are master spawn siblings so that they can register */
   } else if (mype > 0) {
   } else if (!smp_spawn(argc, argv)) {
      fprintf(stderr, "error[%d]:  spawning siblings\n", mype);
      error = TRUE;
   }
/* are we alive? */
   alive = !error;
   return(!error);
}
