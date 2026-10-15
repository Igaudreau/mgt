#undef MAIN
#include <fstream>	// added by BM for to allow nohup execution
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <new>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <signal.h>
#include "iop.h"
extern void broadcast(int *value);
extern void broadcast(float *value);
extern void broadcast(bool *value);
extern void broadcast_db(double *value);
extern "C"
{
void  parc_rank(int *mype, int *limpe, int *ios);
void  parc_broadcast0c(char temp[], int *ios);
}
void activ(char *prompt);
void deactiv();
bool get_toks();
bool get_cmd(char key[]);
bool get_n(char name[]);
bool get_i(int &value);
bool get_r(float &value);
bool get_db(double &value);
bool gen_cmd(char cmd[]);
void sys_cmd(char cmd[]);
bool open_dat(const char tail[], char name[], bool techo);
bool open_log(char name[]);
void report();

//file activ.c
/*
+ activ - activate command interpreter.
     Description

        This activates the command interpreter: sets the local prompt;
        sets the start of token buffer. key must be set to NULL in the
        calling routine.
*/
void activ(char *prompt)
{
  int i, ios, limpe;
  /* Activate */
  ftok = NULL;
  readf = false;
  logon = false;
  idat = -1;
  for (i = 0; i < NDAT; i++)
    datin[i] = NULL;				// BM fix
    //datin[idat] = NULL;			// bug here
  strcpy(lprompt, prompt);
  /* Get process rank */
  ios = 0;
  parc_rank(&tmype, &limpe, &ios);
  return;
}
//end activ.c

//file deactiv.c
/*
+ deactiv - deactivate command interpreter.
     Description:
        This deactivates  the command  interpreter by cleaning  up the
        remainder of the  buffer if any. key has to  be deleted in the
        calling routine.
*/
void deactiv()
{
  struct token *cptr;
  /* Deactivate */
  if (ftok != NULL){
    cptr = ftok;
    while (cptr != NULL){
      ftok = cptr;
      cptr = cptr->nxt;
      delete ftok;
    }
  }
  if (logon)
    fclose(logout);
  return;
}
//end activ.c

//file get_toks.c
/*
+ get_toks - get tokens.
     Description:
        This  gets  a  series  of  tokens from  an  input  line  using
        readline. The tokens are put  in a buffer consisting of a list
        list of structures  of type token. The pointer  ftok points to
        the first of  these. This should only be  called by the master
        processor.
*/
bool get_toks()
{
  size_t lenget;
  int i, lentemp, lenquote;
  //int ios, limpe, mype;
  struct token *cptr;
  static char temp[SZNAM], *line, *tok, *qptr;
  /* Intialize */
  line = (char *)NULL;
  tok = (char *)NULL;
  qptr = (char *)NULL;
  /* Get line */
  if (readf){
    /* From input file */
    if(getline(&line, &lenget, datin[idat]) != -1)
      /* Remove newline */
      line[strlen(line)-1] = '\0';
    else{
      /* End of input file */
      line = NULL;
      fclose(datin[idat]);
      ftok = oftok[idat];
      echo = false;
      idat--;
      readf = (idat >= 0);
    }
  }
  else{
    /* check if auto.dat exists - added by BM for nohup execution */
    ifstream ifile;
    ifile.open("auto.dat");
    if(ifile){
      cout << "loading commands from auto.dat..." << endl;
      line = new char [SZNAM];
      strcpy(line, "$read auto q");
      ifile.close();
    }
    else {
      /* From terminal */
      line = readline(lprompt);
    }
    if (line != NULL)
      add_history (line);
  }
  /* Return if no input */
  if (line == NULL)
    return(false);
  /* Return if line too long */
  if (strlen(line) > SZNAM - 1){
    free(line);
    return(false);
  }
  /* Check echo */
  if (echo)
    cout << line << '\n';
  /* Remove comment (if any) */
  if (line[0] == '!'){
    /* Whole line is a comment -- ignore line */
    free(line);
    return(get_toks());
  }
  else{
    /* Ignore line after comment */
    if (strstr(line, "!") != NULL){
      tok = strtok(line, "!");
      strcpy(temp, tok);
      strcat(temp, " ");
    }
    else
      strcpy(temp,line);
  }
  if (line != NULL)
    free(line);
  /* Check for log file */
  if (logon && (tmype == 0))
    fprintf(logout, "%s \n", line);
  /* Replace tabs with spaces */
  if (strchr(temp, '\t') != NULL){
    lentemp = strlen(temp);
    for (i = 0; i < lentemp; i++){
      if (temp[i] == '\t')
	temp[i] = ' ';
    }
  }
  /* Check for quoted strings - mark spaces with '|' */
  if (strchr(temp, '|') != NULL){
    cout << "Warning: | is used as a delimiter --"; 
    cout << " it will be replaced by space \n";
    return(false);
  }
  qptr = strchr(temp, '\"');
  if (qptr != NULL){
    qptr[0] = ' ';
    lentemp = strlen(temp);
    lenquote = strlen(qptr);
    for (i = lentemp - lenquote; i < lentemp; i++){
      if (temp[i] == '\"'){
	temp[i] = ' ';
	break;
      }
      else if (temp[i] == ' ')
	temp[i] = '|';
    }
  }
  /* Tokenize line */
  tok = strtok(temp," ");
  if (tok != NULL){
    ftok = new token;
    ftok->nxt = NULL;
    strcpy(ftok->key, tok);
  }
  else
    return(false);
  /* Put tokens in buffer */
  cptr = ftok;
  do {
    tok = strtok(NULL," ");
    if (tok !=NULL){
      cptr->nxt = new token;
      cptr = cptr->nxt;
      cptr->nxt = NULL;
      strcpy(cptr->key, tok);
    }
  } while (tok);
  /* Check for quoted strings */
  if (qptr != NULL){
    cptr = ftok;
    while (cptr != NULL){
      tok = strchr(cptr->key, '|');
      if (tok != NULL){
	/* Clean quoted string */
	//for (i = 0; i < strlen(tok); i++)
	for (i = 0; i < static_cast<int>(strlen(tok)); i++)	// BM
	  if (tok[i] == '|')
	    tok[i] = ' ';
      }
      cptr = cptr->nxt;
    }
  }
  return(true);
}
//end get_toks.c

//file get_cmd.c
/*
+ get_cmd - get command.
     Description:
        The master  gets a  command string from  the buffer set  up by
        get_toks  and removes  it  from the  buffer.   The command  is
        returned in  key and the  return value is true  if successful,
        false otherwise.   The command is converted to  upper case. If
        successful it is then broadcast to the other processors.
*/
bool get_cmd(char cmd[]) 
{
  int i, ios;
  bool fail;
  struct token *cptr;
  if (tmype == 0){
    fail = false;
    /* Get new line if buffer is empty */
    if (ftok == NULL)
	fail = !get_toks();
    if (!fail){
      /* Extract command from buffer */
      strcpy(cmd, ftok->key);
      //for (i = 0; i < strlen(ftok->key); i++)
      for (i = 0; i < static_cast<int>(strlen(ftok->key)); i++)	// BM
	cmd[i] = toupper(cmd[i]);
      /* Remove start of buffer */
      cptr = ftok;
      ftok = cptr->nxt;
      delete cptr;
    }
  }
  /* Broadcast fail */
  broadcast(&fail);
  if (fail)
    return(false);
  else{
    /* Broadcast command */
    ios = 0;
    parc_broadcast0c(cmd, &ios);
    return(true);
  }
}
//end get_cmd.c

//file get_n.c
/*
+ get_n - get name.
     Description:
        The master  gets a  command string from  the buffer set  up by
        get_toks  and removes  it  from the  buffer.   The command  is
        returned in  key and the  return value is true  if successful,
        false  otherwise. If successful  it is  then broadcast  to the
        other processors.
*/
bool get_n(char name[]) 
{
  int ios;
  bool fail;
  struct token *cptr;
  if (tmype == 0){
    fail = false;
    /* Get new line if buffer is empty */
    if (ftok == NULL)
	fail = !get_toks();
    if (!fail){
      /* Extract command from buffer */
      strcpy(name, ftok->key);
      /* Remove start of buffer */
      cptr = ftok;
      ftok = cptr->nxt;
      delete cptr;
    }
  }
  /* Broadcast fail */
  broadcast(&fail);
  if (fail)
    return(false);
  else{
    /* Broadcast name */
    ios = 0;
    parc_broadcast0c(name, &ios);
    return(true);
  }
}
//end get_n.c

//file get_i.c
/*
+ get_i - get integer value.
     Description:
        The master gets an integer value from the buffer of tokens set
        up by  get_toks. The  token is removed  from the  buffer.  The
        value  of the  integer is  returned in  value and  the routine
        returns  true if successful,  false otherwise.  If successful,
        the value is broadcast to the other processors.
*/
bool get_i(int &value)
{
  struct token *cptr;
  int i, ind ;
  bool fail;
  if (tmype == 0){
    fail = false;
    /* Get new line if buffer is empty */
    if (ftok == NULL)
      fail = !get_toks();
    if (!fail){
      /* Extract value from buffer (check for validity) */
      //for (i = 0; i < strlen(ftok->key);i++){
      for (i = 0; i < static_cast<int>(strlen(ftok->key)); i++){	// BM
	ind = (int)ftok->key[i];
	if (i == 0)
	  fail = !((((ind >= 45) && (ind <= 57)) && (ind != 47)
		   && (ind != 46)) || (ind == 43));
	else
	  fail = fail || ((ind < 48) || (ind > 57));
      }
      if (!fail)
	/* Set value if valid */
	value = atoi(ftok->key);
      /* Remove start of buffer */
      cptr = ftok;
      ftok = cptr->nxt;
      delete cptr;
    }
  }
  /* Broadcast fail */
  broadcast(&fail);
  if (fail)
    return(false);
  else{
    /* Broadcast value */
    broadcast(&value);
    return(true);
  }
}
//end get_i.c

//file get_r.c
/*
+ get_r - get real value.
     Description:
        The master gets a float value from the buffer of tokens set up
        by get_toks. The token is  removed from the buffer.  The value
        of the  integer is returned  in value and the  routine returns
        true if successful, false otherwise.  If successful, the value
        is broadcast to the other processors.
*/
bool get_r(float &value)
{
  struct token *cptr;
  int i, ind, ndot, ne, ie;
  char temp;
  bool fail;
  if (tmype == 0){
    fail = false;
    /* Get new line if buffer is empty */
    if (ftok == NULL)
      fail = !get_toks();
    if (!fail){
      /* Extract value from buffer (check for validity) */
      ndot = 0;
      ne = 0;
      //for (i = 0; i < strlen(ftok->key); i++){
      for (i = 0; i < static_cast<int>(strlen(ftok->key)); i++){	// BM
	temp = ftok->key[i];
	ind = (int)temp;
	if (i == 0){
	  if (((ind >= 45) && (ind <= 57)) && (ind != 47)){
	    if (temp == '.')
	      ndot ++;
	  }
	  else
	    fail = true;
	}
	else{
	  if (temp == '.')
	    ndot++;
	  else if ((temp == 'e') || (temp == 'E')){
	    ne++;
	    ie = i;
	  }
	  else if ((ne != 0) && (i == ie + 1)){
	    fail = fail || !(((ind >= 48) && (ind <= 57)) || (temp == '+')
			     || (temp == '-'));
	  }
	  else
	    fail = fail || !((ind >= 48) && (ind <= 57));
	}
      }
      if (ndot > 1)
	fail = true;
      if (ne > 1)
	fail = true;
      if (!fail)
	/* Set value if valid */
	value = atof(ftok->key);
      /* Remove token from buffer */
      cptr = ftok;
      ftok = cptr->nxt;
      delete cptr;
    }
  }
  /* Broadcast fail */
  broadcast(&fail);
  if (fail)
    return(false);
  else{
    /* Broadcast value */
    broadcast(&value);
    return(true);
  }
}
//end get_r.c


//file gen_cmd.c
/*
+ gen_cmd - generic commands.
     Description:
        This checks  whether the  command given by  key is one  of the
        generic commands: $, $read, $echo, $noecho, $log.
*/
bool gen_cmd(char cmd[])
{
  int i;
  char temp[SZNAM], key[SZNAM];

  if (strcmp(cmd, "$") == 0){
    /* System command */
    if (!get_n(key)){
      cout << "Error: reading system command" << '\n';
      report();
    }
    else{
      for (i = 0; i < SZNAM; i++)
	temp[i] = ' ';
      strcpy(temp, key);
      while (ftok != NULL){
	if (get_n(key)){
	  strcat(temp, " ");
	  strcat(temp, key);
	}
      }
    }
    /* Only master executes system commands */
    if (tmype == 0)
      sys_cmd(temp);
  }
  else if (strcmp(cmd,"$ECHO") == 0){
    /* Switch on echo */
    echo = true;
  }
  else if (strcmp(cmd,"$NOECHO") == 0){
    /* Switch off echo */
    echo = false;
  }
  else if (strcmp(cmd,"$READ") == 0){
    /* Read commands from .dat file */
    if (!get_n(key)){
      cout << "Error: reading filename." << '\n';
      report();
    }
    else{
      /* Only master opens .dat file */
      if (tmype == 0)
	open_dat(".dat", key, true);
    }
  }
  else if (strcmp(cmd,"$LOG") == 0){
    /* Set up log file */
    if (!get_n(key)){
      cout << "Error: reading log filename \n";
      report();
    }
    else{
      /* Only master writes log files */
      if (tmype == 0){
	sys_cmd(temp);
	if (!open_log(key)){
	  cout << "Error: opening log file \n";
	  report();
	}
      }
    }
  }
  else if (strcmp(cmd,"$NOLOG") == 0){
    /* Master closes log file */
    if (tmype == 0){
      if (logon){
	fclose(logout);
	logon = false;
      }
      else
	cout << "Error: no log file \n";
    }
  }
  else if (strcmp(cmd,"$W") == 0){
    /* Print Warranty */
    if (tmype == 0){
      cout << "This program is distributed in the hope that it will be useful,\n";
      cout << "but WITHOUT ANY WARRANTY; without even the implied warranty \n";
      cout << "of MERCHANTABILITY or FITNESS FOR A PARTICULAR \n";
      cout << "PURPOSE. See the GNU General Public License for more details. \n";
    }
  }
  else if (strcmp(cmd,"$C") == 0){
    /* Print Conditions */
    if (tmype == 0){
      cout << "This program is free software: you can redistribute it and/or modify \n";
      cout << "it under the terms of the GNU General Public License as published by \n";
      cout << "the Free Software Foundation, either version 3 of the License, or \n";
      cout << "any later version. \n";
    }
  }    
  else
    /* Command unknown here - pass it back to caller */
    return(false);

  return(true);
}
//end gen_cmd.c

//file sys_cmd.c
/*
+ sys_cmd - system commands.
     Description:
        This spawns a  process to execute system commands  such as ls,
        cp etc. It should only be executed by the master.
*/
void sys_cmd(char cmd[])
{
   pid_t pid, pidr;
   int status;
   if ((pid=fork()) < 0) {
     cout << "pid_fork < 0 \n";
      perror("system.fork");
      fprintf(stdout, "error:  failed to create process\n");
   } else if (pid == 0) {
      const char *arg0,*arg1;
      char line[SZNAM];
      arg0 = "csh";
      arg1 = "-c";
      strcpy(line, cmd);
      execlp("/bin/csh", arg0, arg1, line, (char*)NULL);
      fprintf(stdout, "error: could not execlp\n");
      exit(0);
   } else {
      while (1) {
         if ((pidr=waitpid(-1, &status, 0)) == -1) {
            if (errno != EINTR) {
               printf("pongz_:  error is not EINTR\n");
               break;
            }
         } else if (pidr != pid) {
         } else {
            break;
         }
      }
   }
}
//end sys_cmd.c

//file open_dat.c
/*
+ open_dat - open .dat file for reading.
     Description:
        This opens the file key.dat  for reading. The commands are now
        read from this file until eof. This is done by setting readf =
        true so that get_toks.c knows what to do.
*/
bool open_dat(const char tail[], char name[], bool techo)
{
  bool success;
  char *fname;
  fname = new char[sizeof(char)*strlen(name)+5];
  strcpy(fname, name);
  strcat(fname, tail);

  if (idat > NDAT - 2){
    cout << "error: too many nested input files \n";
    return(false);
  }
  idat++;
  if (!(datin[idat] = fopen(fname, "r"))){
    cout << "Error: unable to open input file \n";
    idat--;
    readf = !(idat <= 0);
    success = false;
  }
  else{
    oftok[idat] = ftok;
    ftok = NULL;
    readf = true;
    success = true;
  }
  delete [] fname;
  echo = techo;
  return(success);
}
//end open_dat.c

//file open_log.c
/*
+ open_log - open log file for writing.
     Description:
        This opens the file key.log  for writing. The commands are now
        written to this  file until the command nolog  is issued. This
        is done by setting log to  true so that get_toks knows what to
        do. An error is returned if the file exists.
*/
bool open_log(char name[])
{
  bool success;
  char *fname;
  fname = new char[sizeof(char)*strlen(name)+5];
  strcpy(fname, name);
  strcat(fname, ".log");
  if ((logout=fopen(fname, "r"))){
    cout << "Error: file exists \n";
    fclose(logout);
    success = false;
  }
  else if (!(logout = fopen(fname, "w"))){
    cout << "Error: unable to open log file \n";
    success = false;
  }
  else{
    logon = true;
    success = true;
  }
  delete [] fname;
  return(success);
}
//end open_log.c

//file report.c
/*
+ report - close input files if error.
     Description:
        This should  be called to close  the input (.dat)  files if an
        error is encountered. Only the master executes this.
*/
void report()
{
  int i;
  struct token *cptr;
  if (tmype != 0)
    return;
  /* Clear input buffer */
  while (ftok != NULL){
    cptr = ftok;
    ftok = cptr->nxt;
    delete cptr;
  }
  ftok = NULL;
  /* Check for open input files */
  if (readf){
    for (i = 0; i <= idat; i++){
    /* Close all input files */
      if (datin[idat] != NULL)
	fclose(datin[idat]);
      datin[idat] = NULL;
      /* Clear old input buffers */
      while (oftok[idat] != NULL){
	cptr = oftok[idat];
	oftok[idat] = cptr->nxt;
	delete cptr;
      }
      oftok[idat] = NULL;
    }
    idat = -1;
    readf = false;
    echo = false;
  }
  return;
}
//end report.c
