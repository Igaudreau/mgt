/* parc_xdataf: fast xdatar */
void parc_xdataf(int *npe, int outs[], float outbuf[], int ins[], 
  float inbuf[], int *ios)
{
  int limpe, ipe, nouts;
  struct proc_list{
    int proc;
    struct proc_list *nxt;
  };
  struct proc_list *plist, *start;
  limpe = *npe;
  /* Create receive processor list and number of elements */
  start = NULL;
  for (ipe = 0; ipe < limpe; ipe++){
    if (outs[ipe] != 0){
      nouts = outs[ipe];
      plist = (struct proc_list *)
        malloc(sizeof(struct proc_list));
      plist->proc = ipe;
      plist->nxt = start;
      start = plist;
    }
  }

}
