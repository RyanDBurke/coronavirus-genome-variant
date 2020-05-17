




FASTAFILE *
OpenFASTA(char *seqfile)
{
  FASTAFILE *ffp;

  ffp = malloc(sizeof(FASTAFILE));
  ffp->fp = fopen(seqfile, "r");              /* Assume seqfile exists & readable!   */
  if (ffp->fp == NULL) { free(ffp); return NULL; } 
  if ((fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)) == NULL)
    { free(ffp); return NULL; }
  return ffp;
}

int
ReadFASTA(FASTAFILE *ffp, char **ret_seq, char **ret_name, int *ret_L)
{
  char *s;
  char *name;
  char *seq;
  int   n;
  int   nalloc;
  
  /* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline.
   */
  if (ffp->buffer[0] != '>') return 0;    

  /* Parse out the name: the first non-newline token after the >
   */
  s  = strtok(ffp->buffer+1, "\n");
  name = malloc(sizeof(char) * (strlen(s)+1));
  strcpy(name, s);

  /* Everything else 'til the next descline is the sequence.
   * Note the idiom for dynamic reallocation of seq as we
   * read more characters, so we don't have to assume a maximum
   * sequence length.
   */
  seq = malloc(sizeof(char) * 128);     /* allocate seq in blocks of 128 residues */
  nalloc = 128;
  n = 0;
  while (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp))
    {
      if (ffp->buffer[0] == '>') break;	/* a-ha, we've reached the next descline */

      for (s = ffp->buffer; *s != '\0'; s++)
	{
	  if (! isalpha(*s)) continue;  /* accept any alphabetic character */

	  seq[n] = *s;                  /* store the character, bump length n */
	  n++;
	  if (nalloc == n)	        /* are we out of room in seq? if so, expand */
	    {			        /* (remember, need space for the final '\0')*/
	      nalloc += 128;
	      seq = realloc(seq, sizeof(char) * nalloc);
	    }
	}
    }
  seq[n] = '\0';

  *ret_name = name;
  *ret_seq  = seq;
  *ret_L    = n;
  return 1;
}      

void
CloseFASTA(FASTAFILE *ffp)
{
  fclose(ffp->fp);
  free(ffp);
}
