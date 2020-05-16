public class fmmap {

    /* path to fasta file containing reference sequence */
    public String seq;

    /* path to output file */
    public String out;
    
    /* pass in our paths to the fasta file and output file */
    public fmmap(String seq, String out) {
        this.seq = seq;
        this.out = out;
    }

    /* 
    index command 
        * align:
            * boolean to determine if we will eventually align
        * this takes in a reference sequence and will write it's FMIndex to output
    */
    public void index(boolean align) {

        /* our output file we will write to */
        out = this.out

        /* determine file-type extension */
    }

    /* 
    Program Driver
        * reads command-line input and executes accordingly
        * java fmmap index ref_seq output shouldWeAlign
        * java fmmap align ref_seq reads output
    */
    public static void main(String[] args) {

        /* determine if command was entered with a proper amount of arguments */
        if (!(args.length > 1)) {
            System.out.println("Please enter a valid command: index or align");
        }

        /* read from command-line */
        String command = args[0].toLowerCase();

        /* determine which command to execute */

        /* index */
        if (command.equals("index")) {
            
            /* path to fasta */
            String seq = args[1];

            /* path to output file */
            String output = args[2];

            /* create our fmmap */
            fmmap FM = new fmmap(seq, output);

            /* call index command */
            FM.index(false);
        } 
        
        /* align */
        else if (command.equals("align")) {
            System.out.println("align");
        }
        
        /* no valid commands entered */
        else {
            System.out.println(command + " is not a valid command.");
        }
    }
}