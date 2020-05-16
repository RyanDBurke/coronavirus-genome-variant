public class fmmap {

    /* 
    Program Driver
        * reads command-line input and executes accordingly
    */
    public static void main(String[] args) {

        /* 
        determine if command was entered
        with a proper amount of arguments
        */
        if (!(args.length > 1)) {
            System.out.println("Please enter a valid command: index or align");
        }

        /* read from command-line */
        String command = args[0];

        /* determine which command to execute */
        if (command.equals("index")) {
            System.out.println("index");
        } else if (command.equals("align")) {
            System.out.println("align");
        } else {
            System.out.println(command + " is not a valid command.");
        }
    }
}