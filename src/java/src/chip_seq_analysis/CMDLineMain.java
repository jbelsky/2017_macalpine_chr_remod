package chip_seq_analysis;

import java.io.IOException;
import java.util.Arrays;

public class CMDLineMain {

	public static void main(String[] args) throws IOException {

		// Get the program to execute
		String java_executable_program = args[0];
		
		// Get the command line arguments
		String[] cmd_line_args = Arrays.copyOfRange(args, 1, args.length);
		
		// Choose the program to run based on the args[0];
		if(java_executable_program.equals("ChIPSeqFeatureDensitySignal")){
			System.out.println("Executing ChIPSeqFeatureDensitySignal!");
			ChIPSeqFeatureDensitySignal.main(cmd_line_args);
		}else if(java_executable_program.equals("ChIPSeqFeatureSignal")){
			System.out.println("Executing ChIPSeqFeatureSignal!");
			ChIPSeqFeatureSignal.main(cmd_line_args);
		}else{
			System.out.println("Did not specify either 'ChIPSeqFeatureDensitySignal' " +
							   "or 'ChIPSeqFeatureSignal' as first argument!");
		}

	}

}
