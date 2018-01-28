package chip_seq_analysis;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import genomics_functions.DensityEst;
import genomics_functions.TF;

import net.sf.samtools.SAMFileReader;


public class ChIPSeqFeatureDensitySignal {

	public static String get_output_header(int win_input, int span_input){
	
		// Set up the string
		String header_str = "name,chr,pos,strand,";
	
		String sep = ",";
	
		// Iterate through each position
		for(int i = -win_input; i <= win_input; i += span_input){
			
			if(i == win_input)
				sep = "\n";
			
			// Add the position to the string
			header_str += (i + sep);
	
		}
	
		// Return the header_str
		return(header_str);
		
	}
		
	public static void main(String[] args) throws IOException {

		// Set the bam file
		String bam_file_name = args[0];
		
		// Set the input feature name
		String input_feature_file = args[1];
		
		// Set the output file name
		String output_file_name = args[2];
				
		// Set the window
		int win = Integer.parseInt(args[3]);
		
		// Set the shift
		int shift = Integer.parseInt(args[4]);
		
		// Set the bw
		int bw = Integer.parseInt(args[5]);
		
		// Set the span
		int span = Integer.parseInt(args[6]);
		
		/////////////////////////////////////////////////////////////////////////////
					
		// Open the output buffer
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
				
		// Write the output feature header
		output.write(get_output_header(win, span));
		
		// Read in the ARS
		ArrayList<TF> tf_list = TF.read_in_tf_list(input_feature_file);
		
		// Open the BAM File
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));
		
		// Iterate through each feature
		for(int i = 0; i < tf_list.size(); i++){
				
			// Get the subnucleosome coverage density over the region
			double[] feat_subnuc_cov = 
				DensityEst.get_chip_seq_coverage_density(
					bam_file, tf_list.get(i), win, shift, bw, 1
				    );
	
			// Write the subnuc output
			tf_list.get(i).write_reads_output(output, feat_subnuc_cov, span);
			
		}
			
		// Close the output buffer
		output.close();
		
	}

}
