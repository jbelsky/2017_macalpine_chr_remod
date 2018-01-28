package chip_seq_analysis;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import genomics_functions.TF;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;


public class ChIPSeqFeatureSignal {

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
		
		// Set the span
		int span = Integer.parseInt(args[5]);
		
		/////////////////////////////////////////////////////////////////////////////
					
		// Open the output buffer
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
				
		// Write the output feature header
		output.write(ChIPSeqFeatureDensitySignal.get_output_header(win, span));
		
		// Read in the ARS
		ArrayList<TF> tf_list = TF.read_in_tf_list(input_feature_file);
		
		// Open the BAM File
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));
		
		// Iterate through each feature
		for(int i = 0; i < tf_list.size(); i++){
			
			// Get the TF
			TF feature = tf_list.get(i);
						
			// Set up the storage array
			int[] storage_array = new int[2*win + 1];
			
			// Set up the SAM Record Iterator
			SAMRecordIterator bam_itr = bam_file.queryOverlapping(feature.getChr(),
																  feature.getPos() - win,
																  feature.getPos() + win
															     );
			

			
			// Iterate through the bam_itr
			while(bam_itr.hasNext()){

				// Initialize the read_end_pos
				int read_end_pos = -1;
				
				// Get the read
				SAMRecord read = bam_itr.next();
				
				// Get the read_end_pos
				if(read.getReadNegativeStrandFlag()){
					read_end_pos = read.getAlignmentEnd() - shift;
				}else{
					read_end_pos = read.getAlignmentStart() + shift;
				}
				
				// Shift to get the relative pos
				read_end_pos -= (feature.getPos() - win);
				
				// If the read_end_pos is within the storage_array, add it to the storage_array
				if(read_end_pos >= 0 && read_end_pos < storage_array.length){
					storage_array[read_end_pos]++;
				}
				
			}
				
			// Close the bam_itr
			bam_itr.close();
	
			// Write the subnuc output
			tf_list.get(i).write_reads_output(output, storage_array, span);
			
		}
		
		// Close the bam_file
		bam_file.close();
		
		// Close the output buffer
		output.close();
		
	}

}
