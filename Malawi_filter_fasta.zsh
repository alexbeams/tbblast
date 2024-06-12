awk '
  BEGIN { 
    # Read the names from the file into an array
    while ((getline < "May24/Malawi_TB_with_dates.txt") > 0) {
      names[$1] = 1
    }
    # Reset the record separator to ">" for FASTA format
    RS=">"
    # Set the output record separator
    ORS=""
  }
  # For each sequence entry
  NR>1 {
    # Extract the sequence name
    split($0, header, "\n")
    sequence_name = header[1]
    # If the sequence name is found in the list of desired names, print the sequence
    if (names[sequence_name]) {
      print ">" $0
    }
  }
' Malawi_final.fasta > Malawi_final_filtered.fasta
