# load the trajectory
mol load gro re100_r1.gro
mol addfile re100_r1.xtc type xtc first 0 last -1 waitfor all

# open a file for writing the data to
set OUTPUT [open "CONTACTS.dat" w]
set nf [molinfo top get numframes]

set contacts [atomselect top "name CA and within 5.0 of resname REMP"]

# loop over all frames
for {set i 0} {$i < $nf} {incr i} {

    $contacts frame $i
    $contacts update
    set molnum [$contacts num]
 
    # write out the frame number and number of lipid atoms to file  (the bit in the quotes formats the data nicely)
    puts "\t progress: $i/$nf"
#   puts $OUTPUT " $i $contacts"
    puts $OUTPUT "$i $molnum"
}
puts "\t progress: $nf/$nf...... and Done"
puts "output: CONTACTS.dat"
# close the file        
close $OUTPUT

