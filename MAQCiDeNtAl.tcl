#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" ${1+"$@"}

proc gene2reg2content {chosendb} {
	foreach db $chosendb {
		if {![file isfile /complgen/bin/MastrDesign/db/compressed_reg_$db.tsv.lz4]} {puts "Error: The reference file for $db does not exist or the path/name has been changed."; exit}
		cg select -q "\$type shares \"CDS RNA\"" /complgen/bin/MastrDesign/db/compressed_reg_$db.tsv.lz4 ${db}_content.tsv
		cg select -s - ${db}_content.tsv sorted_${db}_content.tsv } }

proc geneproc {dblist genelist} {
	#Checking the databases and using conversion of identifier in the case of ensGene and knownGene.
	# For use in ensGene and knownGene the gene identifier needs to be substituted by the proper format. A table is available for this purpose. However, a genename maps to multiple IDs.
	# Using the option oneof to select matches out of a comma separated list. Requires quotes.	
	# The output file still requires a name field with unique content.
	set leaderfile {chromosome	begin	end	name}
	set filenameleaderfile "pre.processedtargets.tsv"
	set fileIdleaderfile [open $filenameleaderfile "w"]
	puts $fileIdleaderfile $leaderfile
	close $fileIdleaderfile
	foreach usergenename $genelist {
		set usergenename2 \"$usergenename\"
		foreach db $dblist {
			if {$db in [list ensGene knownGene]} {
				set options [idconversion $db $usergenename [llength $dblist]] 
				cg select -q "oneof(\$transcript,\"$options\")" sorted_${db}_content.tsv selectedgene_${db}_$usergenename.tsv
				set strand [lindex [exec cg select -q "oneof(\$name,\"[lindex $options 0]\")" /complgen/bin/MastrDesign/db/strand_${db}.tsv] end]				}
			if {$db in [list refGene gencode]} {
				cg select -q "\$gene==$usergenename2" sorted_${db}_content.tsv selectedgene_${db}_$usergenename.tsv
				if {[exec wc -l < selectedgene_${db}_$usergenename.tsv] == "1" && $dblist == "1"} {puts "No gene matches $usergenename in ${db}. You might try again in another database."; exit}
				if {[exec wc -l < selectedgene_${db}_$usergenename.tsv] == "1" && $dblist != "1"} {puts "No match with $usergenename in ${db}. " ; break}
				set strand [lindex [exec cg select -q "\$name == \"$usergenename\"" /complgen/bin/MastrDesign/db/strand_${db}.tsv] end]} 	
			exec cg regjoin selectedgene_${db}_$usergenename.tsv > joined_selected_${db}_$usergenename.tsv
			cg select -f {chromosome begin end} joined_selected_${db}_$usergenename.tsv > extended_joined_selected_${db}_$usergenename.tsv
			exec sed -e {s/chr//g} -e {s/omo/chromo/1} extended_joined_selected_${db}_$usergenename.tsv > chrOkay_joined_selected_${db}_$usergenename.tsv }
	joinselecteddb $usergenename $dblist
	if {[file isfile joined_all_$usergenename.tsv]} {
	if {$strand == "+"} {
		cg select -f "* {name=concat(\"$usergenename\",\"_\",\$ROW+1)}" joined_all_$usergenename.tsv > named_joined_all_$usergenename.tsv }
	if {$strand == "-"} {
		set offset [exec wc -l < joined_all_$usergenename.tsv]
		cg select -f "* {name=concat(\"$usergenename\",\"_\",$offset-(\$ROW+1))}" joined_all_$usergenename.tsv > named_joined_all_$usergenename.tsv }
	exec tail -n +2 named_joined_all_$usergenename.tsv >> pre.processedtargets.tsv } } }		

proc idconversion {database gene dbnum} {
	#For use in ensGene and knownGene the gene identifier needs to be substituted by the proper format. A table is available for this purpose. However, a gene name maps to multiple IDs.
	#To be used in the foreach loop of dbquery, with if statement for database ensGene and knownGene with database being a single db, either ensGene or knownGene
	if {![file isfile /complgen/refseq/hg19/extra/trans_hg19_${database}2genename.tsv]} {puts "Error: The identifier conversion file for $database does not exist or the path/name has been changed."; exit}
	cg select -q "\$genename == \"$gene\"" /complgen/refseq/hg19/extra/trans_hg19_${database}2genename.tsv > ${database}_temp1_$gene.tsv
	if {[exec wc -l < ${database}_temp1_$gene.tsv] == "1" && $dbnum == "1"} {puts "Error: No gene matches $gene in $database. You might try again in another database."; exit}
	if {[exec wc -l < ${database}_temp1_$gene.tsv] == "1" && $dbnum != "1"} {puts "No match with $gene in $database. "}
	exec cut -f1 ${database}_temp1_$gene.tsv > ${database}_temp2_$gene.tsv
	exec tail -n +2 ${database}_temp2_$gene.tsv > ${database}_temp3_$gene.tsv
	set fp [open ${database}_temp3_$gene.tsv r]
	set tempvar [read $fp]
	close $fp
	set tempvar2 [concat $tempvar]
	set tempvar3 [join  $tempvar2 "\",\""]
	if {[exec wc -l ${database}_temp1_$gene.tsv] != "1 ${database}_temp1_$gene.tsv"} {puts "When converted to $database transcript IDs, $gene corresponds to \"$tempvar3\""}
	return $tempvar3 }
	
proc joinselecteddb {gene databases}	{
	set filenamelist [list]
	foreach element $databases {
		set file chrOkay_joined_selected_${element}_${gene}.tsv 
		if {[file isfile $file]} {lappend filenamelist $file} }
	if {[llength $filenamelist] > 0} {
	cg cat -c 0 {*}$filenamelist > cat_$gene.tsv
	cg select -s - cat_$gene.tsv sorted_$gene.tsv
	cg regjoin sorted_$gene.tsv > joined_all_$gene.tsv } } 

proc concatsidebyside {tojoinlist} {
	set former 0
	set joinlist []
	foreach seq $tojoinlist {
		if {$former != 0} {
		append former $seq
		lappend joinlist $former
			}
		set former $seq 
		}
	return $joinlist }

proc maketargetfile {joinidlist seqs gene} {
	set joinlist []
	foreach joinid $joinidlist seq $seqs {
		if {$joinid != ""} {
		set targetspec ""
		append targetspec [expr [string length $seq] - 1] ",4"
		lappend joinlist $joinid $targetspec} }
	set outcontent [join $joinlist "\n"]
	set outfile "$gene.exons.fas.target"
	set fileId [open $outfile "w"]
	puts -nonewline $fileId $outcontent
	close $fileId
	}

proc joinneighbours { fasfile gene } {
	set fp [open $fasfile r]
	set fasta [read $fp]
	close $fp
	set fastalist [split $fasta "\n"]
	set ids []
	set seqs []                                                        
	foreach {i j} $fastalist {if {$j != ""} {lappend seqs $j}}
	foreach {i j} $fastalist {if {$i != ""} {lappend ids $i}}
	set joinseqlist [concatsidebyside $seqs]
	set joinidlist [concatsidebyside $ids]
	maketargetfile $joinidlist $seqs $gene
	set joinlist []
	foreach seq $joinseqlist id $joinidlist {lappend joinlist $id $seq}
	set outcontent [join $joinlist "\n"]
	set outfile "$gene.$fasfile"
	set fileId [open $outfile "w"]
	puts -nonewline $fileId $outcontent
	close $fileId
	}
	
# Still need to fix declaration of variables
set dblist [list refGene]
set gene [lindex $args 1]
set snpfreq 0.01

file mkdir temp_MAQCiDeNtAl
cd temp_MAQCiDeNtAl

gene2reg2content $dblist
geneproc $dblist $gene

cg genome_seq -n $snpfreq -p Common -r "N" -i "name" pre.processedtargets.tsv /complgen/refseq/hg19/ exons.fas.masked
cg genome_seq -i "name" pre.processedtargets.tsv /complgen/refseq/hg19/ exons.fas

foreach file [list exons.fas exons.fas.masked] {joinneighbours $file $gene}

file copy $gene.exons.fas.target ../
file copy $gene.exons.fas  ../
file copy $gene.exons.fas.masked ../
cd ../
file delete -force temp_MAQCiDeNtAl

puts "Design finished"