#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" ${1+"$@"}

# MastrDesign script V3
# Code by Wouter De Coster, implementing an R script of Arne De Roeck
# Many thanks to Peter De Rijk & Bart Smets
# Script to completely automate the design of a MASTR assay. The script also provides a file ready for use as a custom track in the UCSC Genome Browser
# The script (in coord mode)is now produced for SNPs and should work fine for small indels. However, I'm not too sure about indels larger than a few bases.

set testproject true
set mastrtype undefined
set kit undefined
set inputfile undefined
set extend undefined
set targetsize undefined
set gccontent undefined
set snpfreq undefined
set varlist ""
set knownGene NA
set refGene NA
set ensGene NA
set gencode NA
set dblist [list]
set content [list]
set user [exec whoami]
set server [info hostname]
set mirnalist3 undefined
set time1 [clock format [clock seconds] -format %T]
set genelist [list]
set covered2 undefined
set secretgene false
set targetnumber undefined
set wcnumber undefined
set wcuniquenames undefined
set highgc NA
set lowgc NA
set elements false
set screenshot false
set settings interactive
set gcchecker off
set mastrtype_options [list c coords g gene m mirna]
set db_options [list yes y no n]

proc checkX {} {exec /complgen/bin/anaconda/bin/python2.7 /complgen/bin/MastrDesign/lib/Xrun.py }

proc helpfunction {} {catch {exec /complgen/bin/anaconda/bin/python2.7 /complgen/bin/MastrDesign/lib/helpopener.py}	}

proc makeproject {} {
	pack [labelframe .questions.project -text "Provide a unique projectname for your design (not containing spaces or strange symbols)."] -side top -fill x -pady 10
	pack [entry .questions.project.projectquestion -textvariable projectquestion] -fill x -side left
	pack [button .questions.project.confirmproject -text "Confirm" -command {
	set projectname $projectquestion
	destroy .questions.project
	destroy .questions.used
	}] -side right
	tkwait variable projectname }

proc prompt {question} {
	puts "$question"
	flush stdout
	gets stdin}

proc gene2reg2content {chosendb chosencontent} {
	foreach db $chosendb {
		if {![file isfile /complgen/bin/MastrDesign/db/compressed_reg_$db.tsv.lz4]} {puts "Error: The reference file for $db does not exist or the path/name has been changed."; exit}
		cg select -q "\$type shares \"$chosencontent\"" /complgen/bin/MastrDesign/db/compressed_reg_$db.tsv.lz4 ${db}_content.tsv
		cg select -s - ${db}_content.tsv sorted_${db}_content.tsv } }

proc joinselecteddb {gene databases}	{
	set filenamelist [list]
	foreach element $databases {
		set file chrOkay_joined_selected_${element}_${gene}.tsv
		if {[file isfile $file]} {lappend filenamelist $file} }
	if {[llength $filenamelist] > 0} {
	cg cat -c 0 {*}$filenamelist > cat_$gene.tsv
	cg select -s - cat_$gene.tsv sorted_$gene.tsv
	cg regjoin sorted_$gene.tsv > joined_all_$gene.tsv } }

proc idconversion {database gene dbnum} {
	#For use in ensGene and knownGene the gene identifier needs to be substituted by the proper format. A table is available for this purpose. However, a gene name maps to multiple IDs.
	#To be used in the foreach loop of dbquery, with if statement for database ensGene and knownGene with database being a single db, either ensGene or knownGene
	if {![file isfile /complgen/refseq-0.11.0/hg19/extra/trans_hg19_${database}2genename.tsv]} {puts "Error: The identifier conversion file for $database does not exist or the path/name has been changed."; exit}
	cg select -q "\$genename == \"$gene\"" /complgen/refseq-0.11.0/hg19/extra/trans_hg19_${database}2genename.tsv > ${database}_temp1_$gene.tsv
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

proc handlevars {inputfile} {
	exec sed {/^@/d} $inputfile > cleanedfile.tsv
	exec sed {/^@/!d} $inputfile | tr -d {@} > secretgenefile.tsv
	if {[exec wc -l < secretgenefile.tsv] != 0} {set secretgene true} {set secretgene false}
	#Under the rare condition that someone would like to include both genes and variants in your design, you can smuggle in genes by including them in your targetfile preceded by an @
	if {$secretgene == "true"} {
		set fp3 [open "secretgenefile.tsv" r]
		set presecretgenelist [read $fp3]
		close $fp3
		set genelist [join $presecretgenelist]
		return [list $secretgene $genelist] } }

proc geneproc {dblist genelist} {
	#Checking the databases and using conversion of identifier in the case of ensGene and knownGene.
	#To include splice sites, the targets are extended in both directions with 20bp.
	#For use in ensGene and knownGene the gene identifier needs to be substituted by the proper format. A table is available for this purpose. However, a genename maps to multiple IDs.
	#Using the option oneof to select matches out of a comma separated list. Requires quotes.
	#The output file still requires a name field with unique content.
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
			cg select -f {chromosome {begin=$begin - 20} {end=$end + 20}} joined_selected_${db}_$usergenename.tsv > extended_joined_selected_${db}_$usergenename.tsv
			exec sed -e {s/chr//g} -e {s/omo/chromo/1} extended_joined_selected_${db}_$usergenename.tsv > chrOkay_joined_selected_${db}_$usergenename.tsv }
	joinselecteddb $usergenename $dblist
	if {[file isfile joined_all_$usergenename.tsv]} {
	if {$strand == "+"} {
		cg select -f "* {name=concat(\"$usergenename\",\"_\",\$ROW+1)}" joined_all_$usergenename.tsv > named_joined_all_$usergenename.tsv }
	if {$strand == "-"} {
		set offset [exec wc -l < joined_all_$usergenename.tsv]
		cg select -f "* {name=concat(\"$usergenename\",\"_\",$offset-(\$ROW+1))}" joined_all_$usergenename.tsv > named_joined_all_$usergenename.tsv }
	exec tail -n +2 named_joined_all_$usergenename.tsv >> pre.processedtargets.tsv } } }

set Xactive [checkX]

#Passing 'h' as argument to the script will create a python script that will open the help page. The rest of the script will be aborted.
#Passing 'd' as an argument to the script will use default settings and create a random projectnumber.
#Passing 's' as argument to the script will create a python script and activate the screenshot part at the end of the script.
#Passing 'g' as an argument to the script will use a list of possible values to split on for gc content, as such an optimal distribution can be chosen.
if {"h" in $args} {if {$Xactive == "Active"} {helpfunction} {puts "Please enable X11 tunnelling by starting Xming"}}
if {"s" in $args} {set screenshot True}
if {"d" in $args} {
	set settings default
	set testproject default
	set projectname [expr int(rand()*1000)]
	file mkdir temp_$projectname
	cd temp_$projectname
	set mastrtype g
	set kit 3
	set elements Okay
	set genelist [list PGAM5 BTBD10]
	set content CDS
	set dblist [list refGene ensGene]
	set refgene y
	set knownGene n
	set ensGene y
	set gencode n
	set targetsize 350
	set gcsplit off
	set gccontent no
	set extend "full"
	set snpfreq 0.01
	puts "Using default settings. Your randomly generated projectname is $projectname." }
if {"g" in $args} {
	set gcchecker on
	set gcsplit on
	set gccontent [list 50 55 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80]
	puts "Checking the amount of amplicons in highgc and lowgc using the following cut-offs: $gccontent"}
if {"i" in $args || $Xactive == "Inactive"} {

if {$settings != "default"} {
	puts "Hi $user, I think this is a great day to design a MASTR assay.\nThis is MastrDesign V3.\nType 'help' as project name to view the latest documentation page.\nPlease report problems, ideas, hopes and dreams to wouter.decoster@molgen.vib-ua.be"}

### QUESTION BLOCK ###
# Without the d argument, the user will be asked interactively to provide input for the script.
# The project name should be unique and this is tested by checking whether the Genome Browser Track file already exists in this directory.
# All temporary files and scripts are localized in a temp folder which is deleted at the very end after copying the output files, unless the 't' argument is passed to the script.
while {$testproject == "true" } {
	set projectname [prompt "Please give a unique project name (not containing spaces or special characters)."]
	set dir [pwd]
	if {[file isfile $dir/Genome_Browser_Track_$projectname.txt] && ![string is alnum $projectname]} {puts "Projectname invalid or already used in this directory!"} {set testproject false} }
if {$projectname == "help"} {if {$Xactive == "Active"} {helpfunction} {puts "Please enable X11 tunnelling by starting Xming"} ; exit}
file mkdir temp_$projectname
cd temp_$projectname

while {$mastrtype ni $mastrtype_options} {
	set mastrtype [prompt "What type of MASTR would you like to design today?\nBased on a gene or on genomic coordinates? (Options: \[c\]oords \[g\]ene \[m\]irna)"] }

while {$kit ni [list 2 3]} {
	set kit [prompt "Which MiSeq kit are you going to use? (Options: 2 3)"]}
	if {$kit == "2"} {set maxsize 400}
	if {$kit == "3"} {set maxsize 500}

if { $mastrtype == "coords" || $mastrtype == "c"}	{
	while {![file isfile $inputfile]} {
		set inputfile [prompt "Which are the variants for your MASTR design?\nYour file is tab separated and requires the fields 'chromosome begin end name'.\nNames need to be a unique value. Provide absolute path."]
		set file [file normalize $inputfile] }
	# Extension of target region desired? Full means entire exon disregarding size, limited means only one amplicon per target and as such perhaps an incomplete exon.
	while { $extend ni [list full f limited l no n]} {
		set extend [prompt "Should the target around coding variants be extended to the exon? (Options: \[f\]ull \[l\]imited \[n\]o)"]
		if { $extend in [list full f]} 						{puts "The tool will create amplicons for the entire coding exon. If necessary amplicons will be split."
		} 	elseif { $extend in [list limited l]} 			{puts "The tool will create extend targets to coding sequences around the variants.\nIf the exon is too large it will be included partially."
		} 	elseif { $extend in [list no n]} 				{puts "The tool will process the targets for a SNP MASTR and not include additional coding sequence."; set targetsize 32 } }
		lassign [handlevars $file] secretgene genelist
		if {$secretgene == "true"} {puts "Activating mixed design with $genelist"} }

if {$mastrtype in [list gene g]} {
	set extend "full"
	while {[llength $genelist] == 0} {
		set genelist [prompt "What is/are your favourite gene(s) to include in your gene-MASTR?\nProvide HGNC name(s), case sensitive, space separated. (e.g. LRRK2) or specify path to gene file."]
		if {[file isfile $genelist]} {
			set fp4 [open "$genelist" r]
			set pregenelist [read $fp4]
			close $fp4
			set genelist [join $pregenelist]
			puts "You selected $genelist"
			} else {
			puts "You selected $genelist"} }
	}
if {$mastrtype in [list gene g] || $secretgene == "true"} {
	while {$elements == "false"} {
		set contentoptions [list CDS UTR RNA intron upstream downstream]
		set content [prompt "Which gene elements would you like to include? Provide space separated list, case sensitive!.\n(Options: $contentoptions) (upstream and downstream means 2kb extension)"]
		set elementokay 0
		foreach element $content {
			if {$element in $contentoptions} {incr elementokay} {puts "I don't recognize $element"} }
		if {$elementokay == [llength $content]} {set elements Okay} }

	while {[llength $dblist] == "0"} {
		while {$refGene ni $db_options} {
			set refGene [prompt "Would you like to take gene annotation from RefGene?. (Options: \[y\]es \[n\]o)"]
			if { $refGene == "yes" || $refGene == "y"} {lappend dblist refGene}  }
		while {$knownGene ni $db_options} {
			set knownGene [prompt "Would you like to take gene annotation from knownGene (UCSC)?. (Options: \[y\]es \[n\]o)"]
			if {$knownGene == "yes" || $knownGene == "y"} {lappend dblist knownGene}  }
		while {$ensGene  ni $db_options} {
			set ensGene [prompt "Would you like to take gene annotation from ensGene (Ensembl)?. (Options: \[y\]es \[n\]o)"]
			if { $ensGene == "yes" || $ensGene == "y" } {lappend dblist ensGene}  }
		while {$gencode  ni $db_options} {
			set gencode [prompt "Would you like to take gene annotation from gencode (Encode)?. (Options: \[y\]es \[n\]o)"]
			if { $gencode == "yes" || $gencode == "y"} {lappend dblist gencode}  } }
	if {$settings != "default"} {puts "You selected $dblist to take gene annotation from."} }

if {$mastrtype in [list miRNA m]} {
	set extend "full"
	if {![file isfile /complgen/refseq/hg19/reg_hg19_mirbase.tsv]} {puts "Error: The reference file for mirbase does not exist or the path/name has been changed."; exit}
	set mirnalist [prompt "What is/are your favourite mirna gene(s) to include in your mirna-MASTR?\nProvide mirbase name(s), case sensitive, space separated. (e.g. mir-1273d)"]
	set mirnalist2 [join $mirnalist "\",\"hsa-"]
	set mirnalist3 \"hsa-$mirnalist2\"
	puts "You selected $mirnalist3 for you mirna-MASTR." }

if {$extend ni [list no n] || $secretgene == "true"} {
	while {$targetsize > $maxsize || $targetsize < 100} {
		set targetsize [prompt "How long (bp) do you want your target sequences to be maximally? Limited between 100bp and ${maxsize}bp. (Provide number)"]
		if { $targetsize >= 100 && $targetsize <= $maxsize} {puts "You have chosen a maximum target size of $targetsize bp." } } }

	if {$gcchecker == "off"} {
	while { $gccontent !="no" && ($gccontent > 100 || $gccontent <= 0)} {
		set gccontent [prompt "On which gc percentage would you like to split your plexes? (Options: 'number' no)"]}
		if { $gccontent == "no"} {set gcsplit off }	{puts "The tool will split outputfiles on a gc content of $gccontent %."; set gcsplit on} }

	while { $snpfreq > 1 || $snpfreq < 0} {
		set snpfreq [prompt "On which frequency would you like to mask your fasta files? (Options: decimal number e.g. 0.01)"] }

if {$settings != "default"} {puts "Fasta files will be masked for snps with a frequency of $snpfreq.\nThe script now starts the processing procedure.\nMeanwhile, help yourself to a nice cup of coffee or tea." }
}

if {"i" ni $args && "d" ni $args && $Xactive == "Active"} {
package require Tk
wm protocol . WM_DELETE_WINDOW {
    if {[tk_messageBox -message "Quit?" -type yesno] eq "yes"} { exit } }
frame .questions -padx 50
pack .questions -side top
wm title . "MastrDesign"

pack [labelframe .questions.main] -side top
pack [button .questions.main.exit -text "Help" -command {helpfunction}] -side right  -padx 20 -pady 20
pack [label .questions.main.welcome -text "Hi $user, I think this is a great day to design a MASTR assay.\nPlease report problems, ideas, hopes and dreams to wouter.decoster@molgen.vib-ua.be"] -side left

while {$testproject == "true" } {
	makeproject
	set dir [pwd]
	if {[file isfile [pwd]/Genome_Browser_Track_$projectname.txt] || ![string is alnum $projectname]} {pack [label .questions.used -text "Projectname \"$projectname\" invalid or already used in this directory!" -background red]} {set testproject false} }
file mkdir temp_$projectname
cd temp_$projectname

pack [labelframe .questions.mastrtype -text "What type of input do you have for your MASTR?"] -fill x -pady 10 -side left
pack [radiobutton .questions.mastrtype.g -text gene -variable radiobuttonmastrtype -value gene] -side left
pack [radiobutton .questions.mastrtype.c -text coordinates -variable radiobuttonmastrtype -value coords] -side left
pack [radiobutton .questions.mastrtype.m -text miRNA -variable radiobuttonmastrtype -value miRNA] -side left
.questions.mastrtype.g select
pack [labelframe .questions.kit -text "Which MiSeq kit are you going to use?"] -fill x -pady 10 -side left
pack [radiobutton .questions.kit.2 -text "2" -variable radiobuttonkit -value 2] -side left
pack [radiobutton .questions.kit.3 -text "3" -variable radiobuttonkit -value 3] -side left
.questions.kit.3 select
pack [button .questions.confirm -text "Confirm" -command {
	set mastrtype $radiobuttonmastrtype
	set kit $radiobuttonkit
	if {$kit == "2"} {set maxsize 400}
	if {$kit == "3"} {set maxsize 500}
	destroy .questions.mastrtype
	destroy .questions.kit
	destroy .questions.confirm}] -side right
	tkwait variable mastrtype

if {$mastrtype == "coords"} {
	while {$varlist == ""} {
	pack [labelframe .questions.variants -text "Which are the variants for your MASTR design?"] -fill x -pady 10
	pack [button .questions.variants.getinputfile -text "Click to select a file" -command {set varlist [tk_getOpenFile -initialdir $dir]
	destroy .questions.variants
	}] -side left -fill x -padx 10
	pack [label .questions.variants.varinfo -text "Your file is tab separated and requires the fields 'chromosome begin end name'.\nNames need to be a unique value."] -side right
	tkwait window .questions.variants }
	lassign [handlevars $varlist] secretgene genelist

	pack [labelframe .questions.extend -text "Should the target around coding variants be extended to the exon in which they are localized?"] -pady 10
	pack [labelframe .questions.extend.q ] -fill x -side left -padx 10
	pack [radiobutton .questions.extend.q.f -text "full extension" -variable radiobuttonextend -value full] -side top
	pack [radiobutton .questions.extend.q.l -text "limited extension" -variable radiobuttonextend -value limited] -side top
	pack [radiobutton .questions.extend.q.n -text "no extension" -variable radiobuttonextend -value no -command {set targetsize 32}] -side top
	.questions.extend.q.l select
	pack [button .questions.extend.confirm -text "Confirm" -command {set extend $radiobuttonextend ; destroy .questions.extend}] -side left -padx 10
	pack [labelframe .questions.extend.info] -fill x -side right -padx 10
	pack [label .questions.extend.info.label -text "Full means possibly creating additional amplicons.\nLimited means not creating additional amplicons.\nNo extension means processing the variants like a standard SNP MASTR."]
	tkwait variable extend  }

if {$mastrtype == "gene"} {
	while {[llength $genelist] == 0} {
	set extend "full"
	pack [labelframe .questions.genes -text "Which gene targets should be included? (Provide HGNC names space separated and case sensitive.)"] -fill x -pady 10
	pack [entry .questions.genes.listinput -textvariable ingenelist] -fill x -side left
	pack [button .questions.genes.listconfirm -text "Confirm" -command {set genelist $ingenelist ; destroy .questions.genes}] -side left
	pack [button .questions.genes.getinputfile -text "Or select a file containing your genelist." -command {set genelist [tk_getOpenFile -initialdir $dir]
	destroy .questions.genes
	}] -side right -fill x
	tkwait variable genelist
	if {[file isfile $genelist]} {
	set fp4 [open "$genelist" r]
	set pregenelist [read $fp4]
	close $fp4
	set genelist [join $pregenelist]} } }

if {$mastrtype == "gene" || $secretgene == "true"} {
	pack [labelframe .questions.elements -text "Which gene elements should be included?"] -fill x -pady 10
	foreach i [list CDS UTR RNA intron upstream downstream] {
	pack [checkbutton .questions.elements.elem_$i -text "$i" -onvalue $i] -side left
	.questions.elements.elem_CDS select}
	pack [button .questions.elements.comfirm -text "Confirm" -command {
	foreach i [list $elem_CDS $elem_UTR $elem_RNA $elem_intron $elem_upstream $elem_downstream] {if {$i != "0"} {lappend content $i}}
	destroy .questions.elements}] -side right -padx 10
	pack [label .questions.elements.info -text "Upstream and downstream\nmeans extension of 2kb"] -side right -padx 10
	tkwait variable content

	while {[llength $dblist] == 0} {
	pack [labelframe .questions.db -text "Which databases should be used?"] -fill x -pady 10
	foreach i [list refGene ensGene knownGene gencode] {
	pack [checkbutton .questions.db.$i -text "$i" -onvalue $i] -side left
	.questions.db.$i select }
	pack [button .questions.db.confirm -text "Confirm" -command {
	foreach i [list $refGene $ensGene $knownGene $gencode] {if {$i != "0"} {lappend dblist $i}}
	destroy .questions.db }]
	tkwait variable dblist }}

if {$mastrtype == "miRNA"} {
	pack [labelframe .questions.mirna -text "Which mirna targets should be included?"] -fill x -pady 10
	pack [entry .questions.mirna.listinput -textvariable inlist] -fill x -side left -padx 10
	pack [button .questions.mirna.confirm -text "Confirm" -command {set mirnalist $inlist ; destroy .questions.mirna}] -side left -padx 10
	pack [label .questions.mirna.info -text "Provide mirbase name(s), case sensitive, space separated. (e.g. mir-1273d)"] -side right -padx 10
	tkwait variable mirnalist
	set mirnalist2 [join $mirnalist "\",\"hsa-"]
	set mirnalist3 \"hsa-$mirnalist2\" }

if {$extend ni [list no n] || $secretgene == "true"} {
while {$targetsize > $maxsize || $targetsize < 100} {
	pack [labelframe .questions.ts -text "How long (bp) do you want your target sequences to be maximally?"]
	pack [entry .questions.ts.input -textvariable ints] -fill x -side left -padx 10
	pack [button .questions.ts.confirm -text "Confirm" -command {set targetsize $ints ; destroy .questions.ts}] -side left -padx 10
	pack [label .questions.ts.info -text "Limited between 100bp and ${maxsize}bp."] -side right -padx 10
	tkwait window .questions.ts }}

while { $gccontent !="no" && ($gccontent > 100 || $gccontent <= 0) && $gcchecker != "on" } {
	pack [labelframe .questions.gc -text "On which gc percentage would you like to split your plexes?"]
	pack [entry .questions.gc.input -textvariable ingc] -fill x -side left -padx 10
	pack [button .questions.gc.confirm -text "Confirm" -command {set gccontent $ingc ; set gcsplit on ; destroy .questions.gc}] -side left -padx 10
	pack [button .questions.gc.gcchecker -text "Explore multiple gc values" -command {
			set gcchecker on
			set gcsplit on
			set gccontent [list 50 55 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80]
			destroy .questions.gc}] -side left -padx 10
	pack [button .questions.gc.nosp -text "No splitting" -command {set gcsplit off ; set gccontent no ; destroy .questions.gc}] -side left -padx 10
	tkwait window .questions.gc }

while { $snpfreq > 1 || $snpfreq < 0} {
	pack [labelframe .questions.fr -text "On which frequency would you like to mask your fasta files?"]
	pack [label .questions.fr.info -text "Provide decimal number e.g. 0.01"] -side right
	pack [entry .questions.fr.input -textvariable infr] -fill x -side left -padx 10
	pack [button .questions.fr.confirm -text "Confirm" -command {set snpfreq $infr ; destroy .questions.fr}] -side left -padx 10
	tkwait window .questions.fr }

pack [button .questions.buttonOK -text "Click to start the design." -command {destroy . ; set satisfied true}] -side top -padx 10
tkwait variable satisfied
puts "The script now starts the processing procedure for $projectname.\nMeanwhile, help yourself to a nice cup of coffee or tea."
}

### MAIN CODE ###
if { $mastrtype  in [list coords c]} {
	if {$secretgene == "true"} {
		gene2reg2content $dblist $content
		geneproc $dblist $genelist
		cg select -s - pre.processedtargets.tsv > secretgene_processedtargets.tsv
		set covered1 [cg covered secretgene_processedtargets.tsv]
		set covered2 [lindex $covered1 end]
		puts "Your targets have a total size of ${covered2}bp."
		file copy secretgene_processedtargets.tsv reg_gene.tsv }

## START OF NORMAL 'COORD MODE' PROCESSING
		cg select -s - cleanedfile.tsv sortedfile.tsv
		cg select -f {chromosome begin end name} sortedfile.tsv fieldfile.tsv
		exec sed -e {s/[[:space:]]*$//g} -e {/^[[:space:]]*$/d} sortedfile.tsv > trimmedsortedfile.tsv
		exec sed -e {s/chr//g} -e {s/omosome/chromosome/1} trimmedsortedfile.tsv > reg_name.tsv
		set targetnumber [expr [exec wc -l < reg_name.tsv] - 1]
		puts "You selected $targetnumber variants for your MASTR design."

	# Checks for uniqueness of names used.
		exec awk {!_[$4]++} reg_name.tsv > uniquename.tsv
		set wcuniquenames [expr [exec wc -l < uniquename.tsv] -1]
		if {$targetnumber != $wcuniquenames}	{puts "Error: You did not provide unique names!"; exit}

	# Processing of gene file as present on main server.
		if { $extend in [list full limited f l]} {
			if {![file isfile /complgen/bin/MastrDesign/db/compressed_reg_refGene.tsv.lz4]} {puts "Error: The refseq reference file does not exist or the path/name has been changed. Contact Wouter."; exit}
			cg select -q {$type == "CDS" || $type == "UTR" || $type == "RNA"} /complgen/bin/MastrDesign/db/compressed_reg_refGene.tsv.lz4 refgene_coding.tsv
			cg select -s - refgene_coding.tsv sortedrefgenecoding.tsv
			exec cg regjoin sortedrefgenecoding.tsv > joinedrefgenecoding.tsv
			cg regjoin reg_name.tsv joinedrefgenecoding.tsv > overlapfixed.tsv
			cg annotate overlapfixed.tsv annotatedexons.tsv reg_name.tsv
			cg select -q {$name != ""} annotatedexons.tsv codingtargets.tsv
			exec sed -e {s/chr//g} -e {s/omo/chromo/1} codingtargets.tsv > chrOkay_filteredcodingtargets.tsv

		# To include splice sites, the targets are extended in both directions with 10bp.
			cg select -f {chromosome {begin=$begin - 20} {end=$end + 20} name} chrOkay_filteredcodingtargets.tsv > extended_chrOkay_filteredcodingtargets.tsv
			cg cat -c 0 extended_chrOkay_filteredcodingtargets.tsv reg_name.tsv > catfile.tsv
			exec awk {!_[$4]++} catfile.tsv > uniquecat.tsv
			cg select -s - uniquecat.tsv processedtargets.tsv

		} elseif { $extend in [list no n]} {
			cg select -f {chromosome {begin=$begin - 15} {end=$end + 15} name} reg_name.tsv > processedtargets.tsv }

		if {$secretgene == "true"} {
			file rename processedtargets.tsv variants_processedtargets.tsv
			cg cat -c 0 variants_processedtargets.tsv secretgene_processedtargets.tsv > almost.processedtargets.tsv
			cg select -s - almost.processedtargets.tsv processedtargets.tsv}
			}
if {$mastrtype in [list gene g]} {
	gene2reg2content $dblist $content
	geneproc $dblist $genelist	}

if {$mastrtype in [list miRNA m]} {
	cg select -q "oneof(\$ID,$mirnalist3)" /complgen/refseq/hg19/reg_hg19_mirbase.tsv selectedmirna.tsv
	# To include splice sites, the targets are extended in both directions with 20bp.
	cg select -f {chromosome {begin=$begin - 20} {end=$end + 20} name=$ID} selectedmirna.tsv > extended_selectedmirna.tsv
	exec sed -e {s/chr//g} -e {s/omo/chromo/1} extended_selectedmirna.tsv > pre.processedtargets.tsv }

if 	{$mastrtype in [list miRNA m gene g]} {
	cg select -s - pre.processedtargets.tsv > processedtargets.tsv
	set covered1 [cg covered processedtargets.tsv]
	set covered2 [lindex $covered1 end]
	puts "Your targets have a total size of $covered2\ bp."}

# To specify the extension size in which primers may be located.
set Rarg [list processedtargets.tsv $targetsize 350]
exec Rscript /complgen/bin/MastrDesign/lib/exon_splitting.R $Rarg

# Now processing those exons which were too large for the user-specified target size, and thus splitted.
# The approach used is to annotated with reg_name.tsv to select the fragments containing the desired variant.
	if {$extend == "limited" || $extend == "l"} {
		exec paste first.tsv extended.tsv fas.target > pre.combinedtable.tsv
		cg select -s - pre.combinedtable.tsv > combinedtable.tsv
		file delete extended.tsv
		file delete fas.target
		file delete extendedfragments.tsv
		file delete extendedfragments2.tsv

		if {$secretgene == "true"} {
			cg annotate combinedtable.tsv annot_processedtargets.tsv reg_name.tsv reg_gene.tsv
			cg select -q {$name != "" || $gene != ""} annot_processedtargets.tsv limitedcombinedtable.tsv
		} else {
			cg annotate combinedtable.tsv annot_processedtargets.tsv reg_name.tsv
			cg select -q {$name != ""} annot_processedtargets.tsv limitedcombinedtable.tsv }

		exec cut -f1-5 limitedcombinedtable.tsv > first.tsv
		exec cut -f6-9 limitedcombinedtable.tsv > pre.extended.tsv
		exec sed {1 s/name2/name/} pre.extended.tsv > extended.tsv
		exec cut -f10-11 limitedcombinedtable.tsv > pre.fas.target
		exec tail -n +2 pre.fas.target > pre2.fas.target
		exec sed -e {s/\t/\n/g} pre2.fas.target > fas.target
	 }

if { $extend == "full" || $extend == "no" || $extend == "f" || $extend == "n"} {
	file rename extended.tsv pre.extended.tsv
	exec sed {1 s/name2/name/} pre.extended.tsv > extended.tsv
	exec tail -n +2 fas.target > pre2.fas.target
	exec sed -e {s/\t/\n/g} pre2.fas.target > fas.target }

exec awk {{print "chr"$1":"$2"-"$3}}  < extended.tsv > pre.positionlist.tsv
exec sed 1d pre.positionlist.tsv > positionlist.tsv

file rename fas.target $projectname.fas.target

# A custom UCSC browser track to make the results visible. The input for the track is a bed file construct.
cg select -f {chromosome begin end name1} first.tsv | tail -n +2 > pre.track.tsv

if {$mastrtype in [list coords c]} {exec bash /complgen/bin/MastrDesign/lib/trackisolatorc.sh}
if {$mastrtype  in [list gene g miRNA m]} {exec bash /complgen/bin/MastrDesign/lib/trackisolatorg.sh}
file rename Genome_Browser_Track.txt Genome_Browser_Track_$projectname.txt
file copy ./Genome_Browser_Track_$projectname.txt ../

# Optionally, the script can create the Custom Genome Browser Track using a Python script resulting in screenshots as output.
if {$screenshot == "True"} {
	set UCSCurl [prompt "You wanted screenshots of the Genome Browser Track. Now calling the Python script. Provide a valid public url to the Custom Genome Browser Track, for example using dropbox:"]
	exec /complgen/bin/anaconda/bin/python2.7 /complgen/bin/MastrDesign/lib/GenomeBrowserScreenshots.py $projectname $UCSCurl
	file copy screenshots$projectname ../ }

if {$gcchecker == "on"} {
foreach gcnumber $gccontent	{
	cg genome_seq -i "name" -g 0 -gs $gcnumber extended.tsv /complgen/refseq/hg19/ $projectname-$gcnumber.fas
	set highgc [expr [exec wc -l < $projectname-$gcnumber-highgc.fas] / 2]
	set lowgc [expr [exec wc -l < $projectname-$gcnumber-lowgc.fas] / 2]
	set freqhigh [expr round((100 * double($highgc)) / double($highgc + $lowgc))]
	set freqlow [expr round((100 * double($lowgc)) / double($highgc + $lowgc))]
	puts "Cut-off $gcnumber: $lowgc amplicons in LOW (${freqlow}%) and $highgc amplicons in HIGH (${freqhigh}%)"
	}
set chosencutoff [prompt "Which cut-off do you prefer? (Provide number)"]
set gccontent $chosencutoff }

	# The resulting target region files are processed to cg genome_seq, this produces a fasta_masked file and a fasta file of the targets.
	# The output files are copied to the main directory and the temporary directory is deleted.
	if {"e" ni $args} {
	if { $gcsplit == "off" } {
		cg genome_seq -n $snpfreq -r "N" -i "name" extended.tsv /complgen/refseq/hg19/ $projectname.fas.masked
		cg genome_seq -i "name" extended.tsv /complgen/refseq/hg19/ $projectname.fas
		cd ../
		file copy ./temp_$projectname/$projectname.fas.masked ./
		file copy ./temp_$projectname/$projectname.fas ./
		file copy ./temp_$projectname/$projectname.fas.target ./
	} elseif { $gcsplit =="on"} {
		cg genome_seq -n $snpfreq -r "N" -i "name" -g 0 -gs $gccontent extended.tsv /complgen/refseq/hg19/ pre.$projectname.fas.masked
		cg genome_seq -i "name" -g 0 -gs $gccontent extended.tsv /complgen/refseq/hg19/ pre.$projectname.fas

		exec awk {{print $1, $2}} pre.$projectname.fas-highgc.masked > $projectname-highgc.fas.masked
		exec awk {{print $1, $2}} pre.$projectname.fas-lowgc.masked > $projectname-lowgc.fas.masked
		exec awk {{print $1, $2}} pre.$projectname-highgc.fas > $projectname-highgc.fas
		exec awk {{print $1, $2}} pre.$projectname-lowgc.fas > $projectname-lowgc.fas

		set arg1high $projectname-highgc.fas
		set arg1low $projectname-lowgc.fas
		set arg2 $projectname.fas.target
		set arg3high $projectname-highgc.fas.target
		set arg3low $projectname-lowgc.fas.target

		if {[exec wc -l < $projectname-highgc.fas]  > 1 } {exec bash /complgen/bin/MastrDesign/lib/gcsplitter.sh $arg1high $arg2 $arg3high}
		if {[exec wc -l < $projectname-lowgc.fas]  > 1 } {exec bash /complgen/bin/MastrDesign/lib/gcsplitter.sh $arg1low $arg2 $arg3low}

		if {[exec wc -l < $projectname-highgc.fas] > 1 } {file copy $projectname-highgc.fas.masked ../ }
		if {[exec wc -l < $projectname-lowgc.fas] > 1 } {file copy $projectname-lowgc.fas.masked ../ }
		if {[exec wc -l < $projectname-highgc.fas] > 1 } {file copy $projectname-highgc.fas ../ }
		if {[exec wc -l < $projectname-lowgc.fas] > 1 } {file copy $projectname-lowgc.fas ../ }
		if {[exec wc -l < $projectname-highgc.fas] > 1 } {file copy $projectname-highgc.fas.target ../ }
		if {[exec wc -l < $projectname-lowgc.fas] > 1 } {file copy $projectname-lowgc.fas.target ../ }
		set highgc [expr [exec wc -l < $projectname-highgc.fas] / 2]
		set lowgc [expr [exec wc -l < $projectname-lowgc.fas] / 2]
		cd ../ }
		}
	if {"e" in $args} {
		puts "Not masking on repeats for Eline."
	if { $gcsplit == "off" } {
		cg genome_seq -n $snpfreq -i "name" extended.tsv /complgen/refseq/hg19/ $projectname.fas.masked
		cg genome_seq -i "name" extended.tsv /complgen/refseq/hg19/ $projectname.fas
		cd ../
		file copy ./temp_$projectname/$projectname.fas.masked ./
		file copy ./temp_$projectname/$projectname.fas ./
		file copy ./temp_$projectname/$projectname.fas.target ./
	} elseif { $gcsplit =="on"} {
		cg genome_seq -n $snpfreq -i "name" -g 0 -gs $gccontent extended.tsv /complgen/refseq/hg19/ pre.$projectname.fas.masked
		cg genome_seq -i "name" -g 0 -gs $gccontent extended.tsv /complgen/refseq/hg19/ pre.$projectname.fas

		exec awk {{print $1, $2}} pre.$projectname.fas-highgc.masked > $projectname-highgc.fas.masked
		exec awk {{print $1, $2}} pre.$projectname.fas-lowgc.masked > $projectname-lowgc.fas.masked
		exec awk {{print $1, $2}} pre.$projectname-highgc.fas > $projectname-highgc.fas
		exec awk {{print $1, $2}} pre.$projectname-lowgc.fas > $projectname-lowgc.fas

		set arg1high $projectname-highgc.fas
		set arg1low $projectname-lowgc.fas
		set arg2 $projectname.fas.target
		set arg3high $projectname-highgc.fas.target
		set arg3low $projectname-lowgc.fas.target

		if {[exec wc -l < $projectname-highgc.fas]  > 1 } {exec bash /complgen/bin/MastrDesign/lib/gcsplitter.sh $arg1high $arg2 $arg3high}
		if {[exec wc -l < $projectname-lowgc.fas]  > 1 } {exec bash /complgen/bin/MastrDesign/lib/gcsplitter.sh $arg1low $arg2 $arg3low}

		if {[exec wc -l < $projectname-highgc.fas] > 1 } {file copy $projectname-highgc.fas.masked ../ }
		if {[exec wc -l < $projectname-lowgc.fas] > 1 } {file copy $projectname-lowgc.fas.masked ../ }
		if {[exec wc -l < $projectname-highgc.fas] > 1 } {file copy $projectname-highgc.fas ../ }
		if {[exec wc -l < $projectname-lowgc.fas] > 1 } {file copy $projectname-lowgc.fas ../ }
		if {[exec wc -l < $projectname-highgc.fas] > 1 } {file copy $projectname-highgc.fas.target ../ }
		if {[exec wc -l < $projectname-lowgc.fas] > 1 } {file copy $projectname-lowgc.fas.target ../ }
		set highgc [expr [exec wc -l < $projectname-highgc.fas] / 2]
		set lowgc [expr [exec wc -l < $projectname-lowgc.fas] / 2]
		cd ../ }
		}
# Write log file
set time2 [clock format [clock seconds] -format %T]
if {$mastrtype == "c"} {set mastrtype "coords"}
if {$mastrtype == "g"} {set mastrtype "gene"}
if {$mastrtype == "gene"} {set extend "NA"}
if {$mastrtype == "m"} {set mastrtype "mirna"}
if {$extend == "f"} {set extend "full"}
if {$extend == "l"} {set extend "limited"}
if {$extend == "n"} {set extend "no"}

set logfilecontent "
------------------------------------------------------------------------------------------
Mastrdesign of $user on $server started at $time1 until $time2
Projectname chosen: $projectname
Specific modes: $args

Mastrtype $mastrtype
MiSeq kit V$kit

Targetfile $inputfile with $targetnumber variants
Smuggling $secretgene of $genelist

Genetargets are $genelist covering $covered2 bp
Chosen gene elements $content
Selected database(s) $dblist

miRNA-targets were $mirnalist3
Extension of coords to $extend coding sequence of exon
Maximum size of target $targetsize bp
Split fasta files on $gccontent % (or don't if 'no')
 Number of targets in highgc: $highgc
 Number of targets in lowgc: $lowgc
Mask fasta files on snp frequency of $snpfreq "

set logfile "$projectname.log"
set fileIdlogfile [open $logfile "w"]
puts -nonewline $fileIdlogfile $logfilecontent
close $fileIdlogfile

if {"t" ni $args} {file delete -force temp_$projectname}
puts "The design is finished and your files are ready to be used in mpcr.\nCheck the results using the Genome Browser and the primers using cg primercheck.\nGood luck with your optimization and have a nice day!"
