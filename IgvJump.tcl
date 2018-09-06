#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" ${1+"$@"}
# To examine the igv track for carriers of a certain variant based on a list and a compar file
# For sure works for MASTR, dunno about other project designs. If they follow the same 'project folder architecture', it should be fine.
# Script selects $number carriers 
#wdecoster

proc killlist listName {
    set tributenum [expr {int(rand()*[llength $listName])}]
	set list [lreplace $listName $tributenum $tributenum]
	return $list}
	
proc helpfunction {} {
	frame .main.files.helptext -padx 20 -pady 20
	destroy .main.files.help
	pack .main.files.helptext -side top
	pack [label .main.files.helptext.explanation -text "Welcome to this limited help function. For unexplained questions, ask Wouter.\nProvide an original compar file to the script containing your positions and samples of interest.\nOriginal means that you should use the link that was placed in your working directory initially.\nAlternatively, provide the original compar file as found in the finished project folder.\nYour list of variants of interest should be a tab separated file,\ncontaining the fields 'chromosome' 'begin' 'end' 'name' with a unique identifier in 'name'.\nSelect the number of carriers that maximally should be loaded in Igv.\nAt least one number should be > 0. When less are available, less are loaded.\nWhen more carriers are available,the script will randomly select the required number.\nPostions which are not present in the compar file are not shown.\nPlease inform me if it behaves unexpectedly. Feel free to suggest additional features." -justify left] }

proc tracelink compar {
	while {[file type $compar] == "link"} {set compar [file readlink $compar]}
	return $compar }

proc makearguments {target number projectpath time} {
	cg select -f "chromosome begin {bam-*= (\$zyg-gatk-crsbwa-*== \"$target\") ? concat(\"$projectpath\",\"*\",\"/\",\"map-crsbwa-*.bam\") : \"\"} name" ${time}_positionsofintrest.tsv > ${time}_samplenamereplace.tsv
	cg select -f "{p%%osition=concat(\"\$chromosome\",\":\",\"\$begin\",\"%%\")} bam-* {identifier=concat(\"%%@\",\"\$name\")}" ${time}_samplenamereplace.tsv > ${time}_BaseCompar.tsv
	set fp [open "${time}_BaseCompar.tsv" r]
	set basecompar [read $fp]
	close $fp
	set linebyline [split $basecompar "\n"]

	set arguments [dict create]
	set identifierdict [dict create]
	foreach line $linebyline {
		set trimmed [regsub -all {@} [regsub -all {%%@} [regsub -all {\s} $line {}] { }] {,}]
		if {![string match *%%* $trimmed]} {
			dict set arguments [lindex $trimmed 0] [lindex $trimmed 1]
			dict set identifierdict [lindex $trimmed 0] [lindex $trimmed 2]} }
	dict for {location bam} $arguments {
		set bamlist [split $bam ,]
		while {[llength $bamlist] > $number} {
			dict unset arguments $location
			set bamlist [killlist $bamlist]}
		set bam [join $bamlist ,]
		dict set arguments $location $bam }
	return [list $arguments $identifierdict]}

lassign {"" "" 0 0 0 0}	compar positions numberm numbert numberu numberr
package require Tk

pack [frame .main -padx 10 -pady 10] -side top -fill both
pack [frame .main.files] -side top
wm protocol . WM_DELETE_WINDOW { if {[tk_messageBox -message "Are you sure you want to quit?" -type yesno] eq "yes"} { exit } }
wm title . "Examination of variants in igv"
pack [labelframe .main.files.frame -width 55 -bd 0] -fill x
pack [button .main.files.frame.inputfileMASTR -text "Click to select the original compar file" -command {
	set compar [tk_getOpenFile -initialdir [pwd] -title "Select compar file"]
	set projectpath @/
	append projectpath [join [lrange [split [tracelink $compar] /] 1 end-2] /]
	append projectpath /
	.main.files.frame.inputfilePOS configure -state active	} -width 78 -overrelief raised]
pack [button .main.files.frame.inputfilePOS -text "Click to select the list containing variants of interest" -state disabled -command {
	set positions [tk_getOpenFile -initialdir [pwd] -title "Select variant file"]
	.main.files.buttons.start configure -state active}  -width 78 -overrelief raised ]

while {[expr $numberm + $numbert + $numberr] == 0} {
pack [labelframe .main.files.carriers -width 55 -bd 0] -fill x
pack [labelframe .main.files.carriers.frame2 -width 55 -bd 0] -fill x
pack [label .main.files.carriers.frame2.exp -text "Select the maximum number of homozygous carriers for each position to be included in igv:"] -side left -padx 8
pack [ttk::combobox .main.files.carriers.frame2.num -values [list 0 1 2 3 4 5 6 7 8 9 10] -justify center -textvariable numberm -state readonly -width 3] -side right -padx 10

pack [labelframe .main.files.carriers.frame3 -width 55 -bd 0] -fill x
pack [label .main.files.carriers.frame3.exp -text "Select the maximum number of heterozygous carriers for each position to be included in igv:"] -side left -padx 8
pack [ttk::combobox .main.files.carriers.frame3.num -values [list 0 1 2 3 4 5 6 7 8 9 10] -justify center -textvariable numbert -state readonly -width 3] -side right -padx 10

pack [labelframe .main.files.carriers.frame4 -width 55 -bd 0] -fill x
pack [label .main.files.carriers.frame4.exp -text "Select the maximum number of non-carriers for each position to be included in igv:"] -side left -padx 8
pack [ttk::combobox .main.files.carriers.frame4.num -values [list 0 1 2 3 4 5 6 7 8 9 10] -justify center -textvariable numberr -state readonly -width 3] -side right -padx 10

pack [labelframe .main.files.carriers.frame5 -width 55 -bd 0] -fill x
pack [label .main.files.carriers.frame5.exp -text "Select the maximum number of unsequenced samples for each position to be included in igv:"] -side left -padx 8
pack [ttk::combobox .main.files.carriers.frame5.num -values [list 0 1 2 3 4 5 6 7 8 9 10] -justify center -textvariable numberu -state readonly -width 3] -side right -padx 10

pack [labelframe .main.files.buttons -bd 0]
pack [button .main.files.buttons.start -width 67 -state disabled -overrelief raised -text "Click to start" -command {set satisfied true}] -side left
if {$positions != {}} {.main.files.buttons.start configure -state active}
pack [button .main.files.buttons.help -text "Help" -command {helpfunction} -width 8] -side right
tkwait variable satisfied 
if {[expr $numberm + $numbert + $numberr + $numberu] == 0} {tk_messageBox -message "Invalid number of selected individuals: minimally take 1 of one category."
	destroy .main.files.carriers 
	destroy .main.files.buttons}}
destroy .main.files
puts "Preparing arguments for igv.\nDepending on size of variantlist and compar file, this might take a few minutes."

set samplelist [split [cg select -h $compar]]
set blancolist [list]
foreach i $samplelist { if {[string match -nocase *blanco* $i] } { lappend blancolist $i} }
set blancofields [join $blancolist]
set time [clock format [clock seconds] -format %m%d_%H%M]
cg select -rf "$blancofields" $compar > ${time}_noblancs.tsv
cg select -s - $positions > reg_name.tsv
cg annotate ${time}_noblancs.tsv ${time}_positionannotated.tsv reg_name.tsv
cg select -q {$name != ""} ${time}_positionannotated.tsv > ${time}_positionsofintrest.tsv

set keylist [list]
set idlist [list]
foreach i {m t r u} { 
	set arg number$i
	if {[set $arg] > 0} {
		lassign [makearguments $i [set $arg] $projectpath $time] arguments${i} identifierdict${i}
		lappend idlist [set identifierdict${i}]
		lappend keylist {*}[dict keys [set arguments${i}]]} }
set identifierdict [dict merge {*}$idlist]
file delete ${time}_noblancs.tsv ${time}_sortedpositions.tsv ${time}_positionannotated.tsv ${time}_positionsofintrest.tsv ${time}_samplenamereplace.tsv ${time}_BaseCompar.tsv reg_name.tsv

tk_messageBox -message "Preparation finished!"
frame .main.next -width 100 -height 60
pack .main.next -side top -fill both
wm title . "Jumping through the bam files."
pack [button .main.next.exit -text "Click to exit" -command {exit} -width 15] -side right  -fill x

foreach loca [lsort -unique $keylist] {
	if {$loca != ""} {
		set identifier [dict get $identifierdict $loca]
		set bam ""
		array set argray [list m argumentsm t argumentst r argumentsr u argumentsu]
		foreach i {m t r u} {
			set arg number$i
			if {[set $arg] > 0} {
				if {[dict exists [set $argray($i)] $loca] == 1} {
					set bam2 [dict get [set $argray($i)] $loca]
					if {[string length $bam] > 0} {append bam ","}
					set bam $bam$bam2}}
				}
		pack [button .main.next.pos -text "Click to load next position: $identifier" -width [expr [string length $identifier] + 20] -command {
			set next $loca
			destroy .main.next.pos
			destroy .main.next.id
			exec /complgen/bin/IGV_2.3.50/igv.sh $bam $loca &}]  -fill x -side left
		tkwait variable next } }